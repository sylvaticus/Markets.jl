using Test
using Markets
using DataFrames

const JuMP = Markets.JuMP

# helpers to read single values out of the result tables
price(res, r, p) = only(res.prices[(res.prices.region .== r) .& (res.prices.product .== p), :price])
getq(df, r, p)   = (rows = df[(df.region .== r) .& (df.product .== p), :quantity]; isempty(rows) ? 0.0 : only(rows))
solved(res)      = JuMP.termination_status(res.model) in (JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.OPTIMAL)

@testset "Markets.jl" begin

@testset "Constructors" begin
    n = leontief(:a, 2)
    @test n.products == [:a] && n.composite == 2.0
    c = ces(1.5, [:a, :b], [3, 1]; sigma = 2)
    @test c.shares ≈ [0.75, 0.25]
    @test c.sigma == 2.0 && c.phi == 1.0
    @test_throws AssertionError ces(1.0, [:a, :b], [1.0]; sigma = 2)
    @test_throws AssertionError ces(1.0, [:a, :b], [1.0, 1.0]; sigma = 0.5)
    p = Process(:mill, [leontief(:a, 1)], [:b => 1, :c => 0.5]; vacost = 10)
    @test p.outputs == [:b => 1.0, :c => 0.5] && p.regions === :all
end

@testset "Spatial price equilibrium (one product, two regions)" begin
    τ = 5.0
    d = MarketData(
        regions   = [:A, :B],
        products  = [:x],
        demand    = [DemandSpec(:x, :A; p0 = 100, q0 = 10, elasticity = 1.5),
                     DemandSpec(:x, :B; p0 = 100, q0 = 10, elasticity = 1.5)],
        supply    = [SupplySpec(:x, :A; p0 = 50, q0 = 30, elasticity = 1.0),   # cheap
                     SupplySpec(:x, :B; p0 = 150, q0 = 5, elasticity = 1.0)],  # expensive
        processes = Process[],
        tradable  = [:x],
        transport = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ),
    )
    res = solve_market(d)
    @test solved(res)
    # A exports to B, and the price gap equals the transport cost
    @test only(res.trade.from) == :A && only(res.trade.to) == :B
    @test price(res, :B, :x) - price(res, :A, :x) ≈ τ rtol = 1e-4
    # prices equal willingness to pay (inverse demand) and marginal supply cost
    for s in d.demand
        D = getq(res.consumption, s.region, :x)
        @test price(res, s.region, :x) ≈ s.p0 * (D / s.q0)^(-1 / s.elasticity) rtol = 1e-4
    end
    S_A = getq(res.production, :A, :x)
    @test price(res, :A, :x) ≈ 50 * (S_A / 30)^(1 / 1.0) rtol = 1e-4
    # material balance in each region
    for r in (:A, :B)
        nt = res.net_trade[res.net_trade.region .== r, :net]
        @test getq(res.production, r, :x) - getq(res.consumption, r, :x) ≈ only(nt) atol = 1e-5
    end
    # for_region slices consistently
    a = for_region(res, :A)
    @test all(a.prices.region .== :A) && nrow(a.trade) == 1
end

@testset "Leontief and CES processes" begin
    σ, δ = 2.0, [0.6, 0.4]
    d = MarketData(
        regions   = [:R],
        products  = [:in1, :in2, :mid, :out],
        demand    = [DemandSpec(:out, :R; p0 = 500, q0 = 10, elasticity = 1.5)],
        supply    = [SupplySpec(:in1, :R; p0 = 50, q0 = 10, elasticity = 0.8),
                     SupplySpec(:in2, :R; p0 = 50, q0 = 10, elasticity = 0.8)],
        processes = [Process(:stage1, [ces(2.0, [:in1, :in2], δ; sigma = σ)], [:mid => 1.0]; vacost = 20),
                     Process(:stage2, [leontief(:mid, 1.25)], [:out => 1.0]; vacost = 30)],
        tradable  = Symbol[],
        transport = Dict{Tuple{Symbol,Symbol,Symbol},Float64}(),
    )
    res = solve_market(d)
    @test solved(res)
    p(x) = price(res, :R, x)
    # Leontief: output price = coeff * input price + value-added cost
    @test p(:out) ≈ 1.25 * p(:mid) + 30 rtol = 1e-4
    # CES cost-minimisation FOC: q1/q2 = (δ1/δ2) * (p2/p1)^σ
    q1, q2 = getq(res.production, :R, :in1), getq(res.production, :R, :in2)
    @test q1 / q2 ≈ (δ[1] / δ[2]) * (p(:in2) / p(:in1))^σ rtol = 1e-3
    # CES zero profit: mid price = composite unit cost + value-added cost
    unitcost = 2.0 * (δ[1] * p(:in1)^(1 - σ) + δ[2] * p(:in2)^(1 - σ))^(1 / (1 - σ))
    @test p(:mid) ≈ unitcost + 20 rtol = 1e-3
end

@testset "Forest example" begin
    ex = Module(:ForestExample)
    Base.include(ex, joinpath(@__DIR__, "..", "examples", "forest", "example_data.jl"))
    d   = ex.example_market
    res = solve_market(d)
    @test solved(res)
    @test nrow(res.prices) == length(d.regions) * length(d.products)
    @test all(res.prices.price .> 0)
    # papermill is a Leontief process: price(paper) = 1.1 price(pulp) + 200
    for r in d.regions
        @test price(res, r, :paper) ≈ 1.1 * price(res, r, :pulp) + 200 rtol = 1e-4
    end
    # no arbitrage: price gaps never exceed transport costs, and match them on used routes
    for ((p, from, to), τ) in d.transport
        @test price(res, to, p) - price(res, from, p) <= τ + 1e-3
    end
    for row in eachrow(res.trade)
        row.quantity > 1e-3 || continue
        @test price(res, row.to, row.product) - price(res, row.from, row.product) ≈
              d.transport[(row.product, row.from, row.to)] rtol = 1e-3
    end
    # world trade balances: total exports = total imports for each product
    for g in groupby(res.net_trade, :product)
        @test sum(g.net) ≈ 0 atol = 1e-5
    end
end

end
