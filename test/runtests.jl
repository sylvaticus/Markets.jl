using Test
using Markets
using DataFrames

const JuMP = Markets.JuMP

# helpers to read single values out of the result tables
price(res, r, p) = only(res.prices[(res.prices.region .== r) .& (res.prices.product .== p), :price])
pprice(res, r, p) = only(res.prices[(res.prices.region .== r) .& (res.prices.product .== p), :producer_price])
getq(df, r, p)   = (rows = df[(df.region .== r) .& (df.product .== p), :quantity]; isempty(rows) ? 0.0 : only(rows))
flow(res, p, from, to) = (rows = res.trade[(res.trade.product .== p) .& (res.trade.from .== from) .&
                                           (res.trade.to .== to), :quantity]; isempty(rows) ? 0.0 : only(rows))
solved(res)      = JuMP.termination_status(res.model) in (JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.OPTIMAL)

@testset "Markets.jl" begin

@testset "Constructors" begin
    n = leontief(product = :a, coeff = 2)
    @test n.products == [:a] && n.composite == 2.0 && n.shares == [1.0]
    c = ces(composite = 1.5, products = [:a, :b], shares = [3, 1], sigma = 2)
    @test c.shares ≈ [0.75, 0.25]
    @test c.composite == 1.5 && c.sigma == 2.0 && c.phi == 1.0
    # shares default to equal ones, sigma must be > 1, shares must align
    @test ces(composite = 1.0, products = [:a, :b, :c], sigma = 2).shares ≈ [1/3, 1/3, 1/3]
    @test_throws AssertionError ces(composite = 1.0, products = [:a, :b], shares = [1.0], sigma = 2)
    @test_throws AssertionError ces(composite = 1.0, products = [:a, :b], shares = [1.0, 1.0], sigma = 0.5)
    # every field is a keyword, and the optional ones have defaults
    nest = Nest(composite = 2.0, products = [:a, :b], sigma = 3)
    @test nest.shares ≈ [0.5, 0.5] && nest.phi == 1.0
    @test_throws UndefKeywordError Nest(composite = 1.0, products = [:a])
    p = Process(name = :mill, inputs = [leontief(product = :a, coeff = 1)],
                outputs = [:b => 1, :c => 0.5], vacost = 10)
    @test p.outputs == [:b => 1.0, :c => 0.5] && p.regions === :all && p.vacost == 10.0
    @test_throws UndefKeywordError Process(name = :mill, inputs = Nest[])
    d = MarketData(regions = [:A], products = [:x])
    @test isempty(d.demand) && isempty(d.supply) && isempty(d.processes) &&
          isempty(d.tradable) && isempty(d.transport)
    @test_throws UndefKeywordError DemandSpec(product = :x, region = :A, p0 = 1, q0 = 1)
end

@testset "Spatial price equilibrium (one product, two regions)" begin
    τ = 5.0
    d = MarketData(
        regions   = [:A, :B],
        products  = [:x],
        demand    = [DemandSpec(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 1.5),
                     DemandSpec(product = :x, region = :B, p0 = 100, q0 = 10, elasticity = 1.5)],
        supply    = [SupplySpec(product = :x, region = :A, p0 = 50, q0 = 30, elasticity = 1.0),   # cheap
                     SupplySpec(product = :x, region = :B, p0 = 150, q0 = 5, elasticity = 1.0)],  # expensive
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
        demand    = [DemandSpec(product = :out, region = :R, p0 = 500, q0 = 10, elasticity = 1.5)],
        supply    = [SupplySpec(product = :in1, region = :R, p0 = 50, q0 = 10, elasticity = 0.8),
                     SupplySpec(product = :in2, region = :R, p0 = 50, q0 = 10, elasticity = 0.8)],
        processes = [Process(name = :stage1, vacost = 20,
                             inputs  = [ces(composite = 2.0, products = [:in1, :in2], shares = δ, sigma = σ)],
                             outputs = [:mid => 1.0]),
                     Process(name = :stage2, vacost = 30,
                             inputs  = [leontief(product = :mid, coeff = 1.25)],
                             outputs = [:out => 1.0])],
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

@testset "Armington trade" begin
    τ = 5.0
    # A is the cheap producer, B the expensive one
    economy(arm) = MarketData(
        regions   = [:A, :B],
        products  = [:x],
        demand    = [DemandSpec(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 1.5),
                     DemandSpec(product = :x, region = :B, p0 = 100, q0 = 10, elasticity = 1.5)],
        supply    = [SupplySpec(product = :x, region = :A, p0 = 50,  q0 = 30, elasticity = 1.0),
                     SupplySpec(product = :x, region = :B, p0 = 150, q0 = 5,  elasticity = 1.0)],
        tradable  = [:x],
        transport = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ),
        armington = arm)

    hom = solve_market(economy(Armington[]))

    @testset "σ = Inf is the homogeneous model" begin
        inf = solve_market(economy([Armington(product = :x, sigma = Inf)]))
        @test solved(inf)
        for r in (:A, :B)
            @test price(inf, r, :x) ≈ price(hom, r, :x) rtol = 1e-8
            @test price(inf, r, :x) ≈ pprice(inf, r, :x)          # no wedge
            @test getq(inf.consumption, r, :x) ≈ getq(hom.consumption, r, :x) rtol = 1e-8
        end
    end

    @testset "σ → ∞ converges to the spatial equilibrium" begin
        err(σ) = begin
            r = solve_market(economy([Armington(product = :x, sigma = σ)]))
            @test solved(r)
            maximum(abs(price(r, g, :x) - price(hom, g, :x)) / price(hom, g, :x) for g in (:A, :B))
        end
        e10, e50, e200, e1000 = err(10.0), err(50.0), err(200.0), err(1000.0)
        @test e10 > e50 > e200 > e1000        # monotone convergence
        @test e1000 < 1e-3
    end

    @testset "imperfect substitution: cross-hauling and price wedges" begin
        σ = 3.0
        res = solve_market(economy([Armington(product = :x, sigma = σ)]))
        @test solved(res)
        # both directions are traded at once, which perfect substitutes never do
        @test flow(res, :x, :A, :B) > 1e-3
        @test flow(res, :x, :B, :A) > 1e-3
        @test nrow(hom.trade) == 1                     # ... unlike the homogeneous case
        # producer prices are no longer tied together by the transport cost
        @test pprice(res, :B, :x) - pprice(res, :A, :x) > τ
        # users pay the CES price index of the delivered prices of both origins
        for r in (:A, :B), o in (:A, :B)
            q(o, r) = pprice(res, o, :x) + (o == r ? 0.0 : τ)
            @test price(res, r, :x) ≈
                  sum(0.5 * q(o, r)^(1 - σ) for o in (:A, :B))^(1 / (1 - σ)) rtol = 1e-4
        end
    end

    @testset "calibrated shares drive the mix" begin
        σ, δ = 4.0, Dict((:A, :A) => 0.8, (:B, :A) => 0.2,   # region A is home-biased
                         (:A, :B) => 0.3, (:B, :B) => 0.7)
        res = solve_market(economy([Armington(product = :x, sigma = σ, shares = δ)]))
        @test solved(res)
        # CES demand: the quantity ratio follows the share and price ratios
        dom, imp = getq(res.production, :A, :x) - flow(res, :x, :A, :B), flow(res, :x, :B, :A)
        qd, qi = pprice(res, :A, :x), pprice(res, :B, :x) + τ
        @test dom / imp ≈ (δ[(:A, :A)] / δ[(:B, :A)]) * (qi / qd)^σ rtol = 1e-3
        # and the composite price is the share-weighted CES index
        @test price(res, :A, :x) ≈
              (δ[(:A, :A)] * qd^(1 - σ) + δ[(:B, :A)] * qi^(1 - σ))^(1 / (1 - σ)) rtol = 1e-4
    end

    @testset "an origin with no share is left out" begin
        res = solve_market(economy([Armington(product = :x, sigma = 4,
                                              shares = Dict((:A, :A) => 1.0,     # A: domestic only
                                                            (:A, :B) => 0.5, (:B, :B) => 0.5))]))
        @test solved(res)
        @test flow(res, :x, :B, :A) == 0.0                    # no route B → A is used
        @test flow(res, :x, :A, :B) > 1e-3
        @test price(res, :A, :x) ≈ pprice(res, :A, :x) rtol = 1e-6   # single variety: no wedge
    end

    @testset "prices link across regions through the elasticity" begin
        σ = 4.0
        spec = [Armington(product = :x, sigma = σ)]
        base  = solve_market(economy(spec))
        # a 30% outward shift of region A's supply curve
        shocked = solve_market(MarketData(
            regions   = [:A, :B],
            products  = [:x],
            demand    = economy(spec).demand,
            supply    = [SupplySpec(product = :x, region = :A, p0 = 50, q0 = 39, elasticity = 1.0),
                         SupplySpec(product = :x, region = :B, p0 = 150, q0 = 5, elasticity = 1.0)],
            tradable  = [:x],
            transport = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ),
            armington = spec))
        @test solved(shocked)
        # it reaches region B: cheaper there too, and B buys more from A
        @test price(shocked, :B, :x) < price(base, :B, :x)
        @test flow(shocked, :x, :A, :B) > flow(base, :x, :A, :B)
        # B's own producers are squeezed by the competing variety
        @test pprice(shocked, :B, :x) < pprice(base, :B, :x)
    end

    @testset "multi-stage chain with an Armington input" begin
        d = MarketData(
            regions   = [:A, :B],
            products  = [:w, :f],
            demand    = [DemandSpec(product = :f, region = :A, p0 = 400, q0 = 10, elasticity = 1.4),
                         DemandSpec(product = :f, region = :B, p0 = 400, q0 = 10, elasticity = 1.4)],
            supply    = [SupplySpec(product = :w, region = :A, p0 = 60, q0 = 30, elasticity = 0.7),
                         SupplySpec(product = :w, region = :B, p0 = 90, q0 = 15, elasticity = 0.7)],
            processes = [Process(name = :mill, vacost = 50,
                                 inputs  = [leontief(product = :w, coeff = 2.0)],
                                 outputs = [:f => 1.0])],
            tradable  = [:w],
            transport = Dict((:w, :A, :B) => 8.0, (:w, :B, :A) => 8.0),
            armington = [Armington(product = :w, sigma = 5)])
        res = solve_market(d)
        @test solved(res)
        # the mill pays the composite price of its input and gets the producer
        # price of its output, so zero profit links the two
        for r in (:A, :B)
            @test pprice(res, r, :f) ≈ 2.0 * price(res, r, :w) + 50 rtol = 1e-4
        end
        # :f is homogeneous and untradable, so its two prices coincide
        @test price(res, :A, :f) ≈ pprice(res, :A, :f) rtol = 1e-8
    end

    @testset "specification errors" begin
        @test_throws ArgumentError Armington(product = :x, sigma = 0.5)
        @test_throws ArgumentError Armington(product = :x, sigma = 1.0)
        @test_throws ArgumentError Armington(product = :x, sigma = 2,
                                             shares = Dict((:A, :A) => -1.0))
        @test_throws ArgumentError solve_market(economy([Armington(product = :nope, sigma = 2)]))
        @test_throws ArgumentError solve_market(economy([Armington(product = :x, sigma = 2),
                                                         Armington(product = :x, sigma = 3)]))
        # no origin left with a positive share for destination A
        @test_throws ArgumentError solve_market(economy([Armington(product = :x, sigma = 2,
                                                   shares = Dict((:A, :B) => 1.0))]))
    end
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
