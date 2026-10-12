using Test
using Markets
using DataFrames
using JuMP: termination_status, objective_value, all_variables, name
import JuMP

# helpers to read single values out of the result tables
price(res, r, p) = only(res.prices[(res.prices.region .== r) .& (res.prices.product .== p), :price])
pprice(res, r, p) = only(res.prices[(res.prices.region .== r) .& (res.prices.product .== p), :producer_price])
getq(df, r, p)   = (rows = df[(df.region .== r) .& (df.product .== p), :quantity]; isempty(rows) ? 0.0 : only(rows))
flow(res, p, from, to) = (rows = res.trade[(res.trade.product .== p) .& (res.trade.from .== from) .&
                                           (res.trade.to .== to), :quantity]; isempty(rows) ? 0.0 : only(rows))
solved(res)      = res.status in (Markets.MOI.LOCALLY_SOLVED, Markets.MOI.OPTIMAL)

@testset "Markets.jl" begin

@testset "Constructors" begin
    n = leontief(product = :a, coeff = 2)
    @test n.products == [:a] && n.composite == 2.0 && n.shares == [1.0]
    c = ces(composite = 1.5, products = [:a, :b], shares = [3, 1], sigma = 2)
    @test c.shares ≈ [0.75, 0.25]
    @test c.composite == 1.5 && c.sigma == 2.0 && c.phi == 1.0
    # shares default to equal ones, sigma must be > 1, shares must align
    @test ces(composite = 1.0, products = [:a, :b, :c], sigma = 2).shares ≈ [1/3, 1/3, 1/3]
    @test_throws ArgumentError ces(composite = 1.0, products = [:a, :b], shares = [1.0], sigma = 2)
    @test_throws ArgumentError ces(composite = 1.0, products = [:a, :b], shares = [1.0, 1.0], sigma = 0.5)
    # a Nest built directly validates like one built by `ces`
    @test_throws ArgumentError Nest(1.0, [:a, :b], [1.0], 2.0, 1.0)      # shares mismatch
    @test_throws ArgumentError Nest(1.0, [:a, :b], [1.0, 1.0], 0.5, 1.0) # sigma too low
    @test_throws ArgumentError Nest(1.0, [:a, :a], [1.0, 1.0], 2.0, 1.0) # repeated product
    @test_throws ArgumentError Nest(0.0, [:a], [1.0], 0.0, 1.0)          # nothing needed
    @test leontief(product = :a, coeff = 2).sigma == 0.0                 # still fine
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
          isempty(d.transport_costs)
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
        transport_costs = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ),
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

@testset "A product may be supplied, consumed and used as an input at once" begin
    # :logs come out of the forest, and are either burnt as they are or sawn
    d = MarketData(
        regions   = [:A], products = [:logs, :board],
        demand    = [DemandSpec(product = :logs,  region = :A, p0 = 40,  q0 = 20, elasticity = 0.4),
                     DemandSpec(product = :board, region = :A, p0 = 300, q0 = 10, elasticity = 1.2)],
        supply    = [SupplySpec(product = :logs, region = :A, p0 = 40, q0 = 35, elasticity = 0.6)],
        processes = [Process(name = :mill, vacost = 120,
                             inputs  = [leontief(product = :logs, coeff = 2.0)],
                             outputs = [:board => 1.0])])
    res = solve_market(d)
    @test solved(res)
    burnt  = getq(res.consumption, :A, :logs)
    milled = 2 * only(res.activity.level)
    # one price clears both uses, and the harvest covers exactly the two of them
    @test getq(res.production, :A, :logs) ≈ burnt + milled rtol = 1e-6
    @test burnt > 1e-3 && milled > 1e-3        # both uses are actually served
    @test price(res, :A, :board) ≈ 2 * price(res, :A, :logs) + 120 rtol = 1e-6
    # the direct use competes with the mill: a higher willingness to pay for
    # firewood takes wood away from the board
    hotter = solve_market(MarketData(regions = d.regions, products = d.products,
        demand = [DemandSpec(product = :logs, region = :A, p0 = 70, q0 = 20, elasticity = 0.4),
                  d.demand[2]],
        supply = d.supply, processes = d.processes))
    @test getq(hotter.consumption, :A, :logs) > burnt
    @test 2 * only(hotter.activity.level) < milled
    @test price(hotter, :A, :logs) > price(res, :A, :logs)
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

@testset "Demand elasticities of any sign of (η-1)" begin
    economy(η) = MarketData(
        regions  = [:A, :B], products = [:x],
        demand   = [DemandSpec(product = :x, region = r, p0 = 100, q0 = 10, elasticity = η)
                    for r in (:A, :B)],
        supply   = [SupplySpec(product = :x, region = :A, p0 = 50,  q0 = 30, elasticity = 1.0),
                    SupplySpec(product = :x, region = :B, p0 = 150, q0 = 5,  elasticity = 1.0)],
        transport_costs = Dict((:x, :A, :B) => 5.0, (:x, :B, :A) => 5.0))

    # η ≤ 1 is not a restriction: the benefit term is an antiderivative of the
    # inverse demand, which stays concave and leaves the equilibrium defined
    for η in (0.3, 0.6, 1.0, 1.5, 3.0)
        res = solve_market(economy(η))
        @test solved(res)
        for r in (:A, :B)
            D = getq(res.consumption, r, :x)
            @test price(res, r, :x) ≈ 100 * (D / 10)^(-1 / η) rtol = 1e-4
        end
    end
    # η = 1 is the logarithmic case, and the switch to it is continuous
    just_above = solve_market(economy(1 + 1e-7))
    unit       = solve_market(economy(1.0))
    @test price(just_above, :A, :x) ≈ price(unit, :A, :x) rtol = 1e-5
    # less elastic demand means a supply shock moves prices more
    shock(η) = begin
        base = solve_market(economy(η))
        more = MarketData(regions = [:A, :B], products = [:x],
                 demand = economy(η).demand,
                 supply = [SupplySpec(product = :x, region = :A, p0 = 50, q0 = 45, elasticity = 1.0),
                           SupplySpec(product = :x, region = :B, p0 = 150, q0 = 5, elasticity = 1.0)],
                 transport_costs = Dict((:x, :A, :B) => 5.0, (:x, :B, :A) => 5.0))
        1 - price(solve_market(more), :A, :x) / price(base, :A, :x)
    end
    @test shock(0.4) > shock(1.5) > 0

    @test_throws ArgumentError DemandSpec(product = :x, region = :A, p0 = 100, q0 = 10,
                                          elasticity = 0.0)
    @test_throws ArgumentError DemandSpec(product = :x, region = :A, p0 = -1, q0 = 10,
                                          elasticity = 1.5)
    @test_throws ArgumentError SupplySpec(product = :x, region = :A, p0 = 50, q0 = 0,
                                          elasticity = 1.0)
    @test_throws ArgumentError SupplySpec(product = :x, region = :A, p0 = 50, q0 = 10,
                                          elasticity = -0.5)
end

@testset "Exogenous drivers of the curves" begin
    economy(; supply_kw = (), demand_kw = ()) = MarketData(
        regions = [:A], products = [:x],
        demand  = [DemandSpec(product = :x, region = :A, p0 = 100, q0 = 10,
                              elasticity = 1.5; demand_kw...)],
        supply  = [SupplySpec(product = :x, region = :A, p0 = 50, q0 = 20,
                              elasticity = 0.8; supply_kw...)])
    on_curve(s, S) = s.p0 * (S / (s.shift * s.q0))^(1 / s.elasticity)

    @testset "a shifter moves the curve, not its elasticity" begin
        base = solve_market(economy())
        @test solved(base)
        @test price(base, :A, :x) ≈ on_curve(economy().supply[1],
                                             getq(base.production, :A, :x)) rtol = 1e-6
        for m in (0.6, 1.5, 3.0)
            res = solve_market(economy(supply_kw = (shift = m,)))
            @test solved(res)
            s = economy(supply_kw = (shift = m,)).supply[1]
            S = getq(res.production, :A, :x)
            # the solution still sits on the shifted curve ...
            @test price(res, :A, :x) ≈ on_curve(s, S) rtol = 1e-6
            # ... and at any given price it offers `m` times the base quantity
            @test m * base.production.quantity[1] * (price(res, :A, :x) /
                  price(base, :A, :x))^s.elasticity ≈ S rtol = 1e-6
        end
        # more supply is cheaper, less is dearer
        @test price(solve_market(economy(supply_kw = (shift = 1.5,))), :A, :x) <
              price(base, :A, :x) <
              price(solve_market(economy(supply_kw = (shift = 0.6,))), :A, :x)
        # a demand shifter works the same way, in the other direction
        @test price(solve_market(economy(demand_kw = (shift = 1.4,))), :A, :x) >
              price(base, :A, :x)
        # shift = 1 is the calibrated curve, i.e. the default
        @test getq(solve_market(economy(supply_kw = (shift = 1.0,))).production, :A, :x) ≈
              getq(base.production, :A, :x) rtol = 1e-9
    end

    @testset "a ceiling truncates the curve and earns a rent" begin
        free = solve_market(economy())
        cap  = 0.6 * getq(free.production, :A, :x)
        res  = solve_market(economy(supply_kw = (capacity = cap,)))
        @test solved(res)
        S = getq(res.production, :A, :x)
        @test S ≈ cap rtol = 1e-6                       # the bound binds
        @test price(res, :A, :x) > price(free, :A, :x)  # scarcity raises the price
        # the price now exceeds the marginal cost on the curve: that gap is the
        # scarcity rent of the resource
        @test price(res, :A, :x) > on_curve(economy(supply_kw = (capacity = cap,)).supply[1], S)
        # a ceiling above what would be produced anyway changes nothing
        loose = solve_market(economy(supply_kw = (capacity = 10 * cap,)))
        @test getq(loose.production, :A, :x) ≈ getq(free.production, :A, :x) rtol = 1e-6
    end

    @testset "shifters from drivers" begin
        @test exogenous_shift(stock = (level = 118.0, reference = 100.0, elasticity = 0.6),
                              roads = (level = 1.05, reference = 1.0, elasticity = 0.25)) ≈
              1.18^0.6 * 1.05^0.25
        # a driver at its reference level leaves the curve alone
        @test exogenous_shift(stock = (level = 7.0, reference = 7.0, elasticity = 0.9)) == 1
        # an elasticity of zero means the driver does not matter
        @test exogenous_shift(x = (level = 3.0, reference = 1.0, elasticity = 0.0)) == 1
        @test_throws ArgumentError exogenous_shift()
        @test_throws ArgumentError exogenous_shift(a = (level = 1.0, reference = 1.0))
        @test_throws ArgumentError exogenous_shift(a = (level = -1.0, reference = 1.0,
                                                        elasticity = 1.0))
        @test_throws ArgumentError SupplySpec(product = :x, region = :A, p0 = 1, q0 = 1,
                                              elasticity = 1, shift = 0)
        @test_throws ArgumentError SupplySpec(product = :x, region = :A, p0 = 1, q0 = 1,
                                              elasticity = 1, capacity = -1)
        @test_throws ArgumentError DemandSpec(product = :x, region = :A, p0 = 1, q0 = 1,
                                              elasticity = 2, shift = -0.5)
    end
end

@testset "Validation: the economy must say what it means" begin
    ## a small, correct economy that each case below breaks in one way
    good(; kw...) = MarketData(;
        regions = [:A, :B], products = [:logs, :board],
        demand  = [DemandSpec(product = :board, region = :A, p0 = 300, q0 = 10,
                              elasticity = 1.2)],
        supply  = [SupplySpec(product = :logs, region = :A, p0 = 40, q0 = 35,
                              elasticity = 0.6)],
        processes = [Process(name = :mill, vacost = 120,
                             inputs = [leontief(product = :logs, coeff = 2.0)],
                             outputs = [:board => 1.0])],
        kw...)
    @test validate(good()) === nothing
    @test solved(solve_market(good()))

    @testset "undeclared names" begin
        # a mistyped input used to make the mill produce out of nothing
        @test_throws ArgumentError solve_market(good(processes =
            [Process(name = :mill, vacost = 120,
                     inputs = [leontief(product = :logz, coeff = 2.0)],
                     outputs = [:board => 1.0])]))
        @test_throws ArgumentError solve_market(good(demand =
            [DemandSpec(product = :nope, region = :A, p0 = 1, q0 = 1, elasticity = 2)]))
        @test_throws ArgumentError solve_market(good(supply =
            [SupplySpec(product = :logs, region = :Mars, p0 = 1, q0 = 1, elasticity = 1)]))
        @test_throws ArgumentError solve_market(good(processes =
            [Process(name = :mill, vacost = 1, regions = [:Mars],
                     inputs = [leontief(product = :logs, coeff = 1)],
                     outputs = [:board => 1.0])]))
        @test_throws ArgumentError solve_market(good(processes =
            [Process(name = :mill, vacost = 1,
                     inputs = [leontief(product = :logs, coeff = 1)],
                     outputs = [:nope => 1.0])]))
        @test_throws ArgumentError solve_market(good(
            transport_costs = Dict((:board, :A, :Mars) => 1.0)))
        @test_throws ArgumentError solve_market(good(
            transport_costs = Dict((:nope, :A, :B) => 1.0)))
        @test_throws ArgumentError solve_market(good(armington =
            [Armington(product = :nope, sigma = 3)]))
    end

    @testset "duplicates" begin
        @test_throws ArgumentError solve_market(good(regions = [:A, :A]))
        @test_throws ArgumentError solve_market(good(products = [:logs, :logs, :board]))
        @test_throws ArgumentError solve_market(good(processes =
            [Process(name = :mill, vacost = 1, inputs = [leontief(product = :logs, coeff = 1)],
                     outputs = [:board => 1.0]),
             Process(name = :mill, vacost = 9, inputs = [leontief(product = :logs, coeff = 1)],
                     outputs = [:board => 1.0])]))
        @test_throws ArgumentError solve_market(good(supply =
            [SupplySpec(product = :logs, region = :A, p0 = 40, q0 = 35, elasticity = 0.6),
             SupplySpec(product = :logs, region = :A, p0 = 90, q0 = 35, elasticity = 0.6)]))
        @test_throws ArgumentError solve_market(good(demand =
            [DemandSpec(product = :board, region = :A, p0 = 300, q0 = 10, elasticity = 1.2),
             DemandSpec(product = :board, region = :A, p0 = 300, q0 = 10, elasticity = 1.2)]))
    end

    @testset "routes that go nowhere" begin
        @test_throws ArgumentError solve_market(good(
            transport_costs = Dict((:logs, :A, :A) => 5.0)))   # a route to itself
        # a route is all it takes to make a product tradable
        @test solved(solve_market(good(
            transport_costs = Dict((:logs, :A, :B) => 5.0, (:logs, :B, :A) => 5.0))))
    end

    @testset "selectors that cannot bite" begin
        # a tax on trade ignores `region`, so setting one would tax every route
        @test_throws ArgumentError solve_market(good(
            taxes = [Tax(on = :trade, region = :A, rate = 10.0)]))
        @test_throws ArgumentError solve_market(good(
            taxes = [Tax(on = :process, product = :logs, rate = 10.0)]))
        @test_throws ArgumentError solve_market(good(
            taxes = [Tax(on = :supply, product = :nope, rate = 10.0)]))
        @test_throws ArgumentError solve_market(good(
            taxes = [Tax(on = :process, process = :nomill, rate = 10.0)]))
        @test_throws ArgumentError solve_market(good(trade_limits =
            [TradeLimit(name = :l, products = :nope, capacity = 1.0)]))
        @test_throws ArgumentError solve_market(good(trade_limits =
            [TradeLimit(name = :l, from = :Mars, capacity = 1.0)]))
    end

    @testset "economies with no finite answer" begin
        # a process with no inputs and no cost is a free source of its output
        @test_throws ArgumentError solve_market(MarketData(regions = [:A], products = [:g],
            demand = [DemandSpec(product = :g, region = :A, p0 = 300, q0 = 10,
                                 elasticity = 1.2)],
            processes = [Process(name = :magic, vacost = 0.0, inputs = Nest[],
                                 outputs = [:g => 1.0])]))
        # the same process is legitimate once it costs something ...
        @test solved(solve_market(MarketData(regions = [:A], products = [:g],
            demand = [DemandSpec(product = :g, region = :A, p0 = 300, q0 = 10,
                                 elasticity = 1.2)],
            processes = [Process(name = :backstop, vacost = 50.0, inputs = Nest[],
                                 outputs = [:g => 1.0])])))
        # ... or once it is bounded
        @test solved(solve_market(MarketData(regions = [:A], products = [:g],
            demand = [DemandSpec(product = :g, region = :A, p0 = 300, q0 = 10,
                                 elasticity = 1.2)],
            processes = [Process(name = :windfall, vacost = 0.0, inputs = Nest[],
                                 capacity = 4.0, outputs = [:g => 1.0])])))
        @test_throws ArgumentError solve_market(MarketData(regions = Symbol[],
            products = [:g]))
        @test_throws ArgumentError solve_market(MarketData(regions = [:A],
            products = Symbol[]))
    end

    @testset "a solver verdict is not a suggestion" begin
        # demand with no source anywhere is infeasible; it used to be returned
        # as LOCALLY_SOLVED because the violation sat under Ipopt's tolerance
        starved = MarketData(regions = [:A], products = [:x],
            demand = [DemandSpec(product = :x, region = :A, p0 = 100, q0 = 10,
                                 elasticity = 1.5)])
        @test_throws ErrorException solve_market(starved)
        quiet = solve_market(starved; on_failure = :ignore)
        @test !solved(quiet)                       # the status says so plainly
        @test quiet.status == Markets.MOI.LOCALLY_INFEASIBLE
        @test_logs (:warn,) solve_market(starved; on_failure = :warn)
        @test_throws ArgumentError solve_market(starved; on_failure = :nonsense)
    end
end

@testset "Bounds on transformation and on trade" begin
    economy(; cap = Inf, limits = TradeLimit[]) = MarketData(
        regions   = [:A, :B], products = [:w, :g],
        demand    = [DemandSpec(product = :g, region = r, p0 = 300, q0 = 10, elasticity = 1.2)
                     for r in (:A, :B)],
        supply    = [SupplySpec(product = :w, region = :A, p0 = 40, q0 = 40, elasticity = 0.8),
                     SupplySpec(product = :w, region = :B, p0 = 90, q0 = 5,  elasticity = 0.8)],
        processes = [Process(name = :mill, vacost = 50, capacity = cap,
                             inputs  = [leontief(product = :w, coeff = 2.0)],
                             outputs = [:g => 1.0])],
        transport_costs = Dict((p, o, r) => (p === :w ? 6.0 : 10.0)
                         for p in (:w, :g), o in (:A, :B), r in (:A, :B) if o != r),
        trade_limits = limits)
    level(res, r) = only(res.activity[res.activity.region .== r, :level])

    @testset "a capacity holds a process back" begin
        free = solve_market(economy())
        @test solved(free)
        @test all(isinf, free.activity.capacity) && all(iszero, free.activity.rent)
        cap = 0.6 * level(free, :A)
        res = solve_market(economy(cap = cap))
        @test solved(res)
        @test level(res, :A) ≈ cap rtol = 1e-5          # the bound binds in A ...
        @test level(res, :B) > level(free, :B)          # ... and work moves to B
        @test only(res.activity[res.activity.region .== :A, :rent]) > 0   # a quasi-rent
        # a per-region dictionary bounds only the regions it names
        split = solve_market(economy(cap = Dict(:A => cap)))
        @test level(split, :A) ≈ cap rtol = 1e-5
        @test isinf(only(split.activity[split.activity.region .== :B, :capacity]))
        @test_throws ArgumentError Process(name = :p, inputs = Nest[], outputs = [:g => 1.0],
                                           capacity = 0.0)
    end

    @testset "a trade ceiling binds and earns a congestion rent" begin
        free  = solve_market(economy())
        moved = only(free.trade[free.trade.product .== :g, :quantity])
        quota = 0.55 * moved
        res = solve_market(economy(limits = [TradeLimit(name = :quota, products = :g,
                                             from = :A, to = :B, capacity = quota)]))
        @test solved(res)
        @test flow(res, :g, :A, :B) ≈ quota rtol = 1e-4
        row = only(eachrow(res.limits))
        @test row.limit === :quota
        @test row.used ≈ quota rtol = 1e-4
        @test row.rent > 0
        # with the ceiling binding, the price gap is transport *plus* the rent
        @test price(res, :B, :g) - price(res, :A, :g) ≈ 10.0 + row.rent rtol = 1e-4
    end

    @testset "a ceiling over several products uses its weights" begin
        res = solve_market(economy(limits = [TradeLimit(name = :port, products = [:w, :g],
                                             from = :A, capacity = 6.0,
                                             weights = Dict(:w => 0.5, :g => 1.0))]))
        @test solved(res)
        used = 0.5 * flow(res, :w, :A, :B) + flow(res, :g, :A, :B)
        @test only(res.limits.used) ≈ used rtol = 1e-5
        @test used ≈ 6.0 rtol = 1e-4
        # a loose ceiling changes nothing
        loose = solve_market(economy(limits = [TradeLimit(name = :port, from = :A,
                                                          capacity = 1e4)]))
        @test flow(loose, :g, :A, :B) ≈ flow(solve_market(economy()), :g, :A, :B) rtol = 1e-5
        @test only(loose.limits.rent) < 1e-6
        @test_throws ArgumentError TradeLimit(capacity = -1.0)
    end
end

@testset "Substitution in final demand" begin
    nest = DemandNest(sigma = 4.0, products = [:pellets, :firewood], shares = [0.6, 0.4])
    economy(sh) = MarketData(
        regions = [:A], products = [:pellets, :firewood],
        demand  = [DemandSpec(product = nest, region = :A, p0 = 100, q0 = 20,
                              elasticity = 0.5)],
        supply  = [SupplySpec(product = :pellets,  region = :A, p0 = 100, q0 = 12,
                              elasticity = 0.8),
                   SupplySpec(product = :firewood, region = :A, p0 = 100, q0 = 8,
                              elasticity = 0.8, shift = sh)])

    @testset "the group splits by its shares at equal prices" begin
        res = solve_market(economy(1.0))
        @test solved(res)
        @test price(res, :A, :pellets) ≈ price(res, :A, :firewood) rtol = 1e-6
        @test getq(res.consumption, :A, :pellets) / getq(res.consumption, :A, :firewood) ≈
              0.6 / 0.4 rtol = 1e-5
    end

    @testset "a price change moves demand between the substitutes" begin
        base, tight = solve_market(economy(1.0)), solve_market(economy(0.5))
        @test solved(tight)
        # firewood gets scarcer and dearer, so consumers move to pellets
        @test price(tight, :A, :firewood) > price(base, :A, :firewood)
        @test getq(tight.consumption, :A, :firewood) < getq(base.consumption, :A, :firewood)
        @test getq(tight.consumption, :A, :pellets)  > getq(base.consumption, :A, :pellets)
        # and the CES demand condition holds exactly
        for res in (base, tight)
            @test getq(res.consumption, :A, :pellets) / getq(res.consumption, :A, :firewood) ≈
                  (0.6 / 0.4) * (price(res, :A, :firewood) /
                                 price(res, :A, :pellets))^4 rtol = 1e-5
        end
    end

    @testset "a nest of one product is an ordinary curve" begin
        plain  = MarketData(regions = [:A], products = [:x],
                   demand = [DemandSpec(product = :x, region = :A, p0 = 100, q0 = 10,
                                        elasticity = 1.2)],
                   supply = [SupplySpec(product = :x, region = :A, p0 = 50, q0 = 20,
                                        elasticity = 0.8)])
        nested = MarketData(regions = plain.regions, products = plain.products,
                   demand = [DemandSpec(product = DemandNest(sigma = 3.0, products = [:x]),
                                        region = :A, p0 = 100, q0 = 10, elasticity = 1.2)],
                   supply = plain.supply)
        @test getq(solve_market(nested).consumption, :A, :x) ≈
              getq(solve_market(plain).consumption, :A, :x) rtol = 1e-6
    end

    @testset "nesting rules" begin
        @test_throws ArgumentError DemandNest(sigma = 0.5, products = [:a])
        @test_throws ArgumentError DemandNest(sigma = 2.0)                       # empty
        @test_throws ArgumentError DemandNest(sigma = 2.0, products = [:a, :b],
                                              shares = [1.0])                     # mismatch
        @test_throws ArgumentError DemandNest(sigma = 4.0, products = [:a],
                   nests = [DemandNest(sigma = 2.0, products = [:b])])             # loose inside
        # a product may not have two demand curves in one region
        twice = MarketData(regions = [:A], products = [:x],
            demand = [DemandSpec(product = :x, region = :A, p0 = 1, q0 = 1, elasticity = 2),
                      DemandSpec(product = :x, region = :A, p0 = 1, q0 = 1, elasticity = 2)],
            supply = [SupplySpec(product = :x, region = :A, p0 = 1, q0 = 1, elasticity = 1)])
        @test_throws ArgumentError solve_market(twice)
    end
end

@testset "Taxes and subsidies" begin
    economy(taxes) = MarketData(
        regions  = [:A, :B], products = [:x],
        demand   = [DemandSpec(product = :x, region = r, p0 = 100, q0 = 10, elasticity = 1.2)
                    for r in (:A, :B)],
        supply   = [SupplySpec(product = :x, region = :A, p0 = 50,  q0 = 30, elasticity = 0.8),
                    SupplySpec(product = :x, region = :B, p0 = 150, q0 = 5,  elasticity = 0.8)],
        transport_costs = Dict((:x, :A, :B) => 5.0, (:x, :B, :A) => 5.0),
        taxes    = taxes)
    base = solve_market(economy(Tax[]))

    @testset "a tax on production is a wedge under the price" begin
        res = solve_market(economy([Tax(on = :supply, region = :A, rate = 6.0)]))
        @test solved(res)
        s = economy(Tax[]).supply[1]
        q = getq(res.production, :A, :x)
        # the market price covers the marginal cost plus the tax
        @test price(res, :A, :x) ≈ s.p0 * (q / s.q0)^(1 / s.elasticity) + 6.0 rtol = 1e-5
        @test price(res, :A, :x) > price(base, :A, :x)
        @test only(res.taxes.revenue) ≈ 6.0 * q rtol = 1e-5
        # the relative form multiplies the curve instead
        rel = solve_market(economy([Tax(on = :supply, region = :A, rate = 0.2, relative = true)]))
        qr  = getq(rel.production, :A, :x)
        @test price(rel, :A, :x) ≈ 1.2 * s.p0 * (qr / s.q0)^(1 / s.elasticity) rtol = 1e-5
    end

    @testset "a tax on consumption is a wedge above the price" begin
        res = solve_market(economy([Tax(on = :demand, region = :A, rate = 0.25,
                                        relative = true)]))
        @test solved(res)
        q = getq(res.consumption, :A, :x)
        # consumers pay the curve; the market price is that divided by 1 + rate
        @test 1.25 * price(res, :A, :x) ≈ 100 * (q / 10)^(-1 / 1.2) rtol = 1e-5
        @test q < getq(base.consumption, :A, :x)
        @test only(res.taxes.revenue) ≈ 0.25 * price(res, :A, :x) * q rtol = 1e-4
    end

    @testset "a tariff, specific and ad valorem" begin
        spec = solve_market(economy([Tax(on = :trade, from = :A, to = :B, rate = 10.0)]))
        @test solved(spec)
        @test price(spec, :B, :x) - price(spec, :A, :x) ≈ 15.0 rtol = 1e-5  # 5 freight + 10
        @test flow(spec, :x, :A, :B) < flow(base, :x, :A, :B)
        # an ad valorem rate is a fraction of the price at the origin, found by
        # iteration; at the fixed point the gap must be freight + rate × price
        adv = solve_market(economy([Tax(on = :trade, from = :A, to = :B, rate = 0.2,
                                        relative = true)]))
        @test solved(adv)
        @test price(adv, :B, :x) - price(adv, :A, :x) ≈
              5.0 + 0.2 * price(adv, :A, :x) rtol = 1e-5
        @test only(adv.taxes.revenue) ≈ 0.2 * price(adv, :A, :x) *
              flow(adv, :x, :A, :B) rtol = 1e-4
    end

    @testset "a subsidy is a negative tax" begin
        res = solve_market(economy([Tax(on = :supply, region = :B, rate = -20.0)]))
        @test solved(res)
        @test getq(res.production, :B, :x) > getq(base.production, :B, :x)
        @test only(res.taxes.revenue) < 0
    end

    @testset "taxing a process" begin
        d = MarketData(regions = [:A], products = [:w, :g],
            demand    = [DemandSpec(product = :g, region = :A, p0 = 300, q0 = 10,
                                    elasticity = 1.2)],
            supply    = [SupplySpec(product = :w, region = :A, p0 = 40, q0 = 40,
                                    elasticity = 0.8)],
            processes = [Process(name = :mill, vacost = 50,
                                 inputs = [leontief(product = :w, coeff = 2.0)],
                                 outputs = [:g => 1.0])])
        plain = solve_market(d)
        taxed = solve_market(MarketData(regions = d.regions, products = d.products,
            demand = d.demand, supply = d.supply, processes = d.processes,
            taxes = [Tax(name = :levy, on = :process, process = :mill, rate = 0.3,
                         relative = true)]))
        @test solved(taxed)
        # the levy raises the conversion cost by 30%, so the margin widens by 15
        @test price(taxed, :A, :g) - 2 * price(taxed, :A, :w) ≈ 65.0 rtol = 1e-5
        @test price(plain, :A, :g) - 2 * price(plain, :A, :w) ≈ 50.0 rtol = 1e-5
        @test only(taxed.taxes.revenue) ≈ 15.0 * only(taxed.activity.level) rtol = 1e-5
    end

    @testset "specification errors" begin
        @test_throws ArgumentError Tax(on = :nonsense, rate = 1.0)
        @test_throws ArgumentError Tax(on = :demand, rate = -1.5, relative = true)
        # a group of substitutes shares one curve, so it must share its tax
        nest = DemandNest(sigma = 3.0, products = [:a, :b])
        d = MarketData(regions = [:A], products = [:a, :b],
            demand = [DemandSpec(product = nest, region = :A, p0 = 10, q0 = 2, elasticity = 2)],
            supply = [SupplySpec(product = p, region = :A, p0 = 10, q0 = 1, elasticity = 1)
                      for p in (:a, :b)],
            taxes  = [Tax(on = :demand, product = :a, rate = 0.2, relative = true)])
        @test_throws ArgumentError solve_market(d)
    end
end

@testset "Prices keep their sign" begin
    # :junk is an unavoidable by-product nobody wants; the balance is an
    # equality, so somebody must take it and pay to burn it. Its equilibrium
    # price is therefore negative, and must be reported as such.
    d = MarketData(
        regions   = [:A], products = [:w, :g, :junk],
        demand    = [DemandSpec(product = :g, region = :A, p0 = 300, q0 = 10, elasticity = 1.5)],
        supply    = [SupplySpec(product = :w, region = :A, p0 = 50, q0 = 20, elasticity = 1.0)],
        processes = [Process(name = :mill, vacost = 20,
                             inputs  = [leontief(product = :w, coeff = 1.0)],
                             outputs = [:g => 1.0, :junk => 0.5]),
                     Process(name = :incinerator, vacost = 15,
                             inputs  = [leontief(product = :junk, coeff = 1.0)],
                             outputs = Pair{Symbol,Float64}[])])
    res = solve_market(d)
    @test solved(res)
    @test price(res, :A, :junk) ≈ -15 rtol = 1e-6     # the cost of burning it
    @test price(res, :A, :w) > 0 && price(res, :A, :g) > 0
    # zero-cost disposal is how a user asks for free disposal, and then the
    # price of the by-product is zero rather than negative
    free = solve_market(MarketData(regions = d.regions, products = d.products,
             demand = d.demand, supply = d.supply,
             processes = [d.processes[1],
                          Process(name = :dump, vacost = 0.0,
                                  inputs  = [leontief(product = :junk, coeff = 1.0)],
                                  outputs = Pair{Symbol,Float64}[])]))
    @test solved(free)
    @test abs(price(free, :A, :junk)) < 1e-6
end

@testset "Numerical robustness" begin
    @testset "the quantity floor is set per variable, not per economy" begin
        # one obscure product with a tiny reference quantity must not drag the
        # floor of every other variable down with it
        d = MarketData(regions = [:A], products = [:x, :dust],
            demand = [DemandSpec(product = :x,    region = :A, p0 = 100, q0 = 50,
                                 elasticity = 1.3),
                      DemandSpec(product = :dust, region = :A, p0 = 5, q0 = 0.001,
                                 elasticity = 1.3)],
            supply = [SupplySpec(product = :x,    region = :A, p0 = 40, q0 = 120,
                                 elasticity = 0.8),
                      SupplySpec(product = :dust, region = :A, p0 = 5, q0 = 0.001,
                                 elasticity = 0.8)])
        res = solve_market(d)
        @test solved(res)
        floors = Dict(name(v) => JuMP.lower_bound(v) for v in all_variables(res.model)
                      if JuMP.has_lower_bound(v))
        @test floors["D[A,x]"] > 1e3 * floors["D[A,dust]"]   # each to its own scale
        @test floors["D[A,x]"] < 1e-3 * 50                   # and negligible for it
        # an explicit qfloor overrides the rule with one absolute floor
        flat = solve_market(d; qfloor = 1e-7)
        @test all(≈(1e-7), [JuMP.lower_bound(v) for v in all_variables(flat.model)
                            if JuMP.has_lower_bound(v) && JuMP.lower_bound(v) > 0])
    end

    @testset "an ad valorem tariff settles even when it is enormous" begin
        # homogeneous trade, so that the price gap on a used route is exactly
        # freight plus the wedge and the fixed point can be read off it
        economy(rate, eta) = MarketData(regions = [:A, :B], products = [:x],
            demand = [DemandSpec(product = :x, region = r, p0 = 100, q0 = 10,
                                 elasticity = eta) for r in (:A, :B)],
            supply = [SupplySpec(product = :x, region = :A, p0 = 30,  q0 = 40,
                                 elasticity = 3.0),
                      SupplySpec(product = :x, region = :B, p0 = 200, q0 = 2,
                                 elasticity = 3.0)],
            transport_costs = Dict((:x, :A, :B) => 1.0, (:x, :B, :A) => 1.0),
            taxes = [Tax(on = :trade, product = :x, from = :A, to = :B, rate = rate,
                         relative = true)])
        for rate in (0.2, 0.5, 1.0, 2.0, 4.0)        # up to a 400% tariff
            res = solve_market(economy(rate, 2.0))
            @test solved(res)
            @test flow(res, :x, :A, :B) > 1e-3            # the route is still used
            # at the fixed point the wedge is the rate times the origin price
            @test price(res, :B, :x) - price(res, :A, :x) ≈
                  1.0 + rate * price(res, :A, :x) rtol = 1e-5
        end
        # beyond that the tariff is prohibitive: the route empties, and the gap
        # is then bounded by the wedge rather than equal to it
        choked = solve_market(economy(40.0, 2.0))
        @test solved(choked)
        @test flow(choked, :x, :A, :B) <= 1e-3
        @test price(choked, :B, :x) - price(choked, :A, :x) <=
              1.0 + 40.0 * price(choked, :A, :x) + 1e-6
    end
end

@testset "Disposal of an unwanted by-product" begin
    economy(extra) = MarketData(regions = [:A], products = [:w, :g, :junk],
        demand = [DemandSpec(product = :g, region = :A, p0 = 300, q0 = 10,
                             elasticity = 1.5)],
        supply = [SupplySpec(product = :w, region = :A, p0 = 50, q0 = 20,
                             elasticity = 1.0)],
        processes = [Process(name = :mill, vacost = 20,
                             inputs = [leontief(product = :w, coeff = 1.0)],
                             outputs = [:g => 1.0, :junk => 0.5]); extra])
    # with only a costly route out, the by-product is a bad worth minus that cost
    burn = solve_market(economy([Process(name = :burn, vacost = 15,
                                 inputs = [leontief(product = :junk, coeff = 1.0)],
                                 outputs = Pair{Symbol,Float64}[])]))
    @test price(burn, :A, :junk) ≈ -15 rtol = 1e-6
    # `disposal` is the same thing, written once
    free = solve_market(economy([disposal(product = :junk)]))
    @test abs(price(free, :A, :junk)) < 1e-6              # free disposal: price 0
    paid = solve_market(economy([disposal(product = :junk, cost = 12)]))
    @test price(paid, :A, :junk) ≈ -12 rtol = 1e-6
    @test :dispose_junk in solve_market(economy([disposal(product = :junk)])).activity.process
    # and the good products are unaffected in sign
    @test price(free, :A, :g) > 0 && price(free, :A, :w) > 0
end

@testset "Armington share calibration" begin
    flows  = Dict((:EU, :EU) => 80.0,  (:NA, :EU) => 15.0,  (:AS, :EU) => 5.0,
                  (:NA, :NA) => 70.0,  (:EU, :NA) => 20.0,  (:AS, :NA) => 10.0)
    prices = Dict((:EU, :EU) => 210.0, (:NA, :EU) => 235.0, (:AS, :EU) => 250.0,
                  (:NA, :NA) => 200.0, (:EU, :NA) => 240.0, (:AS, :NA) => 245.0)
    σ = 4
    δ = armington_shares(; flows, prices, sigma = σ)
    # shares sum to one per destination ...
    for dest in (:EU, :NA)
        @test sum(v for (k, v) in δ if k[2] == dest) ≈ 1
    end
    # ... and reproduce the observed sourcing through the CES demand condition
    for dest in (:EU, :NA)
        keys_d  = [k for k in keys(δ) if k[2] == dest]
        implied = Dict(k => δ[k] * prices[k]^(-σ) for k in keys_d)
        tot, ftot = sum(values(implied)), sum(flows[k] for k in keys_d)
        for k in keys_d
            @test implied[k] / tot ≈ flows[k] / ftot rtol = 1e-8
        end
    end
    # they are not the observed shares themselves, which is the point
    @test !isapprox(δ[(:EU, :EU)], 0.8; rtol = 1e-3)
    @test_throws ArgumentError armington_shares(; flows, prices, sigma = 0.5)
    @test_throws ArgumentError armington_shares(flows = Dict((:EU, :EU) => 1.0),
                                                prices = Dict{Tuple{Symbol,Symbol},Float64}(),
                                                sigma = 3)
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
        transport_costs = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ),
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
            transport_costs = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ),
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
            transport_costs = Dict((:w, :A, :B) => 8.0, (:w, :B, :A) => 8.0),
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

@testset "Origin groups (nested substitution)" begin
    regions   = [:EU, :NA, :AS]
    transport_costs = Dict((:x, o, r) => 6.0 for o in regions, r in regions if o != r)
    shares    = Dict((:EU, :EU) => 0.5,  (:NA, :EU) => 0.3,  (:AS, :EU) => 0.2,
                     (:NA, :NA) => 0.6,  (:EU, :NA) => 0.25, (:AS, :NA) => 0.15,
                     (:AS, :AS) => 0.55, (:EU, :AS) => 0.2,  (:NA, :AS) => 0.25)
    economy(arm) = MarketData(
        regions  = regions, products = [:x],
        demand   = [DemandSpec(product = :x, region = r, p0 = 100, q0 = 10, elasticity = 1.5)
                    for r in regions],
        supply   = [SupplySpec(product = :x, region = :EU, p0 = 60, q0 = 20, elasticity = 1.0),
                    SupplySpec(product = :x, region = :NA, p0 = 55, q0 = 25, elasticity = 1.0),
                    SupplySpec(product = :x, region = :AS, p0 = 70, q0 = 15, elasticity = 1.0)],
        transport_costs = transport_costs, armington = arm)

    @testset "a group as substitutable as its parent changes nothing" begin
        flat   = solve_market(economy([Armington(product = :x, sigma = 3, shares = shares)]))
        redund = solve_market(economy([Armington(product = :x, sigma = 3, shares = shares,
                                       nests = [OriginNest(sigma = 3, origins = [:EU, :NA])])]))
        @test solved(redund)
        for r in regions
            @test price(redund, r, :x) ≈ price(flat, r, :x) rtol = 1e-5
            @test pprice(redund, r, :x) ≈ pprice(flat, r, :x) rtol = 1e-5
        end
    end

    @testset "the price index is the two-level CES one" begin
        σt, σg = 2.0, 9.0
        res = solve_market(economy([Armington(product = :x, sigma = σt, shares = shares,
                                    nests = [OriginNest(sigma = σg, origins = [:EU, :NA])])]))
        @test solved(res)
        for r in regions
            q(o)  = pprice(res, o, :x) + (o == r ? 0.0 : 6.0)
            raw   = Dict(o => shares[(o, r)] for o in regions)
            group = raw[:EU] + raw[:NA]
            # price of the EU/NA bundle, then of the composite over bundle and AS
            Pg = sum(raw[o] / group * q(o)^(1 - σg) for o in (:EU, :NA))^(1 / (1 - σg))
            P  = (group * Pg^(1 - σt) + raw[:AS] * q(:AS)^(1 - σt))^(1 / (1 - σt))
            @test price(res, r, :x) ≈ P rtol = 1e-6
        end
    end

    @testset "a tightening group converges to a single pooled good" begin
        gap(σg) = begin
            res = solve_market(economy([Armington(product = :x, sigma = 2.0, shares = shares,
                                        nests = [OriginNest(sigma = σg, origins = [:EU, :NA])])]))
            @test solved(res)
            abs(pprice(res, :EU, :x) - pprice(res, :NA, :x))
        end
        @test gap(3.0) > gap(20.0) > gap(200.0)      # the two prices are pulled together
    end

    @testset "structures and elasticities may differ by destination" begin
        res = solve_market(economy([
            # the EU market keeps AS at arm's length but treats EU and NA alike
            Armington(product = :x, destination = :EU, sigma = 1.5, shares = shares,
                      nests = [OriginNest(sigma = 12, origins = [:EU, :NA])]),
            # the others substitute freely between all three
            Armington(product = :x, destination = [:NA, :AS], sigma = 8, shares = shares)]))
        @test solved(res)
        # the EU price index uses σ = 1.5 over {EU,NA bundle, AS} ...
        q(o, r) = pprice(res, o, :x) + (o == r ? 0.0 : 6.0)
        group = shares[(:EU, :EU)] + shares[(:NA, :EU)]
        Pg = sum(shares[(o, :EU)] / group * q(o, :EU)^(1 - 12) for o in (:EU, :NA))^(1 / (1 - 12))
        @test price(res, :EU, :x) ≈
              (group * Pg^(1 - 1.5) + shares[(:AS, :EU)] * q(:AS, :EU)^(1 - 1.5))^(1 / (1 - 1.5)) rtol = 1e-6
        # ... while the AS price index is the flat one with σ = 8
        @test price(res, :AS, :x) ≈
              sum(shares[(o, :AS)] * q(o, :AS)^(1 - 8) for o in regions)^(1 / (1 - 8)) rtol = 1e-6
    end

    @testset "an :all specification is the fallback for the other destinations" begin
        res = solve_market(economy([
            Armington(product = :x, destination = :EU, sigma = 1.8, shares = shares),
            Armington(product = :x, sigma = 5, shares = shares)]))
        @test solved(res)
        for (r, σ) in ((:EU, 1.8), (:NA, 5), (:AS, 5))
            q(o) = pprice(res, o, :x) + (o == r ? 0.0 : 6.0)
            @test price(res, r, :x) ≈
                  sum(shares[(o, r)] * q(o)^(1 - σ) for o in regions)^(1 / (1 - σ)) rtol = 1e-6
        end
    end

    @testset "an infinitely substitutable group is pooled into one good" begin
        σ = 2.0
        res = solve_market(economy([Armington(product = :x, sigma = σ, shares = shares,
                                    nests = [OriginNest(sigma = Inf, origins = [:EU, :NA])])]))
        @test solved(res)
        # inside the pool the solution is the homogeneous one: no cross-hauling,
        # and no price gap wider than the freight between them
        @test min(flow(res, :x, :EU, :NA), flow(res, :x, :NA, :EU)) <= 1e-3
        @test abs(pprice(res, :EU, :x) - pprice(res, :NA, :x)) <= 6.0 + 1e-6
        # and the pool enters the composite at the cheapest delivered price of
        # its members, carrying their combined share
        for r in regions
            q(o)  = pprice(res, o, :x) + (o == r ? 0.0 : 6.0)
            pool  = min(q(:EU), q(:NA))
            δpool = shares[(:EU, r)] + shares[(:NA, r)]
            tot   = δpool + shares[(:AS, r)]
            @test price(res, r, :x) ≈ ((δpool / tot) * pool^(1 - σ) +
                  (shares[(:AS, r)] / tot) * q(:AS)^(1 - σ))^(1 / (1 - σ)) rtol = 1e-6
        end
    end

    @testset "a region that cannot produce a product has no variety of it" begin
        # :y is manufactured from :x, but only EU and NA have the process
        d = MarketData(
            regions   = regions, products = [:x, :y],
            demand    = [DemandSpec(product = :y, region = r, p0 = 200, q0 = 5, elasticity = 1.4)
                         for r in regions],
            supply    = [SupplySpec(product = :x, region = r, p0 = 60, q0 = 20, elasticity = 0.8)
                         for r in regions],
            processes = [Process(name = :mill, regions = [:EU, :NA], vacost = 30,
                                 inputs  = [leontief(product = :x, coeff = 1.5)],
                                 outputs = [:y => 1.0])],
            transport_costs = merge(transport_costs, Dict((:y, o, r) => 9.0 for o in regions, r in regions if o != r)),
            armington = [Armington(product = :y, sigma = 3)])
        res = solve_market(d)
        @test solved(res)
        # AS makes no :y, so it has no producer price for it and exports none
        @test ismissing(pprice(res, :AS, :y))
        @test !ismissing(pprice(res, :EU, :y)) && !ismissing(pprice(res, :NA, :y))
        @test all(r -> flow(res, :y, :AS, r) == 0.0, regions)
        @test price(res, :AS, :y) > 0            # it still buys, at a price
    end

    @testset "a region that uses none of a product needs no composite" begin
        # C neither consumes :y nor processes it, but it does produce it
        d = MarketData(
            regions   = regions, products = [:x],
            demand    = [DemandSpec(product = :x, region = r, p0 = 100, q0 = 10, elasticity = 1.5)
                         for r in (:EU, :NA)],                      # nothing in AS
            supply    = [SupplySpec(product = :x, region = r, p0 = 60, q0 = 15, elasticity = 1.0)
                         for r in regions],
            transport_costs = transport_costs,
            armington = [Armington(product = :x, sigma = 3)])
        res = solve_market(d)
        @test solved(res)                       # used to be infeasible
        @test ismissing(price(res, :AS, :x))    # nothing there to price a composite
        @test pprice(res, :AS, :x) > 0          # but its own output has a value
        @test getq(res.production, :AS, :x) > 1e-3
        @test !ismissing(price(res, :EU, :x))
    end

    @testset "a group with one available origin collapses" begin
        res = solve_market(economy([Armington(product = :x, sigma = 3, shares = shares,
                                    nests = [OriginNest(sigma = 9, origins = [:NA])])]))
        @test solved(res)
        flat = solve_market(economy([Armington(product = :x, sigma = 3, shares = shares)]))
        @test price(res, :EU, :x) ≈ price(flat, :EU, :x) rtol = 1e-5
    end

    @testset "specification errors" begin
        # a group must hold together at least as tightly as it does to the outside
        @test_throws ArgumentError Armington(product = :x, sigma = 5,
                                             nests = [OriginNest(sigma = 2, origins = [:EU, :NA])])
        @test_throws ArgumentError OriginNest(sigma = 4,
                                              nests = [OriginNest(sigma = 2, origins = [:EU, :NA])])
        @test_throws ArgumentError OriginNest(sigma = 0.5, origins = [:EU])
        @test_throws ArgumentError OriginNest(sigma = 4)                     # empty group
        @test_throws ArgumentError Armington(product = :x, sigma = 2,        # origin in two groups
                                             nests = [OriginNest(sigma = 9, origins = [:EU, :NA]),
                                                      OriginNest(sigma = 9, origins = [:NA, :AS])])
        # unknown region, duplicates, and incomplete destination coverage
        @test_throws ArgumentError solve_market(economy([
            Armington(product = :x, destination = :XX, sigma = 2)]))
        @test_throws ArgumentError solve_market(economy([
            Armington(product = :x, destination = :EU, sigma = 2),
            Armington(product = :x, destination = [:EU, :NA], sigma = 3)]))
        @test_throws ArgumentError solve_market(economy([        # nothing for NA and AS
            Armington(product = :x, destination = :EU, sigma = 2)]))
        @test_throws ArgumentError solve_market(economy([        # homogeneous in AS only
            Armington(product = :x, destination = [:EU, :NA], sigma = 2),
            Armington(product = :x, destination = :AS, sigma = Inf)]))
    end
end

@testset "Minimal example" begin
    # the README and the first documentation page quote these numbers; running
    # the script is what keeps them true
    ex = Module(:MinimalExample)
    Base.include(ex, joinpath(@__DIR__, "..", "examples", "minimal", "minimal_market.jl"))
    res = ex.result
    @test solved(res)
    @test nrow(res.prices) == 4
    # the two conditions the page checks by hand
    @test price(res, :Pasture, :cheese) ≈ 10 * price(res, :Pasture, :milk) + 1300 rtol = 1e-6
    @test price(res, :City, :milk) - price(res, :Pasture, :milk) ≈ 45 rtol = 1e-6
    @test price(res, :City, :cheese) - price(res, :Pasture, :cheese) ≈ 45 rtol = 1e-6
    # the data are a consistent benchmark, so the solver reproduces them
    @test price(res, :Pasture, :milk) ≈ 400 rtol = 1e-2
    @test price(res, :City, :cheese) ≈ 5350 rtol = 1e-2
    # the City makes no cheese at all: a tonne of it costs 45 to carry, the ten
    # tonnes of milk behind it 450
    @test isempty(res.production[(res.production.region .== :City) .&
                                 (res.production.product .== :cheese), :])
    @test getq(res.production, :Pasture, :cheese) > 100
    # milk crosses, but only the milk that is drunk
    @test flow(res, :milk, :Pasture, :City) > 1
    @test flow(res, :cheese, :Pasture, :City) > 1
    @test flow(res, :milk, :City, :Pasture) <= 1e-6
    @test getq(res.consumption, :City, :milk) ≈
          getq(res.production, :City, :milk) + flow(res, :milk, :Pasture, :City) rtol = 1e-6
    # the milk splits itself between the glass and the dairy
    milk_made  = sum(skipmissing(res.production.quantity[res.production.product .== :milk]))
    milk_drunk = sum(skipmissing(res.consumption.quantity[res.consumption.product .== :milk]))
    @test milk_made ≈ milk_drunk + 10 * getq(res.production, :Pasture, :cheese) rtol = 1e-6
    @test 0.3 < milk_drunk / milk_made < 0.4
    # a thirstier City makes cheese dearer everywhere, although no taste for
    # cheese changed, and takes milk away from the dairy
    for r in (:Pasture, :City)
        @test only(ex.thirst.prices[(ex.thirst.prices.region .== r) .&
                                    (ex.thirst.prices.product .== :cheese), :price]) >
              price(res, r, :cheese)
    end
    @test getq(ex.thirst.production, :Pasture, :cheese) <
          getq(res.production, :Pasture, :cheese)
end

@testset "Forest example" begin
    # the documented example is a runnable script: run it, then check that the
    # equilibrium it reports satisfies the conditions it claims
    ex = Module(:ForestExample)
    Base.include(ex, joinpath(@__DIR__, "..", "examples", "forest", "forest_market.jl"))
    d, res, france = ex.example_market, ex.res, ex.france
    @test solved(res)
    @test nrow(res.prices) == length(d.regions) * length(d.products)
    @test all(skipmissing(res.prices.price) .> 0)

    # papermill is Leontief: it buys the pulp composite and sells its own
    # variety of paper, so zero profit ties the producer price to the user one
    for r in d.regions
        @test pprice(res, r, :paper) ≈ 1.1 * price(res, r, :pulp) + 200 rtol = 1e-4
    end
    # no pulp mill in SEF and GEF, hence no pulp of their own to price
    @test ismissing(pprice(res, :SEF, :pulp)) && ismissing(pprice(res, :GEF, :pulp))
    @test !ismissing(pprice(res, :SWF, :pulp))

    # the French regions are perfect substitutes for each other, so between
    # them the solution is the homogeneous one: no cross-hauling ...
    intra = res.trade[in.(res.trade.from, Ref(france)) .& in.(res.trade.to, Ref(france)), :]
    for row in eachrow(intra)
        row.quantity > 1e-3 || continue
        @test flow(res, row.product, row.to, row.from) <= 1e-3
    end
    # ... and no price gap wider than the cost of shipping between them
    for p in d.products, a in france, b in france
        a < b || continue
        (ismissing(pprice(res, a, p)) || ismissing(pprice(res, b, p))) && continue
        @test abs(pprice(res, a, p) - pprice(res, b, p)) <=
              d.transport_costs[(p, a, b)] + 1e-6
    end
    # between blocs that no longer holds: the varieties are different goods
    @test any(abs(pprice(res, :NA, p) - pprice(res, :AS, p)) > d.transport_costs[(p, :NA, :AS)]
              for p in d.products)

    # world trade balances: total exports = total imports for each product
    for g in groupby(res.net_trade, :product)
        @test sum(g.net) ≈ 0 atol = 1e-5
    end

    # the storm scenario reaches France hardest, then the EU, then the far blocs
    swr(r, result) = only(result.prices[(result.prices.region .== r) .&
                                        (result.prices.product .== :swr), :price])
    drop(r) = 1 - swr(r, ex.storm) / swr(r, res)
    @test minimum(drop.(france)) > drop(:EU) > maximum(drop.([:NA, :AS, :RW]))
    @test all(drop(r) > 0 for r in d.regions)        # cheaper everywhere
end

end
