using Test
using Markets
using DataFrames
using JuMP: termination_status, objective_value, all_variables, name
import JuMP

# helpers to read single values out of the result tables.  The tables name
# things with strings, so these take either form and convert.
const S = string
price(res, r, p) = only(res.prices[(res.prices.region .== S(r)) .& (res.prices.product .== S(p)), :price])
pprice(res, r, p) = only(res.prices[(res.prices.region .== S(r)) .& (res.prices.product .== S(p)), :producer_price])
getq(df, r, p)   = (rows = df[(df.region .== S(r)) .& (df.product .== S(p)), :quantity]; isempty(rows) ? 0.0 : only(rows))
flow(res, p, from, to) = (rows = res.trade[(res.trade.product .== S(p)) .& (res.trade.from .== S(from)) .&
                                           (res.trade.to .== S(to)), :quantity]; isempty(rows) ? 0.0 : only(rows))
solved(res)      = res.status in (Markets.MOI.LOCALLY_SOLVED, Markets.MOI.OPTIMAL)

@testset "Markets.jl" begin

@testset "Constructors" begin
    n = leontief(product = :a, coeff = 2)
    @test n.products == ["a"] && n.composite == 2.0 && n.shares == [1.0]
    c = ces(composite = 1.5, products = [:a, :b], shares = [3, 1], sigma = 2)
    @test c.shares ≈ [0.75, 0.25]
    @test c.composite == 1.5 && c.sigma == 2.0 && c.phi == 1.0
    # shares default to equal ones, sigma must be ≥ 0, shares must align
    @test ces(composite = 1.0, products = [:a, :b, :c], sigma = 2).shares ≈ [1/3, 1/3, 1/3]
    @test_throws ArgumentError ces(composite = 1.0, products = [:a, :b], shares = [1.0], sigma = 2)
    @test_throws ArgumentError ces(composite = 1.0, products = [:a, :b], shares = [1.0, 1.0], sigma = -0.5)
    @test_throws ArgumentError ces(composite = 1.0, products = [:a, :b], sigma = NaN)
    # poor substitutes, Cobb–Douglas and fixed proportions are all CES nests
    for σ in (0, 0.5, 1)
        @test ces(composite = 1.0, products = [:a, :b], sigma = σ).sigma == σ
    end
    # a Nest built directly validates like one built by `ces`
    @test_throws ArgumentError Nest(1.0, [:a, :b], [1.0], 2.0, 1.0)      # shares mismatch
    @test_throws ArgumentError Nest(1.0, [:a, :b], [1.0, 1.0], -1.0, 1.0) # sigma negative
    @test_throws ArgumentError Nest(1.0, [:a, :a], [1.0, 1.0], 2.0, 1.0) # repeated product
    @test_throws ArgumentError Nest(0.0, [:a], [1.0], 0.0, 1.0)          # nothing needed
    @test leontief(product = :a, coeff = 2).sigma == 0.0                 # still fine
    # every field is a keyword, and the optional ones have defaults
    nest = Nest(composite = 2.0, products = [:a, :b], sigma = 3)
    @test nest.shares ≈ [0.5, 0.5] && nest.phi == 1.0
    @test_throws UndefKeywordError Nest(composite = 1.0, products = [:a])
    p = Process(name = :mill, inputs = [leontief(product = :a, coeff = 1)],
                outputs = [:b => 1, :c => 0.5], vacost = 10)
    @test p.outputs == ["b" => 1.0, "c" => 0.5] && p.regions isa Colon && p.vacost == 10.0
    @test_throws UndefKeywordError Process(name = :mill, inputs = Nest[])
    d = MarketData(regions = [:A], products = [:x])
    @test isempty(d.demand) && isempty(d.supply) && isempty(d.processes) &&
          isempty(d.transport_costs)
    # a curve needs a reference point and its elasticity, or a line's coefficients
    @test_throws ArgumentError DemandCurve(product = :x, region = :A, p0 = 1, q0 = 1)
end

@testset "Spatial price equilibrium (one product, two regions)" begin
    τ = 5.0
    d = MarketData(
        regions   = [:A, :B],
        products  = [:x],
        demand    = [DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 1.5),
                     DemandCurve(product = :x, region = :B, p0 = 100, q0 = 10, elasticity = 1.5)],
        supply    = [SupplyCurve(product = :x, region = :A, p0 = 50, q0 = 30, elasticity = 1.0),   # cheap
                     SupplyCurve(product = :x, region = :B, p0 = 150, q0 = 5, elasticity = 1.0)],  # expensive
        transport_costs = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ),
    )
    res = solve_market(d)
    @test solved(res)
    # A exports to B, and the price gap equals the transport cost
    @test only(res.trade.from) == "A" && only(res.trade.to) == "B"
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
        nt = res.net_trade[res.net_trade.region .== S(r), :net]
        @test getq(res.production, r, :x) - getq(res.consumption, r, :x) ≈ only(nt) atol = 1e-5
    end
    # for_region slices consistently
    a = for_region(res, :A)
    @test all(a.prices.region .== "A") && nrow(a.trade) == 1
end

@testset "The solver: choosing it, tuning it, and what it reports" begin
    d = MarketData(regions = [:A, :B], products = [:x],
        demand = [DemandCurve(product = :x, region = r, p0 = 100, q0 = 10, elasticity = 1.5)
                  for r in (:A, :B)],
        supply = [SupplyCurve(product = :x, region = :A, p0 = 50, q0 = 30, elasticity = 1.0),
                  SupplyCurve(product = :x, region = :B, p0 = 150, q0 = 5, elasticity = 1.0)],
        transport_costs = Dict((:x, :A, :B) => 5.0, (:x, :B, :A) => 5.0))
    res = solve_market(d)
    @test solved(res)
    @test res.raw_status == "Solve_Succeeded"
    @test res.solve_time > 0
    @test res.iterations isa Int && res.iterations > 0
    @test res.solves == 1
    # options reach the solver, in any of the forms accepted, and override the
    # package's own: one iteration is not enough to solve anything
    for options in ((max_iter = 1,), Dict("max_iter" => 1), [:max_iter => 1])
        cut = solve_market(d; solver_options = options, on_failure = :ignore)
        @test cut.status == Markets.MOI.ITERATION_LIMIT
        @test cut.iterations == 1
    end
    @test_throws ErrorException solve_market(d; solver_options = (max_iter = 1,))
    # an MOI attribute is passed as it is, and a solver may be named
    @test solved(solve_market(d; solver = Markets.Ipopt.Optimizer,
                              solver_options = Dict(Markets.MOI.TimeLimitSec() => 60.0)))
    # the package's own feasibility tolerance can be overridden by name
    loose = solve_market(d; solver_options = (constr_viol_tol = 1e-3,))
    @test JuMP.get_attribute(loose.model, "constr_viol_tol") == 1e-3
    # an Ipopt wrapped the JuMP way is still Ipopt: it gets the package's
    # tolerances, unless it sets its own
    ipopt(options...) = JuMP.optimizer_with_attributes(Markets.Ipopt.Optimizer, options...)
    tight = solve_market(d; solver = ipopt("tol" => 1e-9))
    @test solved(tight)
    @test JuMP.get_attribute(tight.model, "constr_viol_tol") ==
          JuMP.get_attribute(res.model, "constr_viol_tol")
    own = solve_market(d; solver = ipopt("constr_viol_tol" => 1e-3))
    @test JuMP.get_attribute(own.model, "constr_viol_tol") == 1e-3
end

@testset "Looking at an economy before solving it" begin
    d = MarketData(regions = [:A, :B, :C], products = [:log, :board, :chips, :panel],
        supply = [SupplyCurve(product = :log, region = r, p0 = 50, q0 = 10, elasticity = 1.0)
                  for r in (:A, :B)],
        demand = [[DemandCurve(product = :board, region = r, p0 = 200, q0 = 3, elasticity = 1.2)
                   for r in (:A, :B, :C)];
                  [DemandCurve(product = :panel, region = r, p0 = 200, q0 = 3, elasticity = 1.2)
                   for r in (:A, :B)]],
        processes = [Process(name = :sawmill, inputs = [leontief(product = :log, coeff = 2)],
                             outputs = [:board => 1.0, :chips => 0.6], vacost = 30),
                     Process(name = :panelmill, regions = :B,
                             inputs = [ces(composite = 1.2, products = [:chips, :log], sigma = 2)],
                             outputs = [:panel => 1.0], vacost = 20),
                     disposal(product = :chips)],
        transport_costs = Dict((:board, :A, :C) => 5.0, (:board, :C, :A) => 5.0,
                               (:board, :B, :C) => 7.0, (:panel, :B, :A) => 4.0),
        limits = [@limit("panels B", activity(region = :B, process = :panelmill) <= 4.0),
                  @limit("ferry", exports(region = :B) <= 2.0)],
        taxes = [Tax(name = :mill_aid, on = :process, process = :panelmill, rate = -3.0),
                 Tax(name = :duty, on = :trade, product = :board, to = :C, rate = 0.1,
                     relative = true)])

    @testset "the JuMP model is the one solve_market solves" begin
        m = jump_model(d)
        @test m isa JuMP.Model
        @test termination_status(m) == Markets.MOI.OPTIMIZE_NOT_CALLED
        # balances and limits are named, so the printed model can be read
        @test JuMP.constraint_by_name(m, "balance[A,board]") isa JuMP.ConstraintRef
        @test JuMP.constraint_by_name(m, "limit[ferry]") isa JuMP.ConstraintRef
        # a limit on one variable is a bound on it, not a constraint
        @test JuMP.constraint_by_name(m, "limit[panels B]") === nothing
        @test JuMP.upper_bound(JuMP.variable_by_name(m, "z[B,panelmill]")) == 4.0
        @test occursin("balance[C,board]", sprint(print, m))
        # C neither makes, uses nor carries panels: no balance on nothing
        @test !occursin("balance[C,panel]", sprint(print, m))
        # the ad valorem duty is found by iteration, so it is not in the model;
        # without it the model is exactly the first of solve_market's solves
        free = MarketData(; (f => getfield(d, f) for f in fieldnames(MarketData))...,
                          taxes = d.taxes[1:1])
        JuMP.optimize!(m)
        @test objective_value(m) ≈ solve_market(free).objective rtol = 1e-6
        @test JuMP.num_variables(m) == JuMP.num_variables(solve_market(d).model)
        # the solver keywords are those of solve_market, and validation runs first
        @test JuMP.get_attribute(jump_model(d; solver_options = (max_iter = 7,)),
                                 "max_iter") == 7
        @test_throws ArgumentError jump_model(MarketData(regions = [:A], products = [:x],
            demand = [DemandCurve(product = :y, region = :A, p0 = 1, q0 = 1, elasticity = 1)]))
    end

    @testset "the summary says what the economy is" begin
        # typed in the REPL it is one line, and so it is inside a container ...
        @test sprint(show, MIME"text/plain"(), d) == "MarketData(3 regions, 4 products, 3 processes)"
        @test repr(d) == "MarketData(3 regions, 4 products, 3 processes)"
        @test sprint(show, [d]) == "MarketData[MarketData(3 regions, 4 products, 3 processes)]"
        # ... while printed it is the whole summary, as it is when interpolated
        text = sprint(print, d)
        @test string(d) == text && "$d" == text && text != repr(d)
        @test startswith(text, "MarketData: 3 regions, 4 products, 3 processes, 4 trade routes, 2 taxes")
        # processes in the order they chain, with recipe, cost, bound and place
        @test occursin(r"1  sawmill\s+2 log ─+→ board \+ 0.6 chips", text)
        @test occursin(r"2  panelmill\s+1.2 CES\(chips, log; σ 2\) ─+→ panel\s+cost 20 · in B", text)
        @test occursin(r"2  dispose_chips\s+chips ─+→ \(disposed of\)", text)
        # a route priced the same both ways is shown once
        @test occursin("board  A⇄C 5, B→C 7", text)
        @test occursin("panel  B→A 4", text)
        @test occursin(r"panels B\s+≤ 4\s+activity\(region = B, process = panelmill\)", text)
        @test occursin(r"ferry\s+≤ 2\s+exports\(region = B\)", text)
        @test occursin(r"mill_aid\s+subsidy on process\s+3 per unit\s+panelmill", text)
        @test occursin(r"duty\s+tax on trade\s+10% of the origin price\s+board to C", text)
        # in a narrow terminal, a long route list becomes a count and a range
        many = MarketData(regions = [:A, :B, :C], products = [:x],
            transport_costs = Dict((:x, o, r) => 10.0 + (o < r) for o in ("A", "B", "C")
                                   for r in ("A", "B", "C") if o != r))
        wide = sprint(print, many)
        @test occursin("x  A→B 11, A→C 11, B→A 10, B→C 11, C→A 10, C→B 10", wide)
        narrow = sprint(print, many;
                        context = (:limit => true, :displaysize => (24, 40)))
        @test occursin("x  every route among 3 regions", narrow)
        @test all(textwidth(l) <= 40 for l in split(narrow, "\n"))
        # a process in one region is in that region, not in its letters
        @test occursin(r"in B$"m, text)
    end
end

@testset "A product may be supplied, consumed and used as an input at once" begin
    # :logs come out of the forest, and are either burnt as they are or sawn
    d = MarketData(
        regions   = [:A], products = [:logs, :board],
        demand    = [DemandCurve(product = :logs,  region = :A, p0 = 40,  q0 = 20, elasticity = 0.4),
                     DemandCurve(product = :board, region = :A, p0 = 300, q0 = 10, elasticity = 1.2)],
        supply    = [SupplyCurve(product = :logs, region = :A, p0 = 40, q0 = 35, elasticity = 0.6)],
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
        demand = [DemandCurve(product = :logs, region = :A, p0 = 70, q0 = 20, elasticity = 0.4),
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
        demand    = [DemandCurve(product = :out, region = :R, p0 = 500, q0 = 10, elasticity = 1.5)],
        supply    = [SupplyCurve(product = :in1, region = :R, p0 = 50, q0 = 10, elasticity = 0.8),
                     SupplyCurve(product = :in2, region = :R, p0 = 50, q0 = 10, elasticity = 0.8)],
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
        demand   = [DemandCurve(product = :x, region = r, p0 = 100, q0 = 10, elasticity = η)
                    for r in (:A, :B)],
        supply   = [SupplyCurve(product = :x, region = :A, p0 = 50,  q0 = 30, elasticity = 1.0),
                    SupplyCurve(product = :x, region = :B, p0 = 150, q0 = 5,  elasticity = 1.0)],
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
                 supply = [SupplyCurve(product = :x, region = :A, p0 = 50, q0 = 45, elasticity = 1.0),
                           SupplyCurve(product = :x, region = :B, p0 = 150, q0 = 5, elasticity = 1.0)],
                 transport_costs = Dict((:x, :A, :B) => 5.0, (:x, :B, :A) => 5.0))
        1 - price(solve_market(more), :A, :x) / price(base, :A, :x)
    end
    @test shock(0.4) > shock(1.5) > 0

    @test_throws ArgumentError DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10,
                                          elasticity = 0.0)
    @test_throws ArgumentError DemandCurve(product = :x, region = :A, p0 = -1, q0 = 10,
                                          elasticity = 1.5)
    @test_throws ArgumentError SupplyCurve(product = :x, region = :A, p0 = 50, q0 = 0,
                                          elasticity = 1.0)
    @test_throws ArgumentError SupplyCurve(product = :x, region = :A, p0 = 50, q0 = 10,
                                          elasticity = -0.5)
end

@testset "Linear curves" begin
    @testset "two ways to give a line, and what they mean" begin
        # through a reference point with the elasticity there
        d = DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 1.5,
                       form = :linear)
        @test d.form === :linear
        @test d.beta ≈ -1.5 * 10 / 100 && d.alpha ≈ 10 * (1 + 1.5)
        @test Markets.curve_price(d, 10.0) ≈ 100
        s = SupplyCurve(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 0.5,
                       form = "linear")
        @test s.beta ≈ 0.5 * 10 / 100 && s.alpha ≈ 10 * (1 - 0.5)
        @test Markets.curve_price(s, 10.0) ≈ 100
        # by its coefficients
        c = DemandCurve(product = :x, region = :A, alpha = 200, beta = -9)
        @test c.form === :linear && isnan(c.p0) && isnan(c.elasticity)
        @test Markets.curve_price(c, 200 - 9 * 12.0) ≈ 12
        # the default is still the constant elasticity, and a flat curve has one form
        @test DemandCurve(product = :x, region = :A, p0 = 1, q0 = 1, elasticity = 2).form === :isoelastic
        @test SupplyCurve(product = :x, region = :A, p0 = 1, q0 = 1, elasticity = Inf,
                         form = :linear).form === :isoelastic
        # the shifter scales the quantity at every price, in both forms
        shifted = DemandCurve(product = :x, region = :A, alpha = 200, beta = -9, shift = 2)
        @test Markets.curve_price(shifted, 2 * (200 - 9 * 12.0)) ≈ 12
        @test occursin("alpha = 200.0, beta = -9.0", repr(c))
        @test occursin("form = :linear", repr(d))
    end

    @testset "what a line refuses" begin
        @test_throws ArgumentError DemandCurve(product = :x, region = :A, alpha = 200, beta = 9)
        @test_throws ArgumentError DemandCurve(product = :x, region = :A, alpha = -1, beta = -9)
        @test_throws ArgumentError SupplyCurve(product = :x, region = :A, alpha = -50, beta = -1)
        @test_throws ArgumentError SupplyCurve(product = :x, region = :A, alpha = -50)
        @test_throws ArgumentError SupplyCurve(product = :x, region = :A, alpha = 1, beta = 2,
                                              p0 = 1)
        @test_throws ArgumentError SupplyCurve(product = :x, region = :A, alpha = 1, beta = 2,
                                              form = :isoelastic)
        @test_throws ArgumentError SupplyCurve(product = :x, region = :A, p0 = 1, q0 = 1,
                                              elasticity = 1, form = :quadratic)
    end

    # the equilibrium of lines can be found by hand: with A exporting to B,
    # pA + τ = pB and SA(pA) - DA(pA) = DB(pB) - SB(pB)
    τ = 2.0
    lines = MarketData(regions = [:A, :B], products = [:x],
        demand = [DemandCurve(product = :x, region = :A, alpha = 100, beta = -4),
                  DemandCurve(product = :x, region = :B, alpha = 200, beta = -9)],
        supply = [SupplyCurve(product = :x, region = :A, alpha = -50, beta = 20.5),
                  SupplyCurve(product = :x, region = :B, alpha = -50, beta = 10.5)],
        transport_costs = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ))
    res = solve_market(lines)
    @test solved(res)
    # excess supply of A, 24.5pA - 150, meets excess demand of B, 250 - 19.5(pA + τ)
    pA = (150 + 250 - 19.5τ) / (24.5 + 19.5)
    @test price(res, :A, :x) ≈ pA rtol = 1e-6
    @test price(res, :B, :x) ≈ pA + τ rtol = 1e-6
    @test getq(res.consumption, :A, :x) ≈ 100 - 4pA rtol = 1e-6
    @test getq(res.production, :B, :x) ≈ -50 + 10.5 * (pA + τ) rtol = 1e-6
    @test flow(res, :x, :A, :B) ≈ (-50 + 20.5pA) - (100 - 4pA) rtol = 1e-6

    # a line through the reference point agrees with the constant-elasticity
    # curve there, and responds to a shock by the slope it has there: for a
    # small shock the two forms move the price alike, for a large one they part
    market(form, q0) = MarketData(regions = [:A], products = [:x],
        demand = [DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 0.8,
                             form = form)],
        supply = [SupplyCurve(product = :x, region = :A, p0 = 100, q0 = q0, elasticity = 1.2,
                             form = form)])
    @test price(solve_market(market(:linear, 10)), :A, :x) ≈ 100 rtol = 1e-6
    move(form, q0) = price(solve_market(market(form, q0)), :A, :x) / 100 - 1
    @test move(:linear, 10.1) ≈ move(:isoelastic, 10.1) rtol = 1e-2
    @test !isapprox(move(:linear, 15), move(:isoelastic, 15); rtol = 1e-2)
    # a straight demand has a finite price at zero: no floor keeps it alive, so
    # a supply that only starts above the choke price leaves the market empty
    empty = solve_market(MarketData(regions = [:A], products = [:x],
        demand = [DemandCurve(product = :x, region = :A, alpha = 100, beta = -4)],      # chokes at 25
        supply = [SupplyCurve(product = :x, region = :A, alpha = -60, beta = 2)]))      # starts at 30
    @test solved(empty) && getq(empty.consumption, :A, :x) == 0
    # a tax on a straight supply: the market price covers the marginal cost
    # along the line plus the tax
    taxed = MarketData(regions = [:A], products = [:x],
        demand = [DemandCurve(product = :x, region = :A, alpha = 100, beta = -4)],
        supply = [SupplyCurve(product = :x, region = :A, alpha = -10, beta = 3)],
        taxes = [Tax(name = :levy, on = :supply, rate = 2.0)])
    t = solve_market(taxed)
    q = getq(t.production, :A, :x)
    @test price(t, :A, :x) ≈ (q + 10) / 3 + 2 rtol = 1e-5
    @test only(t.taxes.revenue) ≈ 2q rtol = 1e-5
end

@testset "Shifters on straight lines" begin
    # a shifter multiplies the quantity at every price, whatever the form: on a
    # line D = shift·(α + β·p) it scales α and β together, so the price at
    # which demand chokes stays where it was
    line(; kw...) = MarketData(regions = [:A], products = [:x],
        demand = [DemandCurve(product = :x, region = :A, alpha = 100, beta = -4; get(kw, :demand, ())...)],
        supply = [SupplyCurve(product = :x, region = :A, p0 = 10, q0 = 30, elasticity = 0.5,
                             form = :linear; get(kw, :supply, ())...)])
    base = solve_market(line())
    for m in (0.6, 1.5, 3.0)
        res = solve_market(line(supply = (shift = m,)))
        s = line(supply = (shift = m,)).supply[1]
        P, Q = price(res, :A, :x), getq(res.production, :A, :x)
        @test Q ≈ m * (s.alpha + s.beta * P) rtol = 1e-6          # on the shifted line
        @test getq(res.consumption, :A, :x) ≈ 100 - 4P rtol = 1e-6
        @test Markets.curve_price(s, Q) ≈ P rtol = 1e-6
    end
    @test price(solve_market(line(supply = (shift = 1.5,))), :A, :x) < price(base, :A, :x)
    @test price(solve_market(line(demand = (shift = 1.4,))), :A, :x) > price(base, :A, :x)
    shifted = line(demand = (shift = 2.0,)).demand[1]
    @test Markets.curve_price(shifted, 0.0) ≈ Markets.curve_price(line().demand[1], 0.0) ≈ 25
    # exogenous_shift gives the same elasticity of the quantity to each driver
    # at every price, on a line as on a constant-elasticity curve
    z = exogenous_shift(rain = (level = 1.2, reference = 1.0, elasticity = 0.5))
    plain  = SupplyCurve(product = :x, region = :A, alpha = -20, beta = 4)
    wetter = SupplyCurve(product = :x, region = :A, alpha = -20, beta = 4, shift = z)
    for P in (6.0, 10.0, 30.0)
        Q0 = plain.alpha + plain.beta * P
        @test Markets.curve_price(wetter, 1.2^0.5 * Q0) ≈ P
    end
    # a scenario re-shifts a curve without touching its calibration or its form
    moved = SupplyCurve(plain; shift = 2.5)
    @test moved.form === :linear && moved.alpha == plain.alpha && moved.beta == plain.beta &&
          moved.shift == 2.5
    tangent = DemandCurve(product = :x, region = :A, p0 = 5, q0 = 9, elasticity = 0.4,
                         form = :linear)
    @test DemandCurve(tangent; shift = 0.5) ==
          DemandCurve(product = :x, region = :A, p0 = 5, q0 = 9, elasticity = 0.4,
                     form = :linear, shift = 0.5)
    @test_throws ArgumentError SupplyCurve(plain; shift = 0)
    # a flat curve has no quantity to scale: a shifter on it is an error, in
    # the constructor and when re-shifting
    @test_throws ArgumentError SupplyCurve(product = :x, region = :world, p0 = 9, q0 = 1,
                                          elasticity = Inf, shift = 1.2)
    flat = DemandCurve(product = :x, region = :world, p0 = 9, q0 = 1, elasticity = Inf)
    @test_throws ArgumentError DemandCurve(flat; shift = 2)
    @test DemandCurve(flat; shift = 1) == flat
    # where the price is at the reference point, the two forms move alike
    for form in (:isoelastic, :linear)
        s = SupplyCurve(product = :x, region = :A, p0 = 10, q0 = 30, elasticity = 0.5,
                       form = form, shift = z)
        @test Markets.curve_price(s, 1.2^0.5 * 30) ≈ 10
    end
end

@testset "Flat curves: a price given from outside" begin
    # a small country facing the world: its own curves at 200 in autarky, and
    # a "world" region that sells at `buy` and buys at `sell`, whatever
    # the quantity (elasticity Inf), each 10 away by transport
    τ = 10.0
    economy(buy, sell) = MarketData(
        regions  = [:home, :world], products = [:wheat],
        demand   = [DemandCurve(product = :wheat, region = :home, p0 = 200, q0 = 100, elasticity = 0.5),
                    DemandCurve(product = :wheat, region = :world, p0 = sell, q0 = 1, elasticity = Inf)],
        supply   = [SupplyCurve(product = :wheat, region = :home, p0 = 200, q0 = 100, elasticity = 0.8),
                    SupplyCurve(product = :wheat, region = :world, p0 = buy, q0 = 1, elasticity = Inf)],
        transport_costs = Dict((:wheat, :world, :home) => τ, (:wheat, :home, :world) => τ))
    on_curve(res) = price(res, :home, :wheat) ≈
        200 * (getq(res.consumption, :home, :wheat) / 100)^(-1 / 0.5)

    # cheap imports: the price at home is import parity, and nothing goes out
    cheap = solve_market(economy(170, 150))
    @test solved(cheap)
    @test price(cheap, :home, :wheat) ≈ 170 + τ rtol = 1e-6
    @test price(cheap, :world, :wheat) ≈ 170 rtol = 1e-6
    @test flow(cheap, :wheat, :world, :home) > 0
    @test flow(cheap, :wheat, :home, :world) == 0
    @test getq(cheap.consumption, :world, :wheat) == 0     # no floor forces a sale
    @test on_curve(cheap)
    # dear exports: the price at home is export parity
    dear = solve_market(economy(260, 230))
    @test price(dear, :home, :wheat) ≈ 230 - τ rtol = 1e-6
    @test flow(dear, :wheat, :home, :world) > 0
    @test flow(dear, :wheat, :world, :home) == 0
    @test on_curve(dear)
    # world prices on either side of 200: neither trade pays, and home is in
    # autarky, at the price its own curves set
    apart = solve_market(economy(230, 170))
    @test price(apart, :home, :wheat) ≈ 200 rtol = 1e-4
    @test nrow(apart.trade) == 0
    # the world selling for less than it buys: no finite answer
    @test_throws ErrorException solve_market(economy(150, 170))
    # a flat curve has no quantity of its own: its q0 does not set the scale
    # of the economy, and so neither the floors nor the solver's tolerances
    tiny(q0) = MarketData(regions = [:home, :world], products = [:wheat],
        demand = [DemandCurve(product = :wheat, region = :home, p0 = 200, q0 = 100, elasticity = 0.5)],
        supply = [SupplyCurve(product = :wheat, region = :home, p0 = 200, q0 = 100, elasticity = 0.8),
                  SupplyCurve(product = :wheat, region = :world, p0 = 170, q0 = q0, elasticity = Inf)],
        transport_costs = Dict((:wheat, :world, :home) => τ))
    @test Markets.quantity_scale(tiny(1e-9)) == Markets.quantity_scale(tiny(1e9)) == 100
end

@testset "Exogenous drivers of the curves" begin
    economy(; supply_kw = (), demand_kw = (), limits = Limit[]) = MarketData(
        regions = [:A], products = [:x],
        demand  = [DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10,
                              elasticity = 1.5; demand_kw...)],
        supply  = [SupplyCurve(product = :x, region = :A, p0 = 50, q0 = 20,
                              elasticity = 0.8; supply_kw...)],
        limits  = limits)
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
        ceiling(c) = economy(limits = [@limit "cut" supply(region = :A, product = :x) <= c])
        res  = solve_market(ceiling(cap))
        @test solved(res)
        S = getq(res.production, :A, :x)
        @test S ≈ cap rtol = 1e-6                       # the bound binds
        @test price(res, :A, :x) > price(free, :A, :x)  # scarcity raises the price
        # the price now exceeds the marginal cost on the curve: that gap is the
        # scarcity rent of the resource, reported on the supply and the limit alike
        @test price(res, :A, :x) > on_curve(economy().supply[1], S)
        @test only(res.supply.rent) ≈ price(res, :A, :x) - on_curve(economy().supply[1], S) rtol = 1e-5
        @test only(res.supply.capacity) == cap
        @test only(res.limits.rent) ≈ only(res.supply.rent) rtol = 1e-10
        @test only(res.limits.regions) == ["A"] && only(res.limits.variables) == ["supply[A, x]"]
        # a ceiling above what would be produced anyway changes nothing
        loose = solve_market(ceiling(10 * cap))
        @test getq(loose.production, :A, :x) ≈ getq(free.production, :A, :x) rtol = 1e-6
        @test only(loose.supply.rent) == 0
        @test only(loose.limits.rent) == 0
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
        @test_throws ArgumentError SupplyCurve(product = :x, region = :A, p0 = 1, q0 = 1,
                                              elasticity = 1, shift = 0)
        @test_throws ArgumentError DemandCurve(product = :x, region = :A, p0 = 1, q0 = 1,
                                              elasticity = 2, shift = -0.5)
    end
end

@testset "Validation: the economy must say what it means" begin
    ## a small, correct economy that each case below breaks in one way
    good(; kw...) = MarketData(;
        regions = [:A, :B], products = [:logs, :board],
        demand  = [DemandCurve(product = :board, region = :A, p0 = 300, q0 = 10,
                              elasticity = 1.2)],
        supply  = [SupplyCurve(product = :logs, region = :A, p0 = 40, q0 = 35,
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
            [DemandCurve(product = :nope, region = :A, p0 = 1, q0 = 1, elasticity = 2)]))
        @test_throws ArgumentError solve_market(good(supply =
            [SupplyCurve(product = :logs, region = :Mars, p0 = 1, q0 = 1, elasticity = 1)]))
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
            [SupplyCurve(product = :logs, region = :A, p0 = 40, q0 = 35, elasticity = 0.6),
             SupplyCurve(product = :logs, region = :A, p0 = 90, q0 = 35, elasticity = 0.6)]))
        @test_throws ArgumentError solve_market(good(demand =
            [DemandCurve(product = :board, region = :A, p0 = 300, q0 = 10, elasticity = 1.2),
             DemandCurve(product = :board, region = :A, p0 = 300, q0 = 10, elasticity = 1.2)]))
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
        @test_throws ArgumentError solve_market(good(limits =
            [@limit "l" flow(product = :nope) <= 1.0]))
        @test_throws ArgumentError solve_market(good(limits =
            [@limit "l" exports(region = :Mars) <= 1.0]))
    end

    @testset "economies with no finite answer" begin
        # a process with no inputs and no cost is a free source of its output
        @test_throws ArgumentError solve_market(MarketData(regions = [:A], products = [:g],
            demand = [DemandCurve(product = :g, region = :A, p0 = 300, q0 = 10,
                                 elasticity = 1.2)],
            processes = [Process(name = :magic, vacost = 0.0, inputs = Nest[],
                                 outputs = [:g => 1.0])]))
        # the same process is legitimate once it costs something ...
        @test solved(solve_market(MarketData(regions = [:A], products = [:g],
            demand = [DemandCurve(product = :g, region = :A, p0 = 300, q0 = 10,
                                 elasticity = 1.2)],
            processes = [Process(name = :backstop, vacost = 50.0, inputs = Nest[],
                                 outputs = [:g => 1.0])])))
        # ... or once it is bounded
        @test solved(solve_market(MarketData(regions = [:A], products = [:g],
            demand = [DemandCurve(product = :g, region = :A, p0 = 300, q0 = 10,
                                 elasticity = 1.2)],
            processes = [Process(name = :windfall, vacost = 0.0, inputs = Nest[],
                                 outputs = [:g => 1.0])],
            limits = [@limit "plant" activity(process = :windfall) <= 4.0])))
        @test_throws ArgumentError solve_market(MarketData(regions = Symbol[],
            products = [:g]))
        @test_throws ArgumentError solve_market(MarketData(regions = [:A],
            products = Symbol[]))
    end

    @testset "a solver verdict is not a suggestion" begin
        # demand with no source anywhere is infeasible; it used to be returned
        # as LOCALLY_SOLVED because the violation sat under Ipopt's tolerance
        starved = MarketData(regions = [:A], products = [:x],
            demand = [DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10,
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
    economy(; limits = Limit[]) = MarketData(
        regions   = [:A, :B], products = [:w, :g],
        demand    = [DemandCurve(product = :g, region = r, p0 = 300, q0 = 10, elasticity = 1.2)
                     for r in (:A, :B)],
        supply    = [SupplyCurve(product = :w, region = :A, p0 = 40, q0 = 40, elasticity = 0.8),
                     SupplyCurve(product = :w, region = :B, p0 = 90, q0 = 5,  elasticity = 0.8)],
        processes = [Process(name = :mill, vacost = 50,
                             inputs  = [leontief(product = :w, coeff = 2.0)],
                             outputs = [:g => 1.0])],
        transport_costs = Dict((p, o, r) => (p === :w ? 6.0 : 10.0)
                         for p in (:w, :g), o in (:A, :B), r in (:A, :B) if o != r),
        limits = limits)
    level(res, r) = only(res.activity[res.activity.region .== S(r), :level])

    @testset "a capacity holds a process back" begin
        free = solve_market(economy())
        @test solved(free)
        @test all(isinf, free.activity.capacity) && all(iszero, free.activity.rent)
        cap = 0.6 * level(free, :A)
        res = solve_market(economy(limits = [@limit "plant A" activity(region = :A, process = :mill) <= cap]))
        @test solved(res)
        @test level(res, :A) ≈ cap rtol = 1e-5          # the bound binds in A ...
        @test level(res, :B) > level(free, :B)          # ... and work moves to B
        row = only(res.activity[res.activity.region .== "A", :])
        @test row.rent > 0 && row.capacity == cap      # a quasi-rent, on the activity ...
        @test row.rent ≈ only(res.limits.rent) rtol = 1e-10   # ... and on the limit
        @test only(res.limits.variables) == ["activity[A, mill]"]
        # a limit in one region bounds that region only
        @test isinf(only(res.activity[res.activity.region .== "B", :capacity]))
        @test iszero(only(res.activity[res.activity.region .== "B", :rent]))
    end

    @testset "a value-added cost may differ by region" begin
        regional(cost) = MarketData(
            regions   = [:A, :B], products = [:w, :g],
            demand    = [DemandCurve(product = :g, region = r, p0 = 300, q0 = 10, elasticity = 1.2)
                         for r in (:A, :B)],
            supply    = [SupplyCurve(product = :w, region = r, p0 = 40, q0 = 40, elasticity = 0.8)
                         for r in (:A, :B)],
            processes = [Process(name = :mill, vacost = cost,
                                 inputs  = [leontief(product = :w, coeff = 2.0)],
                                 outputs = [:g => 1.0])],
            transport_costs = Dict((p, o, r) => 30.0
                                   for p in (:w, :g), o in (:A, :B), r in (:A, :B) if o != r))
        res = solve_market(regional(Dict(:A => 50, :B => 80)))
        @test solved(res)
        # zero profit, each region at its own cost
        for (r, c) in ((:A, 50), (:B, 80))
            @test price(res, r, :g) ≈ 2 * price(res, r, :w) + c rtol = 1e-6
        end
        @test price(res, :B, :g) > price(res, :A, :g)
        # one number for every region is the same as a dictionary of equal ones
        same = solve_market(regional(Dict(:A => 50, :B => 50)))
        flat = solve_market(regional(50))
        @test price(same, :B, :g) ≈ price(flat, :B, :g) rtol = 1e-8
        @test occursin("cost 50–80", sprint(print, regional(Dict(:A => 50, :B => 80))))
        # a region the process runs in must have a cost, and only declared ones may
        @test_throws ArgumentError validate(regional(Dict(:A => 50)))
        @test_throws ArgumentError validate(regional(Dict(:A => 50, :B => 80, :C => 1)))
    end

    @testset "a limit on any quantity" begin
        free = solve_market(economy())
        supplied(res, r) = only(res.supply[res.supply.region .== S(r), :quantity])
        cap = 0.8 * supplied(free, :A)
        capped(c) = solve_market(economy(limits = [@limit "stands" supply(region = :A, product = :w) <= c]))
        res = capped(cap)
        @test solved(res)
        row = only(eachrow(res.limits))
        @test row.limit == "stands" && row.upper == cap && row.lower == -Inf
        @test row.value ≈ cap rtol = 1e-5
        @test supplied(res, :A) ≈ cap rtol = 1e-5
        # the rent is what one more unit of the bound is worth to the objective
        δ = 1e-3
        @test row.rent > 0
        @test row.rent ≈ (capped(cap + δ).objective - res.objective) / δ rtol = 1e-2
        # built without the macro, the same limit gives the same equilibrium
        same = solve_market(economy(limits = [Limit(name = "stands", upper = cap,
                                    expression = Markets.supply(region = :A, product = :w))]))
        @test same.objective ≈ res.objective rtol = 1e-8

        with(limits...) = solve_market(economy(limits = collect(Limit, limits)))
        # a minimum is met at a cost, and its rent is negative: a subsidy
        floor = 1.5 * supplied(free, :B)
        res = with(@limit "harvest B" supply(region = :B, product = :w) >= floor)
        @test only(res.limits.value) ≈ floor rtol = 1e-5
        @test only(res.limits.rent) < 0
        # an equality fixes the expression
        quota = 0.5 * only(free.trade[free.trade.product .== "g", :quantity])
        res = with(@limit "quota" exports(region = :A, product = :g) == quota)
        @test only(res.limits.value) ≈ quota rtol = 1e-5
        @test only(res.limits.lower) == only(res.limits.upper) == quota
        # quantities on both sides, a constant, and a range
        res = with(@limit("milled at home", activity(region = :A, process = :mill) <=
                                            0.4 * production(region = :A, product = :w) - 1),
                   @limit("band", 2.0 <= supply(region = :B, product = :w) <= 3.0))
        @test solved(res)
        mill = only(res.activity[res.activity.region .== "A", :level])
        @test mill <= 0.4 * supplied(res, :A) - 1 + 1e-6
        @test 2.0 - 1e-6 <= supplied(res, :B) <= 3.0 + 1e-6
        band = only(res.limits[res.limits.limit .== "band", :])
        @test band.lower == 2.0 && band.upper == 3.0
        # a family is one limit per combination of its indices, named after them
        regions_, goods = [:A, :B], [:w]
        fam = @limits "cut" [r in regions_, p in goods] supply(region = r, product = p) <= 30
        @test getfield.(fam, :name) == ["cut[A, w]", "cut[B, w]"]
        @test getfield.(@limits("only", [r ∈ regions_; r == :B], supply(region = r) <= 1), :name) == ["only[B]"]
        @test getfield.(@limits("eq", [r = regions_], supply(region = r) <= 1), :name) == ["eq[A]", "eq[B]"]
        res = with(fam...)
        @test solved(res) && nrow(res.limits[startswith.(res.limits.limit, "cut["), :]) == 2
        # every region or product is `:`, and leaving a selector out means the same
        @test Markets.scope(@limit("w", supply(region = :, product = :w) <= 99), economy()).variables ==
              Markets.scope(@limit("w", supply(product = :w) <= 99), economy()).variables ==
              ["supply[A, w]", "supply[B, w]"]
        # what a limit applies to is known before solving, and reported after
        shared = @limit("shared", activity(region = :A, process = :mill) +
                                  production(region = :B, product = :w) <= 30)
        sc = Markets.scope(shared, economy())
        @test sc.regions == ["A", "B"]
        @test sc.variables == ["activity[A, mill]", "supply[B, w]"]
        res = with(shared)
        @test only(res.limits.regions) == sc.regions && only(res.limits.variables) == sc.variables
        # a group exports only what leaves it
        @test_throws ArgumentError with(@limit "inside" exports(region = [:A, :B], product = :g) <= 1)
        # it is part of the economy's summary
        @test occursin("Limits", sprint(print, economy(limits = [@limit "stands" supply(region = :A) <= 1])))

        # what it refuses: names the economy does not have, quantities that
        # do not exist, two limits of one name, and anything not linear
        @test_throws ArgumentError with(@limit "x" supply(region = :Mars) <= 1)
        @test_throws ArgumentError with(@limit "x" consumption(region = :A, product = :w) <= 1)
        @test_throws ArgumentError with(@limit("x", supply(region = :A) <= 1),
                                        @limit("x", supply(region = :B) <= 1))
        @test_throws ArgumentError Markets.supply(product = :w) * Markets.supply(product = :g)
        @test_throws ArgumentError 1 / Markets.supply(product = :w)
        @test_throws ArgumentError sqrt(Markets.supply(product = :w))
        @test_throws ArgumentError @limit "x" 1 <= 2
        @test_throws ArgumentError Limit(name = "x", expression = Markets.supply(), lower = 2, upper = 1)
        @test_throws ArgumentError Limit(name = "x", expression = Markets.supply())
    end

    @testset "a soft limit can be crossed, at its penalty" begin
        free = solve_market(economy())
        made = level(free, :A)
        cap(penalty) = @limit "mill A" activity(region = :A, process = :mill) <= 0.5made penalty = penalty
        hard = solve_market(economy(limits = [@limit "mill A" activity(region = :A, process = :mill) <= 0.5made]))
        hard_rent = only(hard.limits.rent)
        @test hard_rent > 0
        # a penalty above the rent the hard limit earns is never worth paying:
        # the soft limit holds where the hard one does, at the same rent
        dear = solve_market(economy(limits = [cap(2hard_rent)]))
        @test level(dear, :A) ≈ 0.5made rtol = 1e-6
        @test only(dear.limits.rent) ≈ hard_rent rtol = 1e-4
        @test only(dear.limits.excess) ≈ 0 atol = 1e-6
        @test dear.penalties ≈ 0 atol = 1e-4
        # a cheaper one is paid: the limit is crossed, its rent is the penalty,
        # and the mill runs where its margin covers that much
        cheap = solve_market(economy(limits = [cap(0.25hard_rent)]))
        row = only(eachrow(cheap.limits))
        @test level(cheap, :A) > 0.5made
        @test row.excess ≈ level(cheap, :A) - 0.5made rtol = 1e-5
        @test row.rent ≈ 0.25hard_rent rtol = 1e-4
        @test row.penalty == 0.25hard_rent
        @test cheap.penalties ≈ row.penalty * row.excess rtol = 1e-6
        @test cheap.welfare ≈ cheap.objective + cheap.penalties rtol = 1e-9
        margin = price(cheap, :A, :g) - 2 * price(cheap, :A, :w) - 50
        @test margin ≈ row.rent rtol = 1e-4
        # the rent never exceeds the penalty, and the mill sits between
        @test level(hard, :A) < level(cheap, :A) < made
        # a bound not reached changes nothing
        slack = solve_market(economy(limits = [@limit "mill A" activity(region = :A, process = :mill) <=
                                               2made penalty = 1.0]))
        @test level(slack, :A) ≈ made rtol = 1e-6
        @test only(slack.limits.rent) == 0 && only(slack.limits.excess) == 0
        # a penalty of Inf is the hard limit, and a soft ceiling is no capacity
        @test cap(Inf).penalty == Inf
        @test level(solve_market(economy(limits = [cap(Inf)])), :A) ≈ level(hard, :A) rtol = 1e-9
        @test only(cheap.activity[cheap.activity.region .== "A", :capacity]) == Inf
        @test only(hard.activity[hard.activity.region .== "A", :capacity]) ≈ 0.5made
        # a soft floor: short of it, the penalty is paid on the shortfall, and
        # the rent is negative, as for any lower bound
        floor = solve_market(economy(limits = [@limit "mill B" activity(region = :B, process = :mill) >=
                                               level(free, :B) + 50 penalty = 5.0]))
        r = only(eachrow(floor.limits))
        @test r.rent ≈ -5 rtol = 1e-4
        @test r.excess ≈ level(free, :B) + 50 - r.value rtol = 1e-5
        @test r.value > level(free, :B)
        # a family takes the penalty for each member
        fam = @limits "mills" [r in (:A, :B)] activity(region = r, process = :mill) <= 1.0 penalty = 3
        @test all(l.penalty == 3 for l in fam) && length(fam) == 2
        @test occursin("penalty 3.0", repr(fam[1]))
        @test occursin("a unit past it", sprint(print, economy(limits = fam)))
        # what a penalty must be, and what else a limit takes
        @test_throws ArgumentError cap(0.0)
        @test_throws ArgumentError cap(-1.0)
        @test_throws Exception macroexpand(@__MODULE__,
            :(@limit "x" activity(region = :A) <= 1 cost = 3))
    end

    @testset "a trade ceiling binds and earns a congestion rent" begin
        free  = solve_market(economy())
        moved = only(free.trade[free.trade.product .== S(:g), :quantity])
        quota = 0.55 * moved
        res = solve_market(economy(limits = [@limit "quota" flow(product = :g, from = :A, to = :B) <= quota]))
        @test solved(res)
        @test flow(res, :g, :A, :B) ≈ quota rtol = 1e-4
        row = only(eachrow(res.limits))
        @test row.limit == "quota"
        @test row.value ≈ quota rtol = 1e-4
        @test row.upper == quota && row.lower == -Inf
        @test row.rent > 0
        # with the ceiling binding, the price gap is transport *plus* the rent,
        # which the flow itself reports
        @test price(res, :B, :g) - price(res, :A, :g) ≈ 10.0 + row.rent rtol = 1e-4
        moved_g = only(res.trade[res.trade.product .== "g", :])
        @test moved_g.rent ≈ row.rent rtol = 1e-10
        @test moved_g.capacity == quota
        @test row.regions == ["A", "B"] && row.variables == ["flow[g, A, B]"]
    end

    @testset "a ceiling over several products uses its weights" begin
        res = solve_market(economy(limits = [@limit "port" 0.5 * exports(region = :A, product = :w) +
                                                           exports(region = :A, product = :g) <= 6.0]))
        @test solved(res)
        used = 0.5 * flow(res, :w, :A, :B) + flow(res, :g, :A, :B)
        @test only(res.limits.value) ≈ used rtol = 1e-5
        @test used ≈ 6.0 rtol = 1e-4
        # each flow it covers is charged the rent times its weight — the boards
        # that cross carry all of it — and none has a ceiling of its own
        rent = only(res.limits.rent)
        @test only(res.trade[res.trade.product .== "g", :rent]) ≈ rent rtol = 1e-8
        @test only(res.limits.variables) == ["flow[g, A, B]", "flow[w, A, B]"]
        @test all(isinf, res.trade.capacity)
        # and the region's slice holds it
        @test nrow(for_region(res, :A).limits) == 1
        # a loose ceiling changes nothing
        loose = solve_market(economy(limits = [@limit "port" exports(region = :A) <= 1e4]))
        @test flow(loose, :g, :A, :B) ≈ flow(solve_market(economy()), :g, :A, :B) rtol = 1e-5
        @test only(loose.limits.rent) == 0
    end
end

@testset "Substitution in final demand" begin
    nest = DemandNest(sigma = 4.0, products = [:pellets, :firewood], shares = [0.6, 0.4])
    economy(sh) = MarketData(
        regions = [:A], products = [:pellets, :firewood],
        demand  = [DemandCurve(product = nest, region = :A, p0 = 100, q0 = 20,
                              elasticity = 0.5)],
        supply  = [SupplyCurve(product = :pellets,  region = :A, p0 = 100, q0 = 12,
                              elasticity = 0.8),
                   SupplyCurve(product = :firewood, region = :A, p0 = 100, q0 = 8,
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
                   demand = [DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10,
                                        elasticity = 1.2)],
                   supply = [SupplyCurve(product = :x, region = :A, p0 = 50, q0 = 20,
                                        elasticity = 0.8)])
        nested = MarketData(regions = plain.regions, products = plain.products,
                   demand = [DemandCurve(product = DemandNest(sigma = 3.0, products = [:x]),
                                        region = :A, p0 = 100, q0 = 10, elasticity = 1.2)],
                   supply = plain.supply)
        @test getq(solve_market(nested).consumption, :A, :x) ≈
              getq(solve_market(plain).consumption, :A, :x) rtol = 1e-6
    end

    @testset "nesting rules" begin
        @test_throws ArgumentError DemandNest(sigma = -0.5, products = [:a])
        @test_throws ArgumentError DemandNest(sigma = 2.0)                       # empty
        @test_throws ArgumentError DemandNest(sigma = 2.0, products = [:a, :b],
                                              shares = [1.0])                     # mismatch
        @test_throws ArgumentError DemandNest(sigma = 4.0, products = [:a],
                   nests = [DemandNest(sigma = 2.0, products = [:b])])             # loose inside
        # a product may not have two demand curves in one region
        twice = MarketData(regions = [:A], products = [:x],
            demand = [DemandCurve(product = :x, region = :A, p0 = 1, q0 = 1, elasticity = 2),
                      DemandCurve(product = :x, region = :A, p0 = 1, q0 = 1, elasticity = 2)],
            supply = [SupplyCurve(product = :x, region = :A, p0 = 1, q0 = 1, elasticity = 1)])
        @test_throws ArgumentError solve_market(twice)
    end
end

@testset "Taxes and subsidies" begin
    economy(taxes) = MarketData(
        regions  = [:A, :B], products = [:x],
        demand   = [DemandCurve(product = :x, region = r, p0 = 100, q0 = 10, elasticity = 1.2)
                    for r in (:A, :B)],
        supply   = [SupplyCurve(product = :x, region = :A, p0 = 50,  q0 = 30, elasticity = 0.8),
                    SupplyCurve(product = :x, region = :B, p0 = 150, q0 = 5,  elasticity = 0.8)],
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

    @testset "the scarcity rent is net of the tax" begin
        q0  = getq(base.production, :A, :x)
        mc(q) = 50 * (q / 30)^(1 / 0.8)
        with_cap(taxes, cap) = (e = economy(taxes);
            push!(e.limits, @limit("cut", supply(region = :A, product = :x) <= cap)); e)
        for (tax, taxed) in ((Tax(on = :supply, region = :A, rate = 6.0), q -> mc(q) + 6.0),
                             (Tax(on = :supply, region = :A, rate = 0.2, relative = true),
                              q -> 1.2 * mc(q)))
            # a slack ceiling earns nothing, however the production is taxed
            slack = solve_market(with_cap([tax], 3 * q0))
            @test solved(slack)
            @test getq(slack.supply, :A, :x) < 3 * q0
            @test only(slack.supply[slack.supply.region .== "A", :rent]) == 0
            # a binding one earns the gap between the price and the *taxed* cost
            tight = solve_market(with_cap([tax], 0.5 * q0))
            @test solved(tight)
            q = getq(tight.supply, :A, :x)
            @test q ≈ 0.5 * q0 rtol = 1e-6
            @test only(tight.supply[tight.supply.region .== "A", :rent]) ≈
                  price(tight, :A, :x) - taxed(q) rtol = 1e-5
        end
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
            demand    = [DemandCurve(product = :g, region = :A, p0 = 300, q0 = 10,
                                    elasticity = 1.2)],
            supply    = [SupplyCurve(product = :w, region = :A, p0 = 40, q0 = 40,
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
            demand = [DemandCurve(product = nest, region = :A, p0 = 10, q0 = 2, elasticity = 2)],
            supply = [SupplyCurve(product = p, region = :A, p0 = 10, q0 = 1, elasticity = 1)
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
        demand    = [DemandCurve(product = :g, region = :A, p0 = 300, q0 = 10, elasticity = 1.5)],
        supply    = [SupplyCurve(product = :w, region = :A, p0 = 50, q0 = 20, elasticity = 1.0)],
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
            demand = [DemandCurve(product = :x,    region = :A, p0 = 100, q0 = 50,
                                 elasticity = 1.3),
                      DemandCurve(product = :dust, region = :A, p0 = 5, q0 = 0.001,
                                 elasticity = 1.3)],
            supply = [SupplyCurve(product = :x,    region = :A, p0 = 40, q0 = 120,
                                 elasticity = 0.8),
                      SupplyCurve(product = :dust, region = :A, p0 = 5, q0 = 0.001,
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
            demand = [DemandCurve(product = :x, region = r, p0 = 100, q0 = 10,
                                 elasticity = eta) for r in (:A, :B)],
            supply = [SupplyCurve(product = :x, region = :A, p0 = 30,  q0 = 40,
                                 elasticity = 3.0),
                      SupplyCurve(product = :x, region = :B, p0 = 200, q0 = 2,
                                 elasticity = 3.0)],
            transport_costs = Dict((:x, :A, :B) => 1.0, (:x, :B, :A) => 1.0),
            taxes = [Tax(on = :trade, product = :x, from = :A, to = :B, rate = rate,
                         relative = true)])
        for rate in (0.2, 0.5, 1.0, 2.0, 4.0)        # up to a 400% tariff
            res = solve_market(economy(rate, 2.0))
            @test solved(res)
            # found by repeated solves, whose cost is reported in full
            @test res.solves > 1 && res.iterations >= res.solves
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
        demand = [DemandCurve(product = :g, region = :A, p0 = 300, q0 = 10,
                             elasticity = 1.5)],
        supply = [SupplyCurve(product = :w, region = :A, p0 = 50, q0 = 20,
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
    @test "dispose_junk" in solve_market(economy([disposal(product = :junk)])).activity.process
    # and the good products are unaffected in sign
    @test price(free, :A, :g) > 0 && price(free, :A, :w) > 0
end

@testset "Armington share calibration" begin
    flows  = Dict(("EU", "EU") => 80.0,  ("NA", "EU") => 15.0,  ("AS", "EU") => 5.0,
                  ("NA", "NA") => 70.0,  ("EU", "NA") => 20.0,  ("AS", "NA") => 10.0)
    prices = Dict(("EU", "EU") => 210.0, ("NA", "EU") => 235.0, ("AS", "EU") => 250.0,
                  ("NA", "NA") => 200.0, ("EU", "NA") => 240.0, ("AS", "NA") => 245.0)
    σ = 4
    δ = armington_shares(; flows, prices, sigma = σ)
    # shares sum to one per destination ...
    for dest in ("EU", "NA")
        @test sum(v for (k, v) in δ if k[2] == dest) ≈ 1
    end
    # ... and reproduce the observed sourcing through the CES demand condition
    for dest in ("EU", "NA")
        keys_d  = [k for k in keys(δ) if k[2] == dest]
        implied = Dict(k => δ[k] * prices[k]^(-σ) for k in keys_d)
        tot, ftot = sum(values(implied)), sum(flows[k] for k in keys_d)
        for k in keys_d
            @test implied[k] / tot ≈ flows[k] / ftot rtol = 1e-8
        end
    end
    # they are not the observed shares themselves, which is the point
    @test !isapprox(δ[("EU", "EU")], 0.8; rtol = 1e-3)
    @test_throws ArgumentError armington_shares(; flows, prices, sigma = -0.5)
    @test_throws ArgumentError armington_shares(; flows, prices, sigma = Inf)
    # at σ = 0 the shares are the quantity shares, at σ = 1 the value shares
    for (σ, weight) in ((0, k -> flows[k]), (1, k -> flows[k] * prices[k]))
        δσ = armington_shares(; flows, prices, sigma = σ)
        for k in keys(flows)
            @test δσ[S.(k)] ≈ weight(k) / sum(weight(j) for j in keys(flows) if j[2] == k[2])
        end
    end
    @test_throws ArgumentError armington_shares(flows = Dict((:EU, :EU) => 1.0),
                                                prices = Dict{Tuple{Symbol,Symbol},Float64}(),
                                                sigma = 3)
end

@testset "Armington calibration from a trade matrix" begin
    # A small world solved once with arbitrary shares stands in for the base
    # year: its flows and prices are what a user would observe. Rebuilt from
    # those alone, the economy must give them back.
    regions = [:A, :B, :C, :D]
    τ(o, r) = o == r ? 0.0 : 4.0 + 3.0 * abs(findfirst(==(o), regions) - findfirst(==(r), regions))
    nests = [OriginNest(sigma = 8.0, origins = [:A], nests = [OriginNest(sigma = Inf, origins = [:B, :C])])]
    truth = MarketData(regions = regions, products = [:x],
        demand = [DemandCurve(product = :x, region = r, p0 = 100, q0 = q0, elasticity = 0.8)
                  for (r, q0) in zip(regions, (12.0, 8.0, 15.0, 10.0))],
        supply = [SupplyCurve(product = :x, region = r, p0 = p0, q0 = q0, elasticity = 1.2)
                  for (r, p0, q0) in zip(regions, (60, 55, 58, 70), (10.0, 14.0, 9.0, 12.0))],
        transport_costs = Dict((:x, o, r) => τ(o, r) for o in regions, r in regions if o != r),
        armington = [Armington(product = :x, sigma = 3.0, nests = nests,
                               shares = Dict((o, r) => 1.0 + 0.3 * i + 0.1 * j
                                             for (i, o) in enumerate(regions),
                                                 (j, r) in enumerate(regions)))])
    base = solve_market(truth; solver_options = (tol = 1e-12,))
    @test solved(base)

    made     = Dict(S(r) => getq(base.production, r, :x) for r in regions)
    received = Dict(S(r) => pprice(base, r, :x) for r in regions)
    flows = Dict((S(o), S(r)) => flow(base, :x, o, r) for o in regions, r in regions if o != r)
    for r in S.(regions)
        flows[(r, r)] = made[r] - sum(x for ((o, _), x) in flows if o == r)
    end
    filter!(f -> last(f) > 1e-7, flows)
    delivered = Dict((S(o), S(r)) => received[S(o)] + τ(o, r) for o in regions, r in regions)

    arm = armington_calibration(product = :x, flows = flows, prices = delivered, sigma = 3.0,
                                nests = nests, rtol = 1e-6)
    @test arm.nests == nests && arm.sigma == 3.0
    for r in S.(regions)
        @test sum(v for ((_, d), v) in arm.shares if d == r) ≈ 1
    end
    # the consumer price the base year implies: the average delivered unit value
    bought(r) = sum(x for ((_, d), x) in flows if d == r)
    unit_value(r) = sum(x * delivered[k] for (k, x) in flows if k[2] == r) / bought(r)
    rebuilt = MarketData(regions = regions, products = [:x],
        demand = [DemandCurve(product = :x, region = r, p0 = unit_value(S(r)), q0 = bought(S(r)),
                             elasticity = 0.8) for r in regions],
        supply = [SupplyCurve(product = :x, region = r, p0 = received[S(r)], q0 = made[S(r)],
                             elasticity = 1.2) for r in regions],
        transport_costs = truth.transport_costs, armington = [arm])
    res = solve_market(rebuilt; solver_options = (tol = 1e-12,))
    @test solved(res)
    for r in regions
        @test pprice(res, r, :x) ≈ received[S(r)] rtol = 1e-6
        @test price(res, r, :x) ≈ unit_value(S(r)) rtol = 1e-6
        @test getq(res.production, r, :x) ≈ made[S(r)] rtol = 1e-5
        # the composite is counted in the units of the flows
        @test getq(res.consumption, r, :x) ≈ bought(S(r)) rtol = 1e-5
    end
    # the same prices as the economy it was observed from, which counted its
    # composites in index units
    for r in regions
        @test pprice(res, r, :x) ≈ pprice(base, r, :x) rtol = 1e-6
    end
    # flows from outside the pool come back one by one; inside it only its
    # total to each destination is determined
    for o in (:A, :D), r in regions
        o == r && continue
        @test flow(res, :x, o, r) ≈ flow(base, :x, o, r) rtol = 1e-4 atol = 1e-6
    end

    # the calibration does not depend on the shape of the curves: rebuilt with
    # the tangent lines at the same base year, the economy gives it back too
    lines = MarketData(regions = regions, products = [:x],
        demand = [DemandCurve(product = :x, region = r, p0 = unit_value(S(r)), q0 = bought(S(r)),
                             elasticity = 0.8, form = :linear) for r in regions],
        supply = [SupplyCurve(product = :x, region = r, p0 = received[S(r)], q0 = made[S(r)],
                             elasticity = 1.2, form = :linear) for r in regions],
        transport_costs = truth.transport_costs, armington = [arm])
    straight = solve_market(lines; solver_options = (tol = 1e-12,))
    @test solved(straight)
    for r in regions
        @test pprice(straight, r, :x) ≈ received[S(r)] rtol = 1e-6
        @test price(straight, r, :x) ≈ unit_value(S(r)) rtol = 1e-6
        @test getq(straight.production, r, :x) ≈ made[S(r)] rtol = 1e-5
        @test getq(straight.consumption, r, :x) ≈ bought(S(r)) rtol = 1e-5
    end

    # a pool member that would be the cheapest but shipped nothing contradicts
    # the data, and a flow without a price cannot be calibrated
    cheap = copy(delivered)
    pooled = [k for k in keys(flows) if k[1] in ("B", "C") && k[2] == "A"]
    idle = only(k for k in (("B", "A"), ("C", "A")) if !haskey(flows, k) || flows[k] < 1e-7)
    cheap[idle] = 0.5 * minimum(delivered[k] for k in pooled)
    @test_throws ArgumentError armington_calibration(product = :x, flows = flows, prices = cheap,
                                                     sigma = 3.0, nests = nests, rtol = 1e-6)
    @test_throws ArgumentError armington_calibration(product = :x, flows = flows,
        prices = Dict(k => v for (k, v) in delivered if k != ("D", "A")), sigma = 3.0)
    @test_throws ArgumentError armington_calibration(product = :x, flows = flows,
                                                     prices = delivered, sigma = Inf)
    # the scale is checked like the shares, and a missing destination gets 1
    @test_throws ArgumentError Armington(product = :x, sigma = 3, scale = 0.0)
    @test Markets.scale_of(Armington(product = :x, sigma = 3, scale = Dict(:A => 2.0)), "B") == 1.0
end

@testset "Physical accounts" begin
    # wood is an Armington product, drunk by households and milled together
    # with scrap, a homogeneous one, in a CES bundle
    res = solve_market(MarketData(
        regions = [:A, :B], products = [:wood, :scrap, :board],
        demand  = [DemandCurve(product = :wood, region = :A, p0 = 60, q0 = 8, elasticity = 0.8),
                   DemandCurve(product = :wood, region = :B, p0 = 60, q0 = 6, elasticity = 0.8),
                   DemandCurve(product = :board, region = :A, p0 = 300, q0 = 3, elasticity = 1.2)],
        supply  = [SupplyCurve(product = :wood, region = :A, p0 = 40, q0 = 12, elasticity = 1.0),
                   SupplyCurve(product = :wood, region = :B, p0 = 45, q0 = 10, elasticity = 1.0),
                   SupplyCurve(product = :scrap, region = :A, p0 = 30, q0 = 4, elasticity = 1.0)],
        processes = [Process(name = :mill, regions = [:A], vacost = 50,
                             inputs = [ces(composite = 2.0, products = [:wood, :scrap],
                                           shares = [0.6, 0.4], sigma = 2.0)],
                             outputs = [:board => 1.0])],
        transport_costs = Dict((p, o, r) => 5.0 for p in (:wood, :scrap, :board),
                               o in (:A, :B), r in (:A, :B) if o != r),
        armington = [Armington(product = :wood, sigma = 3.0,
                               shares = Dict((:A, :A) => 0.7, (:B, :A) => 0.3,
                                             (:B, :B) => 0.6, (:A, :B) => 0.4))]))
    @test solved(res)
    comp(r) = only(res.composites[res.composites.region .== S(r), :])
    for r in (:A, :B)
        c = comp(r)
        # received is the sum of the varieties delivered, the region's own included
        @test c.received ≈ getq(res.production, r, :wood) - flow(res, :wood, r, r == :A ? :B : :A) +
                           flow(res, :wood, r == :A ? :B : :A, r) rtol = 1e-6
        @test c.ratio ≈ c.received / c.composite
        @test c.ratio > 1 + 1e-4           # a mixed, unscaled composite counts for less
        # every use is credited with the region's ratio, and the uses add up to
        # what was received
        cons = only(res.consumption[(res.consumption.region .== S(r)) .&
                                    (res.consumption.product .== "wood"), :])
        @test cons.physical ≈ cons.quantity * c.ratio
        drawn = res.inputs[(res.inputs.region .== S(r)) .& (res.inputs.product .== "wood"), :]
        @test cons.physical + sum(drawn.physical; init = 0.0) ≈ c.received rtol = 1e-6
    end
    # a homogeneous product has no composite and no gap of its own
    @test !any(res.composites.product .== "scrap")
    @test only(res.consumption[res.consumption.product .== "board", :physical]) ≈
          only(res.consumption[res.consumption.product .== "board", :quantity])
    # the recipe of the mill against its physical inputs
    rec = only(res.recipes[res.recipes.process .== "mill", :])
    @test rec.inputs == "wood, scrap"
    @test rec.recipe ≈ 2.0 * only(res.activity[res.activity.process .== "mill", :level])
    @test rec.physical ≈ sum(res.inputs[res.inputs.process .== "mill", :physical])
    @test rec.gap ≈ rec.physical - rec.recipe
    @test abs(rec.gap) > 1e-4
    # the region slice carries the physical tables too
    @test nrow(for_region(res, :A).composites) == 1
end

@testset "Material accounts" begin
    # logs, an Armington product harvested in both regions, are sawn into
    # boards and chips in A, and chips and logs pulped there; boards and pulp
    # are homogeneous and traded
    economy = MarketData(
        regions = [:A, :B], products = [:logs, :chips, :board, :pulp],
        demand  = [DemandCurve(product = p, region = r, p0 = p0, q0 = q0, elasticity = 1.0)
                   for (p, p0, q0) in ((:board, 300, 4), (:pulp, 500, 2)) for r in (:A, :B)],
        supply  = [SupplyCurve(product = :logs, region = :A, p0 = 40, q0 = 20, elasticity = 1.0),
                   SupplyCurve(product = :logs, region = :B, p0 = 45, q0 = 12, elasticity = 1.0)],
        processes = [Process(name = :sawmill, regions = [:A], vacost = 30,
                             inputs = [leontief(product = :logs, coeff = 1.0)],
                             outputs = [:board => 0.5, :chips => 0.3]),
                     Process(name = :pulpmill, regions = [:A], vacost = 60,
                             inputs = [ces(composite = 2.0, products = [:chips, :logs],
                                           shares = [0.5, 0.5], sigma = 2.0)],
                             outputs = [:pulp => 1.0])],
        transport_costs = Dict((p, o, r) => 5.0 for p in (:logs, :chips, :board, :pulp),
                               o in (:A, :B), r in (:A, :B) if o != r),
        armington = [Armington(product = :logs, sigma = 4.0,
                               shares = Dict((:A, :A) => 0.7, (:B, :A) => 0.3,
                                             (:B, :B) => 0.7, (:A, :B) => 0.3))])
    res = solve_market(economy)
    @test solved(res)
    carbon = Material(name = :carbon, content = Dict(:logs => 0.5, (:B, :logs) => 0.6),
                      losses = Dict(:sawmill => 0.2, :pulpmill => 0.5),
                      allocation = Dict(:sawmill => Dict(:board => 1.0, :chips => 2.0)))
    acc = material_accounts(res, carbon)
    intensity(r, p, col) = only(acc.intensities[(acc.intensities.region .== S(r)) .&
                                                (acc.intensities.product .== S(p)), col])

    # every row balances, and so does the world: what was extracted was
    # consumed or lost
    b = acc.balances
    @test b.extracted .+ b.made .+ b.imported .- b.exported ≈ b.consumed .+ b.processed
    @test sum(b.extracted) ≈ sum(b.consumed) + sum(acc.processes.lost) rtol = 1e-10
    @test all(acc.processes.output .≈ acc.processes.input .- acc.processes.lost)
    # the content of primary supply is the one given, region by region, and
    # what a region uses is a mix of what it received
    @test intensity(:A, :logs, :produced) ≈ 0.5
    @test intensity(:B, :logs, :produced) ≈ 0.6
    @test 0.5 < intensity(:A, :logs, :used) < 0.6
    @test sum(b[(b.region .== "A") .& (b.product .== "logs"), :extracted]) ≈
          0.5 * only(res.supply[res.supply.region .== "A", :quantity])
    # the sawmill keeps four fifths of what enters, and gives a unit of chips
    # twice the content of a unit of board
    @test intensity(:A, :chips, :produced) ≈ 2 * intensity(:A, :board, :produced) rtol = 1e-8
    saw = only(acc.processes[acc.processes.process .== "sawmill", :])
    @test saw.lost ≈ 0.2 * saw.input
    # homogeneous boards carry the content of where they were made
    @test intensity(:B, :board, :used) ≈ intensity(:A, :board, :produced) rtol = 1e-8

    # what it refuses
    @test_throws ArgumentError Material(name = :c, content = Dict(:logs => -1.0))
    @test_throws ArgumentError Material(name = :c, content = Dict(:logs => 1.0),
                                        losses = Dict(:sawmill => 1.5))
    @test_throws ArgumentError material_accounts(res, Material(name = :c,
                                                 content = Dict(:wood => 1.0)))
    @test_throws ArgumentError material_accounts(res, Material(name = :c,
                                                 content = Dict(:logs => 1.0),
                                                 allocation = Dict(:sawmill => Dict(:pulp => 1.0))))
end

@testset "Imputing a trade matrix" begin
    regions = ["A", "B", "C"]
    pos = Dict("A" => 0.0, "B" => 1.0, "C" => 2.5)
    dist(o, r) = abs(pos[o] - pos[r])
    made = Dict("A" => 10.0, "B" => 6.0, "C" => 4.0)
    used = Dict("A" => 5.0, "B" => 8.0, "C" => 7.0)
    sums(f, side, r) = sum(x for (k, x) in f if k[side] == r)

    # production and use alone: the totals hold, and the flows are the gravity
    # pattern scaled by one factor per origin and one per destination
    f = impute_flows(production = made, use = used, distance = dist, decay = 0.8, home_bias = 2.0)
    for r in regions
        @test sums(f, 1, r) ≈ made[r] rtol = 1e-8
        @test sums(f, 2, r) ≈ used[r] rtol = 1e-8
    end
    g(o, r) = exp(-0.8 * dist(o, r)) * (o == r ? 2.0 : 1.0)
    @test f[("A", "B")] * f[("B", "C")] / (f[("A", "C")] * f[("B", "B")]) ≈
          g("A", "B") * g("B", "C") / (g("A", "C") * g("B", "B")) rtol = 1e-6

    # known totals are met: one flow, and a group over several origins
    known = Dict(("A", "C") => 1.0, (["A", "B"], "B") => 6.5)
    f = impute_flows(production = made, use = used, distance = dist, decay = 0.8,
                     totals = known)
    @test f[("A", "C")] ≈ 1.0 rtol = 1e-8
    @test f[("A", "B")] + f[("B", "B")] ≈ 6.5 rtol = 1e-8
    @test all(sums(f, 1, r) ≈ made[r] for r in regions)

    # uses that do not add up to production are scaled until they do
    f = impute_flows(production = made, use = Dict(r => 2q for (r, q) in used),
                     distance = dist, decay = 0.8)
    @test all(sums(f, 2, r) ≈ used[r] for r in regions)

    # distances as a dictionary, in either order; a region that makes or uses
    # nothing has no row or no column
    d = Dict(("A", "B") => 1.0, ("C", "A") => 2.5, ("B", "C") => 1.5)
    f = impute_flows(production = merge(made, Dict("D" => 0.0)), use = used,
                     distance = d, decay = 0.8)
    @test sums(f, 1, "A") ≈ 10.0 rtol = 1e-8
    @test !any(k -> k[1] == "D", keys(f))

    # contradictions and nonsense
    @test_throws ArgumentError impute_flows(production = made, use = used, distance = dist,
                                            decay = 0.8, totals = Dict(("A", "B") => 20.0))
    @test_throws ArgumentError impute_flows(production = made, use = used, distance = d,
                                            decay = 0.8, totals = Dict(("A", "Z") => 1.0))
    @test_throws ArgumentError impute_flows(production = made, use = used, distance = dist,
                                            decay = -1)
    @test_throws ArgumentError impute_flows(production = Dict("A" => -1.0), use = used,
                                            distance = dist, decay = 1)
    @test_throws ArgumentError impute_flows(production = made, use = used,
                                            distance = Dict(("A", "B") => 1.0), decay = 1)

    # an imputed matrix is what armington_calibration takes
    arm = armington_calibration(product = :x, sigma = 4.0, flows = f,
                                prices = Dict(k => 100.0 + 10dist(k...) for k in keys(f)))
    @test all(sum(v for (k, v) in arm.shares if k[2] == r) ≈ 1 for r in regions)
end

@testset "CES calibration" begin
    q, p = [7.0, 1.4, 3.0], [55.0, 250.0, 120.0]
    spending = sum(q .* p)

    @testset "σ = $σ inverts the CES demand condition" for σ in (0.3, 3.0)
        c = ces_calibration(quantities = q, prices = p, sigma = σ)
        @test sum(c.shares) ≈ 1
        @test c.shares ≈ (q .* p .^ σ) ./ sum(q .* p .^ σ)
        @test c.price ≈ sum(c.shares .* p .^ (1 - σ))^(1 / (1 - σ))
        @test c.price * c.quantity ≈ spending
        # at the base-year prices, the composite buys back the base-year quantities
        @test c.quantity .* c.shares .* (c.price ./ p) .^ σ ≈ q
    end

    @testset "the limits at 0, 1 and Inf" begin
        c0 = ces_calibration(quantities = q, prices = p, sigma = 0)
        @test c0.shares ≈ q ./ sum(q)                      # fixed proportions
        @test c0.price ≈ sum(c0.shares .* p) && c0.quantity ≈ sum(q)
        c1 = ces_calibration(quantities = q, prices = p, sigma = 1)
        @test c1.shares ≈ (q .* p) ./ spending              # Cobb–Douglas: value shares
        @test c1.price ≈ prod(p .^ c1.shares) && c1.price * c1.quantity ≈ spending
        for σ in (1 - 1e-6, 1 + 1e-6)                       # and it is the limit there
            near = ces_calibration(quantities = q, prices = p, sigma = σ)
            @test near.shares ≈ c1.shares rtol = 1e-5
            @test near.price ≈ c1.price rtol = 1e-5
        end
        # perfect substitutes: bought together only at one price
        same = ces_calibration(quantities = q, prices = fill(80.0, 3), sigma = Inf)
        @test same.price == 80.0 && same.quantity ≈ sum(q) && same.shares ≈ q ./ sum(q)
        @test ces_calibration(quantities = [2.0, 0.0], prices = [80.0, 95.0],
                              sigma = Inf).price == 80.0  # a dearer member nobody buys
        @test_throws ArgumentError ces_calibration(quantities = q, prices = p, sigma = Inf)
    end

    @testset "the scale that counts the composite in the members' unit" begin
        for σ in (0.0, 0.5, 1.0, 3.0)
            c = ces_calibration(quantities = q, prices = p, sigma = σ)
            @test c.scale ≈ sum(q) / c.quantity
            # priced per unit of the sum, the composite costs the average unit price
            @test c.price / c.scale ≈ spending / sum(q)
        end
        @test ces_calibration(quantities = q, prices = fill(80.0, 3), sigma = Inf).scale == 1
        # a process bundle built on it uses the base-year inputs per unit of activity
        c = ces_calibration(quantities = q, prices = p, sigma = 2.0)
        n = ces(composite = sum(q), products = [:a, :b, :c], shares = c.shares, sigma = 2.0,
                phi = c.scale)
        @test n.phi * (sum(n.shares .^ (1 / 2) .* q .^ (1 / 2)))^2 ≈ n.composite
        # perfect substitutes bought at prices equal only to their last digits
        near = [80.0, 80.0 * (1 + 1e-7)]
        @test_throws ArgumentError ces_calibration(quantities = [1.0, 2.0], prices = near,
                                                   sigma = Inf)
        @test ces_calibration(quantities = [1.0, 2.0], prices = near, sigma = Inf,
                              rtol = 1e-6).price == 80.0
    end

    @testset "what it refuses" begin
        @test_throws ArgumentError ces_calibration(quantities = q, prices = p, sigma = -1)
        @test_throws ArgumentError ces_calibration(quantities = q, prices = p, sigma = NaN)
        @test_throws ArgumentError ces_calibration(quantities = q, prices = p[1:2], sigma = 2)
        @test_throws ArgumentError ces_calibration(quantities = [1.0, -1.0], prices = [1.0, 1.0], sigma = 2)
        @test_throws ArgumentError ces_calibration(quantities = [0.0, 0.0], prices = [1.0, 1.0], sigma = 2)
        @test_throws ArgumentError ces_calibration(quantities = [1.0, 1.0], prices = [1.0, 0.0], sigma = 2)
        # a member nobody buys gets no share
        @test ces_calibration(quantities = [1.0, 0.0], prices = [1.0, 2.0], sigma = 2).shares == [1.0, 0.0]
    end

    # supply curves through the base-year points, demand from the calibration:
    # the equilibrium must be the base year itself
    supply_at(names, q, p) = [SupplyCurve(product = n, region = :A, p0 = p[i], q0 = q[i],
                                         elasticity = 1.0) for (i, n) in enumerate(names)]
    function reproduces(res, names, q, p)
        @test solved(res)
        for (i, n) in enumerate(names)
            @test getq(res.consumption, :A, n) ≈ q[i] rtol = 1e-5
            @test price(res, :A, n) ≈ p[i] rtol = 1e-5
        end
    end

    # the composite's curve may be either form: they agree at the base year
    @testset "a demand nest at σ = $σ reproduces the base year ($form)" for
            σ in (0.0, 0.5, 1.0, 3.0), form in (:isoelastic, :linear)
        names = [:a, :b, :c]
        c = ces_calibration(quantities = q, prices = p, sigma = σ)
        res = solve_market(MarketData(regions = [:A], products = names,
            demand = [DemandCurve(region = :A, p0 = c.price, q0 = c.quantity, elasticity = 0.6,
                                 form = form,
                                 product = DemandNest(sigma = σ, products = names,
                                                      shares = c.shares))],
            supply = supply_at(names, q, p)))
        reproduces(res, names, q, p)
    end

    @testset "a nested bundle is calibrated from the inside out (σ $σin inside, $σout outside)" for
            (σin, σout) in ((5.0, 1.6), (1.0, 0.0))
        names = [:a, :b, :c]
        inner = ces_calibration(quantities = q[2:3], prices = p[2:3], sigma = σin)
        outer = ces_calibration(quantities = [q[1], inner.quantity],
                                prices = [p[1], inner.price], sigma = σout)
        res = solve_market(MarketData(regions = [:A], products = names,
            demand = [DemandCurve(region = :A, p0 = outer.price, q0 = outer.quantity,
                                 elasticity = 0.6,
                                 product = DemandNest(sigma = σout, products = [:a],
                                                      shares = [outer.shares[1]],
                                                      nests = [DemandNest(sigma = σin,
                                                          share = outer.shares[2],
                                                          products = [:b, :c],
                                                          shares = inner.shares)]))],
            supply = supply_at(names, q, p)))
        reproduces(res, names, q, p)
    end
end

@testset "Armington trade" begin
    τ = 5.0
    # A is the cheap producer, B the expensive one
    economy(arm) = MarketData(
        regions   = [:A, :B],
        products  = [:x],
        demand    = [DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 1.5),
                     DemandCurve(product = :x, region = :B, p0 = 100, q0 = 10, elasticity = 1.5)],
        supply    = [SupplyCurve(product = :x, region = :A, p0 = 50,  q0 = 30, elasticity = 1.0),
                     SupplyCurve(product = :x, region = :B, p0 = 150, q0 = 5,  elasticity = 1.0)],
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
            supply    = [SupplyCurve(product = :x, region = :A, p0 = 50, q0 = 39, elasticity = 1.0),
                         SupplyCurve(product = :x, region = :B, p0 = 150, q0 = 5, elasticity = 1.0)],
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
            demand    = [DemandCurve(product = :f, region = :A, p0 = 400, q0 = 10, elasticity = 1.4),
                         DemandCurve(product = :f, region = :B, p0 = 400, q0 = 10, elasticity = 1.4)],
            supply    = [SupplyCurve(product = :w, region = :A, p0 = 60, q0 = 30, elasticity = 0.7),
                         SupplyCurve(product = :w, region = :B, p0 = 90, q0 = 15, elasticity = 0.7)],
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
        @test_throws ArgumentError Armington(product = :x, sigma = -0.5)
        @test_throws ArgumentError Armington(product = :x, sigma = NaN)
        @test_throws ArgumentError OriginNest(sigma = -1.0, origins = [:A])
        @test_throws ArgumentError Armington(product = :x, sigma = 2,
                                             shares = Dict((:A, :A) => -1.0))
        @test_throws ArgumentError solve_market(economy([Armington(product = :nope, sigma = 2)]))
        @test_throws ArgumentError solve_market(economy([Armington(product = :x, sigma = 2),
                                                         Armington(product = :x, sigma = 3)]))
        # no origin left with a positive share for destination A
        @test_throws ArgumentError solve_market(economy([Armington(product = :x, sigma = 2,
                                                   shares = Dict((:A, :B) => 1.0))]))
    end

    @testset "the scarcity rent is read off the producer price" begin
        with_cap(cap) = (e = economy([Armington(product = :x, sigma = 3)]);
            isfinite(cap) && push!(e.limits, @limit("cut", supply(region = :A, product = :x) <= cap)); e)
        q0 = getq(solve_market(with_cap(Inf)).supply, :A, :x)
        # slack: users pay more than producers get, and that gap is no rent
        slack = solve_market(with_cap(3 * q0))
        @test solved(slack)
        @test price(slack, :A, :x) > pprice(slack, :A, :x) + 1
        @test only(slack.supply[slack.supply.region .== "A", :rent]) == 0
        # binding: the rent is the producer price less the cost on the curve
        tight = solve_market(with_cap(0.5 * q0))
        @test solved(tight)
        q = getq(tight.supply, :A, :x)
        @test only(tight.supply[tight.supply.region .== "A", :rent]) ≈
              pprice(tight, :A, :x) - 50 * q / 30 rtol = 1e-5
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
        demand   = [DemandCurve(product = :x, region = r, p0 = 100, q0 = 10, elasticity = 1.5)
                    for r in regions],
        supply   = [SupplyCurve(product = :x, region = :EU, p0 = 60, q0 = 20, elasticity = 1.0),
                    SupplyCurve(product = :x, region = :NA, p0 = 55, q0 = 25, elasticity = 1.0),
                    SupplyCurve(product = :x, region = :AS, p0 = 70, q0 = 15, elasticity = 1.0)],
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

    @testset "a `:` specification is the fallback for the other destinations" begin
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
            demand    = [DemandCurve(product = :y, region = r, p0 = 200, q0 = 5, elasticity = 1.4)
                         for r in regions],
            supply    = [SupplyCurve(product = :x, region = r, p0 = 60, q0 = 20, elasticity = 0.8)
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
            demand    = [DemandCurve(product = :x, region = r, p0 = 100, q0 = 10, elasticity = 1.5)
                         for r in (:EU, :NA)],                      # nothing in AS
            supply    = [SupplyCurve(product = :x, region = r, p0 = 60, q0 = 15, elasticity = 1.0)
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
        @test_throws ArgumentError OriginNest(sigma = -0.5, origins = [:EU])
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

@testset "Elasticities of substitution of 1 and below" begin
    # the CES price index and its limit at σ = 1 (at σ = 0 the formula is
    # already the Leontief Σδp)
    index(δ, p, σ) = isone(σ) ? prod(p .^ δ) : sum(δ .* p .^ (1 - σ))^(1 / (1 - σ))
    # 0: fixed proportions, 0.3: poor substitutes, 1: Cobb–Douglas
    σs = (0.0, 0.3, 1.0)

    processes(σ) = MarketData(
        regions   = [:R],
        products  = [:in1, :in2, :mid, :out],
        demand    = [DemandCurve(product = :out, region = :R, p0 = 500, q0 = 10, elasticity = 1.5)],
        supply    = [SupplyCurve(product = :in1, region = :R, p0 = 50, q0 = 10, elasticity = 0.8),
                     SupplyCurve(product = :in2, region = :R, p0 = 50, q0 = 3,  elasticity = 0.8)],
        processes = [Process(name = :stage1, vacost = 20,
                             inputs  = [ces(composite = 2.0, products = [:in1, :in2],
                                            shares = [0.6, 0.4], sigma = σ)],
                             outputs = [:mid => 1.0]),
                     Process(name = :stage2, vacost = 30,
                             inputs  = [leontief(product = :mid, coeff = 1.25)],
                             outputs = [:out => 1.0])])

    τ = 5.0
    δ = Dict((:A, :A) => 0.8, (:B, :A) => 0.2, (:A, :B) => 0.3, (:B, :B) => 0.7)
    trade(σ; qA = 30) = MarketData(
        regions   = [:A, :B],
        products  = [:x],
        demand    = [DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 1.5),
                     DemandCurve(product = :x, region = :B, p0 = 100, q0 = 10, elasticity = 1.5)],
        supply    = [SupplyCurve(product = :x, region = :A, p0 = 50,  q0 = qA, elasticity = 1.0),
                     SupplyCurve(product = :x, region = :B, p0 = 150, q0 = 5,  elasticity = 1.0)],
        transport_costs = Dict((:x, :A, :B) => τ, (:x, :B, :A) => τ),
        armington = [Armington(product = :x, sigma = σ, shares = δ)])
    home(res, r) = getq(res.production, r, :x) - flow(res, :x, r, r == :A ? :B : :A)

    @testset "process nest at σ = $σ: cost-minimising mix and zero profit" for σ in σs
        res = solve_market(processes(σ))
        @test solved(res)
        p(x) = price(res, :R, x)
        q1, q2 = getq(res.production, :R, :in1), getq(res.production, :R, :in2)
        @test q1 / q2 ≈ (0.6 / 0.4) * (p(:in2) / p(:in1))^σ rtol = 1e-5
        @test p(:mid) ≈ 2.0 * index([0.6, 0.4], [p(:in1), p(:in2)], σ) + 20 rtol = 1e-5
    end

    @testset "Armington at σ = $σ: variety mix and price index" for σ in σs
        res = solve_market(trade(σ))
        @test solved(res)
        for (r, o) in ((:A, :B), (:B, :A))
            qd, qi = pprice(res, r, :x), pprice(res, o, :x) + τ
            @test home(res, r) / flow(res, :x, o, r) ≈
                  (δ[(r, r)] / δ[(o, r)]) * (qi / qd)^σ rtol = 1e-5
            @test price(res, r, :x) ≈ index([δ[(r, r)], δ[(o, r)]], [qd, qi], σ) rtol = 1e-5
        end
    end

    @testset "demand nest at σ = $σ: the CES demand condition" for σ in σs
        nest = DemandNest(sigma = σ, products = [:pellets, :firewood], shares = [0.6, 0.4])
        res = solve_market(MarketData(
            regions = [:A], products = [:pellets, :firewood],
            demand  = [DemandCurve(product = nest, region = :A, p0 = 100, q0 = 20, elasticity = 0.5)],
            supply  = [SupplyCurve(product = :pellets,  region = :A, p0 = 100, q0 = 12, elasticity = 0.8),
                       SupplyCurve(product = :firewood, region = :A, p0 = 100, q0 = 8,
                                  elasticity = 0.8, shift = 0.5)]))
        @test solved(res)
        qp, qf = getq(res.consumption, :A, :pellets), getq(res.consumption, :A, :firewood)
        pp, pf = price(res, :A, :pellets), price(res, :A, :firewood)
        @test qp / qf ≈ (0.6 / 0.4) * (pf / pp)^σ rtol = 1e-5
        # Cobb–Douglas spends a fixed share on each member, whatever the prices
        isone(σ) && @test pp * qp / (pp * qp + pf * qf) ≈ 0.6 rtol = 1e-5
    end

    @testset "the branches are the limits of the CES" begin
        # σ = 0 and σ = 1 are written with their own formulas: an elasticity
        # just beside each must give nearly the same equilibrium. Beside 0 the
        # CES formula is out of floating-point range below σ ≈ 0.005, so the
        # neighbour is 0.01 and the tolerance follows (the variety mix moves
        # with the price ratio to the power σ)
        for (σ, ε, tol) in ((0.0, 1e-2, 2e-2), (1.0, 1e-4, 1e-3), (1.0, -1e-4, 1e-3))
            at, near = solve_market(trade(σ)), solve_market(trade(σ + ε))
            for r in (:A, :B)
                @test price(near, r, :x) ≈ price(at, r, :x) rtol = tol
                @test pprice(near, r, :x) ≈ pprice(at, r, :x) rtol = tol
                @test home(near, r) ≈ home(at, r) rtol = tol
            end
            at, near = solve_market(processes(σ)), solve_market(processes(σ + ε))
            for x in (:in1, :in2, :mid)
                @test price(near, :R, x) ≈ price(at, :R, x) rtol = tol
            end
        end
    end

    @testset "below 1 the varieties are complements" begin
        # more supply in A makes A consume more. Its users then want more of
        # B's variety too if the two are complements, and less if they are
        # substitutes
        imports(σ, qA) = flow(solve_market(trade(σ; qA)), :x, :B, :A)
        @test imports(0.3, 39) > imports(0.3, 30)
        @test imports(3.0, 39) < imports(3.0, 30)
    end

    @testset "nested groups: Cobb–Douglas inside, Leontief outside" begin
        regions = [:EU, :NA, :AS]
        shares  = Dict((:EU, :EU) => 0.5,  (:NA, :EU) => 0.3,  (:AS, :EU) => 0.2,
                       (:NA, :NA) => 0.6,  (:EU, :NA) => 0.25, (:AS, :NA) => 0.15,
                       (:AS, :AS) => 0.55, (:EU, :AS) => 0.2,  (:NA, :AS) => 0.25)
        res = solve_market(MarketData(
            regions  = regions, products = [:x],
            demand   = [DemandCurve(product = :x, region = r, p0 = 100, q0 = 10, elasticity = 1.5)
                        for r in regions],
            supply   = [SupplyCurve(product = :x, region = :EU, p0 = 60, q0 = 20, elasticity = 1.0),
                        SupplyCurve(product = :x, region = :NA, p0 = 55, q0 = 25, elasticity = 1.0),
                        SupplyCurve(product = :x, region = :AS, p0 = 70, q0 = 15, elasticity = 1.0)],
            transport_costs = Dict((:x, o, r) => 6.0 for o in regions, r in regions if o != r),
            armington = [Armington(product = :x, sigma = 0, shares = shares,
                                   nests = [OriginNest(sigma = 1, origins = [:EU, :NA])])]))
        @test solved(res)
        for r in regions
            q(o)  = pprice(res, o, :x) + (o == r ? 0.0 : 6.0)
            group = shares[(:EU, r)] + shares[(:NA, r)]
            Pg = index([shares[(:EU, r)], shares[(:NA, r)]] ./ group, [q(:EU), q(:NA)], 1.0)
            @test price(res, r, :x) ≈ group * Pg + shares[(:AS, r)] * q(:AS) rtol = 1e-5
        end
    end

    @testset "only what is raised to a power is floored" begin
        lb(m, name) = JuMP.lower_bound(JuMP.variable_by_name(m, name))
        @test lb(jump_model(trade(0.0)), "X[x,B,A]") == 0.0
        @test lb(jump_model(trade(1.0)), "X[x,B,A]") > 0.0
        @test lb(jump_model(processes(0.0)), "x[R,stage1,in1]") == 0.0
        @test lb(jump_model(processes(1.0)), "x[R,stage1,in1]") > 0.0
    end
end

@testset "Names may be strings or symbols" begin
    # the documented API is strings; symbols are accepted everywhere and
    # converted on the way in, so the two spellings must build the same economy
    build(n) = MarketData(
        regions   = [n("A"), n("B")],
        products  = [n("logs"), n("board")],
        supply    = [SupplyCurve(product = n("logs"), region = n("A"), p0 = 40, q0 = 35,
                                elasticity = 0.6),
                     SupplyCurve(product = n("logs"), region = n("B"), p0 = 90, q0 = 35,
                                elasticity = 0.6)],
        demand    = [DemandCurve(product = n("board"), region = n("A"), p0 = 300, q0 = 10,
                                elasticity = 1.2),
                     DemandCurve(product = n("board"), region = n("B"), p0 = 300, q0 = 20,
                                elasticity = 1.2)],
        processes = [Process(name = n("mill"), vacost = 20, regions = [n("A"), n("B")],
                             inputs = [leontief(product = n("logs"), coeff = 2)],
                             outputs = [n("board") => 1.0])],
        limits    = [@limit("mill A", activity(region = n("A"), process = n("mill")) <= 40.0),
                     @limit("port", exports(region = n("A"), product = n("board")) <= 50.0)],
        taxes     = [Tax(name = n("duty"), on = n("trade"), product = n("logs"),
                         from = n("A"), rate = 1.5)],
        transport_costs = Dict((n("logs"), n("A"), n("B")) => 5.0,
                               (n("logs"), n("B"), n("A")) => 5.0,
                               (n("board"), n("A"), n("B")) => 7.0,
                               (n("board"), n("B"), n("A")) => 7.0))

    strung, symbolic = build(identity), build(Symbol)
    # the stored economy is the same either way, down to the field types
    for f in fieldnames(MarketData)      # structs holding arrays compare by
        @test repr(getfield(strung, f)) == repr(getfield(symbolic, f))   # identity
    end
    @test strung.regions isa Vector{String} && symbolic.regions isa Vector{String}
    @test symbolic.processes[1].name === "mill" && symbolic.taxes[1].on == "trade"
    @test symbolic.limits[1].expression.terms[1].first.selection ==
          (region = "A", process = "mill")
    # a field read back compares to the string that built it, which is the whole
    # reason names are not stored as symbols
    @test symbolic.supply[1].region == "A" && symbolic.supply[1].product == "logs"

    # the nested types take either spelling too
    arm(n) = Armington(product = n("board"), sigma = 4,
                       shares = Dict((n("A"), n("B")) => 0.7, (n("B"), n("B")) => 0.3),
                       destination = n("B"),
                       nests = [OriginNest(sigma = 8, origins = [n("A"), n("B")])])
    @test repr(arm(identity)) == repr(arm(Symbol))
    @test repr(ces(composite = 1.3, products = ["a", "b"], sigma = 2)) ==
          repr(ces(composite = 1.3, products = [:a, :b], sigma = 2))
    @test repr(DemandNest(sigma = 2, products = ["a", "b"])) ==
          repr(DemandNest(sigma = 2, products = [:a, :b]))
    @test disposal(product = :chips).name == "dispose_chips"
    @test armington_shares(flows = Dict((:A, :B) => 1.0), prices = Dict((:A, :B) => 2.0),
                           sigma = 3) == Dict(("A", "B") => 1.0)

    # and the two spellings solve to the same numbers, reported with strings
    a, b = solve_market(strung), solve_market(symbolic)
    @test a.objective ≈ b.objective rtol = 1e-10
    @test a.prices.region == b.prices.region
    @test eltype(a.prices.region) == String && eltype(a.trade.from) == String
    @test price(a, "B", "board") ≈ price(b, :B, :board) rtol = 1e-10
    @test nrow(for_region(a, "A").prices) == nrow(for_region(a, :A).prices)
    @test solve_market(strung; on_failure = :warn).objective ≈ a.objective rtol = 1e-10

    # every one is the colon, and "all" is a name like any other
    @test Tax(name = "t", on = "supply", rate = 1.0).product isa Colon
    @test Process(name = "p", inputs = [leontief(product = "x", coeff = 1)],
                  outputs = ["y" => 1.0], vacost = 1, regions = :).regions isa Colon
    @test Process(name = "p", inputs = [leontief(product = "x", coeff = 1)],
                  outputs = ["y" => 1.0], vacost = 1, regions = :all).regions == "all"
    @test Markets.matches(Markets.selector(:), "anything")
    @test !Markets.matches(Markets.selector("all"), "anything")
    # anything that is not a name is refused where a name is wanted
    @test_throws ArgumentError SupplyCurve(product = 7, region = "A", p0 = 1, q0 = 1,
                                          elasticity = 1)
    @test_throws ArgumentError Tax(name = "t", on = "nowhere", rate = 1.0)
    @test_throws ArgumentError solve_market(strung; on_failure = "shrug")
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
    milk_made  = sum(skipmissing(res.production.quantity[res.production.product .== "milk"]))
    milk_drunk = sum(skipmissing(res.consumption.quantity[res.consumption.product .== "milk"]))
    @test milk_made ≈ milk_drunk + 10 * getq(res.production, :Pasture, :cheese) rtol = 1e-6
    @test 0.3 < milk_drunk / milk_made < 0.4
    # a thirstier City makes cheese dearer everywhere, although no taste for
    # cheese changed, and takes milk away from the dairy
    for r in (:Pasture, :City)
        @test only(ex.thirst.prices[(ex.thirst.prices.region .== S(r)) .&
                                    (ex.thirst.prices.product .== "cheese"), :price]) >
              price(res, r, :cheese)
    end
    @test getq(ex.thirst.production, :Pasture, :cheese) <
          getq(res.production, :Pasture, :cheese)
end

@testset "Forest example" begin
    # the documented example is three runnable scripts: run them in order —
    # the last reads what the first two save — and check that the equilibria
    # they report satisfy the conditions they claim
    run(file) = (m = Module(); Base.include(m, joinpath(@__DIR__, "..", "examples", "forest", file)); m)
    swr(r, result) = only(result.prices[(result.prices.region .== S(r)) .&
                                        (result.prices.product .== "swr"), :price])
    drop(r, ex) = 1 - swr(r, ex.storm) / swr(r, ex.res)

    @testset "writing the model" begin
        ex = run("forest_market.jl")
        d, res, france = ex.example_market, ex.res, ex.france
        @test solved(res) && solved(ex.storm) && solved(ex.gas)
        @test nrow(res.prices) == length(d.regions) * length(d.products)
        @test all(skipmissing(res.prices.price) .> 0)
        # papermill is Leontief: it buys the pulp composite and sells its own
        # variety of paper, so zero profit ties the producer price to the user one
        for r in d.regions
            @test pprice(res, r, :paper) ≈ 1.1 * price(res, r, :pulp) + 200 rtol = 1e-4
        end
        # no pulp mill in SEF and GEF, hence no pulp of their own to price
        @test ismissing(pprice(res, :SEF, :pulp)) && ismissing(pprice(res, :GEF, :pulp))
        # the French regions are near-perfect substitutes: wood moves between
        # them both ways, but far less of it than leaves France
        intra = res.trade[in.(res.trade.from, Ref(france)) .& in.(res.trade.to, Ref(france)), :]
        @test any(flow(res, row.product, row.to, row.from) > 1e-3 for row in eachrow(intra))
        @test sum(intra.quantity) < sum(ex.french_trade.crossing_the_border)
        # world trade balances: total exports = total imports for each product
        for g in groupby(res.net_trade, :product)
            @test sum(g.net) ≈ 0 atol = 1e-5
        end
        # the storm reaches France hardest, then the EU, then the far blocs
        @test minimum(drop.(france, Ref(ex))) > drop("EU", ex) >
              maximum(drop.(["NA", "AS", "RW"], Ref(ex)))
        @test all(drop(r, ex) > 0 for r in d.regions)
        # the shared stands bound every combination of the primary products:
        # slack in a normal year, binding once a third of France is set aside
        # every ceiling is a limit: the stands, the pellet plants and the Landes roads
        @test nrow(res.limits) == 7 * length(d.regions) + length(d.regions) + 1
        @test all(iszero, res.limits.rent)
        bound = ex.aside.limits[ex.aside.limits.rent .> 0, :]
        @test sort(bound.limit) == sort(["$r: swr + hwr + smallwood" for r in france])
        @test all(bound.value .≈ bound.upper)
        for r in france
            @test sum(getq(ex.aside.production, r, p) for p in ex.primary) <=
                  ex.available(r, ex.primary, ex.protected) + 1e-6
        end
        # uncalibrated, the reference run misses the one price given per product
        @test any(abs(price(res, r, :plywood) / ex.reference["plywood"][1] - 1) > 0.05
                  for r in d.regions)
    end

    @testset "calibrating with partial data" begin
        ex = run("forest_calibration.jl")
        res = ex.res
        @test solved(res) && solved(ex.storm) && solved(ex.gas)
        # the imputed flows add up to what customs and surveys report
        @test minimum(values(ex.flows)) >= 0
        # built from rounded, unbalanced data, the model lands near every
        # observation: within a per cent, a few per cent for customs flows
        largest(name) = last(ex.relative_gaps(last(only(c for c in ex.comparison
                                                         if first(c) == name))...))
        for name in ("producer prices", "production", "consumption")
            @test largest(name) < 0.01
        end
        @test largest("customs trade") < 0.05
        @test 0 < largest("production")               # near, not on
        # the physical account closes: what each region received of a product
        # is what its households and its mills used, physically
        for c in eachrow(res.composites)
            used = sum(res.consumption[(res.consumption.region .== c.region) .&
                                       (res.consumption.product .== c.product), :physical];
                       init = 0.0) +
                   sum(res.inputs[(res.inputs.region .== c.region) .&
                                  (res.inputs.product .== c.product), :physical]; init = 0.0)
            @test used ≈ c.received rtol = 1e-5
        end
        # the carbon accounts close, and the storm moves the content of pellets
        acc = ex.accounts
        @test sum(acc.balances.extracted) ≈ sum(acc.balances.consumed) +
                                            sum(acc.processes.lost) rtol = 1e-10
        pellets(a) = only(a.intensities[(a.intensities.region .== "SWF") .&
                                        (a.intensities.product .== "pellets"), :produced])
        @test pellets(ex.storm_accounts) > pellets(acc)
        # zero profit at each region's own cost
        papermill = only(p for p in ex.processes if p.name == "papermill")
        for r in papermill.regions
            @test pprice(res, r, :paper) ≈ 1.1 * price(res, r, :pulp) +
                                           Markets.vacost_of(papermill, r) rtol = 1e-5
        end
    end

    @testset "complete data, and the comparison" begin
        ex = run("forest_full_data.jl")
        res = ex.res
        @test solved(res)
        # consistent data come back to the precision of the solver
        for (model, observed) in (
                (ex.modelled(res.prices, :producer_price, :region, :product), ex.πˢ),
                (ex.modelled(res.production, :quantity, :region, :product), ex.made),
                (ex.modelled(res.consumption, :quantity, :region, :product), ex.consumed),
                (ex.modelled(res.trade, :quantity, :product, :from, :to), ex.shipped))
            @test last(ex.relative_gaps(model, observed)) < 1e-5
        end
        # partial data give the scenarios of complete data; uncalibrated ones
        # deepen the storm in the Landes
        row(case) = ex.comparison[startswith.(ex.comparison.case, case), :]
        @test abs(only(row("partial").storm_SWF) - only(row("complete").storm_SWF)) < 0.5
        @test only(row("uncalibrated").storm_SWF) < only(row("complete").storm_SWF) - 3
    end
end

@testset "Ports of published models" begin
    # a page may include a file of its own, as Documenter's modules allow
    function run(dir, file)
        m = Module()
        Core.eval(m, :(include(path) = Base.include($m, path)))
        Base.include(m, joinpath(@__DIR__, "..", "examples", dir, file))
        return m
    end

    @testset "Takayama and Judge (GAMS spatequ)" begin
        ex = run("spatequ", "spatequ.jl")
        @test solved(ex.result)
        # every quantity on its line and every region in balance ...
        @test maximum(abs, Matrix(ex.checks[:, [:demand_gap, :supply_gap, :balance_gap]])) < 1e-6
        # ... no price gap wider than the cost of closing it, and none narrower
        # where goods move
        a = ex.arbitrage
        @test all(a.price_gap .<= a.cost .+ 1e-6)
        @test all(abs.(a.price_gap[a.shipped .> 1e-6] .- a.cost[a.shipped .> 1e-6]) .< 1e-6)
        @test price(ex.result, :Reg2, :Com1) ≈ 9 rtol = 1e-6
    end

    @testset "Egyptian agriculture (GAMS egypt)" begin
        ex = run("egypt", "egypt.jl")
        @test solved(ex.result) && solved(ex.dry_result)
        # the port is the original with exact demand lines: same welfare, same
        # prices, same water rent when water binds
        @test ex.exact.objective - ex.meat_and_milk ≈ ex.result.objective rtol = 1e-7
        @test maximum(abs.(ex.comparison.difference ./ ex.comparison.transcription)) < 1e-5
        @test ex.exact_drier.objective - ex.meat_and_milk ≈ ex.dry_result.objective rtol = 1e-7
        @test ex.water_rent(ex.dry_result) ≈ ex.exact_drier.water_rent rtol = 1e-5
        @test ex.water_rent(ex.result) < 1e-6 < ex.water_rent(ex.dry_result)
        # the grid of the original is near the exact lines, not on them
        gap = maximum(abs(ex.original.price[c] / ex.exact.price[c] - 1) for c in ex.marketed)
        @test 1e-4 < gap < 0.05
    end
end

@testset "One Curve type, two sides" begin
    dem = DemandCurve(product = :x, region = :A, p0 = 100, q0 = 10, elasticity = 1.2)
    sup = SupplyCurve(product = :x, region = :A, alpha = -5, beta = 0.5)
    @test dem isa Markets.Curve{:demand} && sup isa Markets.Curve{:supply}
    @test !(dem isa SupplyCurve) && !(sup isa DemandCurve)
    @test typeof([dem]) == Vector{DemandCurve}
    @test sprint(show, dem) == "DemandCurve(product = \"x\", region = \"A\", p0 = 100.0, q0 = 10.0, elasticity = 1.2)"
    @test sprint(show, sup) == "SupplyCurve(product = \"x\", region = \"A\", alpha = -5.0, beta = 0.5)"
    # the same calibration, moved: the alias keeps its side
    @test DemandCurve(dem; shift = 2.0) isa DemandCurve
    @test Markets.curve_price(DemandCurve(dem; shift = 2.0), 20.0) ≈ 100
    @test Markets.curve_price(sup, 5.0) ≈ 20       # S = -5 + 0.5 p  →  p = (S + 5) / 0.5
    # the sign of the elasticity in the exponent is the whole difference
    @test Markets.curve_price(dem, 20.0) ≈ 100 * 2^(-1 / 1.2)
    @test Markets.curve_price(SupplyCurve(product = :x, region = :A, p0 = 100, q0 = 10,
                                          elasticity = 1.2), 20.0) ≈ 100 * 2^(1 / 1.2)
    # the integral of the inverse curve is what the objective holds
    @test Markets.curve_integral(sup, 5.0) ≈ (5^2 / 2 + 5 * 5) / 0.5
    @test Markets.curve_integral(dem, 20.0) ≈ 100 * 10^(1 / 1.2) / (1 - 1 / 1.2) * 20^(1 - 1 / 1.2)
    # a supply curve describes one product, never a group
    @test_throws ArgumentError SupplyCurve(product = DemandNest(sigma = 2, products = [:a, :b]),
                                           region = :A, p0 = 1, q0 = 1, elasticity = 1)
    @test_throws ArgumentError Markets.Curve{:sideways}(:x, :A, :isoelastic, 1.0, 1.0, 1.0, NaN, NaN, 1.0)
end

@testset "Corner cases of the groups, the balances and the bounds" begin
    two(; kw...) = MarketData(regions = [:A, :B], products = [:w],
        supply = [SupplyCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 1),
                  SupplyCurve(product = :w, region = :B, p0 = 12, q0 = 100, elasticity = 1)],
        demand = [DemandCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 1),
                  DemandCurve(product = :w, region = :B, p0 = 12, q0 = 100, elasticity = 1)],
        transport_costs = Dict((:w, :A, :B) => 1.0, (:w, :B, :A) => 1.0); kw...)

    @testset "every origin in one finite group is a group with shares" begin
        grouped = two(armington = [Armington(product = :w, sigma = 2,
                                             nests = [OriginNest(sigma = 5, origins = [:A, :B])])])
        flat    = two(armington = [Armington(product = :w, sigma = 5)])
        rg, rf = solve_market(grouped), solve_market(flat)
        @test solved(rg)
        # a group holding every origin is the flat CES at the group's elasticity
        @test rg.prices.price ≈ rf.prices.price rtol = 1e-6
        @test rg.trade.quantity ≈ rf.trade.quantity rtol = 1e-6
    end

    @testset "a group of one product keeps the floor of the curve over it" begin
        d = MarketData(regions = [:A], products = [:w],
            supply = [SupplyCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 1)],
            demand = [DemandCurve(product = DemandNest(sigma = Inf, products = [:w]), region = :A,
                                  p0 = 10, q0 = 100, elasticity = 0.5)])
        plain = MarketData(regions = d.regions, products = d.products, supply = d.supply,
            demand = [DemandCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 0.5)])
        floor_of(m) = JuMP.lower_bound(JuMP.variable_by_name(m, "D[A,w]"))
        @test floor_of(jump_model(d)) == floor_of(jump_model(plain)) > 0
        # and no composite or group constraint is written for it
        @test !occursin("demand_group", sprint(print, jump_model(d)))
    end

    @testset "a product nothing touches is allowed, and has no price" begin
        spare = MarketData(regions = [:A, :B], products = [:w, :v],
            supply = [SupplyCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 1)],
            demand = [DemandCurve(product = :w, region = :B, p0 = 12, q0 = 100, elasticity = 1)],
            transport_costs = Dict((:w, :A, :B) => 1.0))
        res = solve_market(spare)
        @test solved(res)
        @test all(ismissing, res.prices.price[res.prices.product .== "v"])
        @test all(!ismissing, res.prices.price[res.prices.product .== "w"])
        @test !occursin("balance[A,v]", sprint(print, res.model))
    end

    @testset "a cap below the floor is refused with a reason" begin
        d = MarketData(regions = [:A], products = [:w],
            supply = [SupplyCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 1)],
            demand = [DemandCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 0.5)],
            limits = [@limit "zero" consumption(region = :A, product = :w) <= 0])
        @test_throws ArgumentError solve_market(d)
        err = try solve_market(d) catch e; e end
        @test occursin("below the floor", err.msg)
    end

    @testset "two CES nests of one process get distinct names" begin
        d = MarketData(regions = [:A], products = [:a, :b, :c],
            supply = [SupplyCurve(product = p, region = :A, p0 = 10, q0 = 100, elasticity = 1) for p in [:a, :b]],
            demand = [DemandCurve(product = :c, region = :A, p0 = 60, q0 = 50, elasticity = 1)],
            processes = [Process(name = :q, outputs = [:c => 1.0], vacost = 1,
                                 inputs = [ces(composite = 1, products = [:a, :b], sigma = 2),
                                           ces(composite = 1, products = [:a, :b], sigma = 3)])])
        m = jump_model(d)
        names = [name(c) for c in JuMP.all_constraints(m; include_variable_in_set_constraints = false)]
        @test allunique(filter(!isempty, names))
        @test "input_nest[A,q,1]" in names && "input_nest[A,q,2]" in names
    end

    @testset "a single-product ces nest cannot carry a scale" begin
        @test_throws ArgumentError ces(composite = 1.0, products = [:a], sigma = 2, phi = 2.0)
        @test ces(composite = 1.0, products = [:a], sigma = 2).phi == 1.0
    end

    @testset "exports and imports with no region are every flow" begin
        d = two(limits = [@limit "world trade" exports(product = :w) <= 5])
        res = solve_market(d)
        @test only(res.limits.value) ≈ sum(res.trade.quantity) rtol = 1e-6
        @test Markets.scope(d.limits[1], d).variables == ["flow[w, A, B]", "flow[w, B, A]"]
        @test Markets.scope(@limit("in", imports(product = :w) <= 5), d).variables ==
              ["flow[w, A, B]", "flow[w, B, A]"]
    end

    @testset "validate is the whole verdict" begin
        # what used to surface only while the model was built
        @test_throws ArgumentError validate(two(armington = [Armington(product = :w, sigma = 2),
                                                             Armington(product = :w, sigma = 3)]))
        @test_throws ArgumentError validate(two(armington = [Armington(product = :w, sigma = 2,
                                                                       destination = :A)]))
        @test_throws ArgumentError validate(two(limits = [@limit "none" activity(process = :mill) <= 1]))
        @test_throws ArgumentError validate(MarketData(regions = [:A, :B], products = [:w],
            supply = two().supply[1:1], demand = two().demand[2:2],
            transport_costs = Dict((:w, :A, :B) => 1.0),
            limits = [@limit "none" supply(region = :B, product = :w) <= 5]))
        # a tax that selects nothing it could fall on
        @test_throws ArgumentError validate(two(taxes = [Tax(on = :supply, product = :w, region = :C, rate = 1)]))
        err = try validate(MarketData(regions = [:A, :B], products = [:w, :v],
            supply = two().supply, demand = two().demand, transport_costs = two().transport_costs,
            taxes = [Tax(name = :t, on = :supply, product = :v, rate = 1)])) catch e; e end
        @test err isa ArgumentError && occursin("selects no supply curve", err.msg)
        @test_throws ArgumentError validate(two(taxes = [Tax(on = :trade, product = :w, from = :A, to = :A, rate = 1)]))
        @test_throws ArgumentError validate(two(taxes = [Tax(on = :process, process = :mill, rate = 1)]))
        # a free process is held back by a ceiling, not by a floor
        free(limit) = MarketData(regions = [:A], products = [:w],
            demand = [DemandCurve(product = :w, region = :A, p0 = 12, q0 = 100, elasticity = 1)],
            processes = [Process(name = :magic, inputs = Nest[], outputs = [:w => 1.0])],
            limits = [limit])
        @test_throws ArgumentError validate(free(@limit "floor" activity(process = :magic) >= 1))
        @test validate(free(@limit "cap" activity(process = :magic) <= 1)) === nothing
        @test validate(free(@limit "cap" -activity(process = :magic) >= -1)) === nothing
        # members of one demand group carry one ad valorem rate
        grouped = MarketData(regions = [:A], products = [:a, :b],
            supply = [SupplyCurve(product = p, region = :A, p0 = 10, q0 = 100, elasticity = 1) for p in [:a, :b]],
            demand = [DemandCurve(product = DemandNest(sigma = 2, products = [:a, :b]), region = :A,
                                  p0 = 10, q0 = 200, elasticity = 1)],
            taxes = [Tax(on = :demand, product = :a, rate = 0.2, relative = true)])
        @test_throws ArgumentError validate(grouped)
    end
end

@testset "jump_model refuses the keywords of the solve" begin
    d = MarketData(regions = [:A], products = [:w],
                   supply = [SupplyCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 1)],
                   demand = [DemandCurve(product = :w, region = :A, p0 = 10, q0 = 100, elasticity = 1)])
    err = try jump_model(d; tol = 1e-3) catch e; e end
    @test err isa ArgumentError && occursin("tol", err.msg) && occursin("qfloor", err.msg)
    @test jump_model(d; qfloor = 1e-9) isa JuMP.Model
end

@testset "validate warns about Armington shares the model will drop" begin
    three(; shares, routes = [(o, r) for o in [:A, :B, :C] for r in [:A, :B, :C] if o != r]) =
        MarketData(regions = [:A, :B, :C], products = [:w],
                   supply = [SupplyCurve(product = :w, region = r, p0 = 10, q0 = 100, elasticity = 1)
                             for r in [:A, :B]],
                   demand = [DemandCurve(product = :w, region = r, p0 = 12, q0 = 60, elasticity = 1)
                             for r in [:A, :B, :C]],
                   transport_costs = Dict((:w, o, r) => 1.0 for (o, r) in routes),
                   armington = [Armington(product = :w, sigma = 2, shares = shares)])
    every = Dict((o, r) => 1.0 for o in [:A, :B, :C] for r in [:A, :B, :C])
    # C makes no w, so its shares fall on a variety that does not exist
    @test_logs (:warn, r"C → A \(C produces no w\).*C → C") validate(three(shares = every))
    # a pair with no route is said too
    @test_logs (:warn, r"A → C \(no route from A to C\)") validate(
        three(shares = Dict(k => v for (k, v) in every if k[1] != :C),
              routes = [(:A, :B), (:B, :A), (:B, :C), (:C, :A), (:C, :B)]))
    # shares on the varieties the model has give no warning, and the economy solves
    fine = three(shares = Dict(k => v for (k, v) in every if k[1] != :C))
    @test_logs validate(fine)
    @test solve_market(fine).status in (Markets.MOI.OPTIMAL, Markets.MOI.LOCALLY_SOLVED)
end

end
