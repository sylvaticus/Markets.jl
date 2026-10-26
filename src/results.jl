# ----------------------------------------------------------------------------
#  Results & retrieval
# ----------------------------------------------------------------------------

"""
$(TYPEDEF)

The outcome of [`solve_market`](@ref): the solved model and the equilibrium
quantities and prices, as tidy `DataFrame`s.

Use [`for_region`](@ref) to slice every table for one region at once.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct Results
    "The economy that was solved"
    data::MarketData
    "The underlying JuMP model, solved"
    model::Model
    """
    What the solver made of it. Anything but `OPTIMAL` or `LOCALLY_SOLVED`
    means the tables below do not describe an equilibrium, whatever they look
    like; [`solve_market`](@ref) throws on those unless told otherwise
    """
    status::MOI.TerminationStatusCode
    "The solver's own account of how it stopped, such as Ipopt's `Solve_Succeeded`"
    raw_status::String
    """
    Seconds spent solving, as the solver reports them, or as measured around
    the solve for a solver that does not. Summed over all `solves`
    """
    solve_time::Float64
    """
    The iterations the solver took, summed over all `solves`: interior-point
    iterations for Ipopt, simplex ones for a simplex solver. `missing` for a
    solver that reports neither
    """
    iterations::Union{Int,Missing}
    """
    How many times the model was solved: 1, unless an ad valorem tax on trade
    had to be found by iteration
    """
    solves::Int = 1
    """
    The value of the maximised objective: net social surplus, **less** any tax
    revenue, which the program counts as a cost rather than a transfer
    """
    objective::Float64
    """
    The objective with the tax revenue added back, which is the surplus the
    economy actually enjoys and the figure to compare across policy scenarios
    """
    welfare::Float64
    "Primary supply + manufactured output: `region`, `product`, `quantity`"
    production::DataFrame
    """
    Primary supply alone: `region`, `product`, `quantity`, the `capacity`
    bounding it and the scarcity `rent` a binding capacity earns per unit
    """
    supply::DataFrame
    """
    Final demand: `region`, `product`, `quantity`, in the units the demand
    curves and prices use, and `physical`, the same in the physical units of
    the flows. The two differ only for an [`Armington`](@ref) product, whose
    users buy a composite: see `composites`
    """
    consumption::DataFrame
    """
    The composite of each [`Armington`](@ref) product in each region that uses
    it, against the physical quantity behind it: `region`, `product`,
    `received` (the sum of the varieties delivered there, its own included),
    `composite` (what its users draw on, final demand and processes together)
    and `ratio` = `received` / `composite`. A CES composite equals the sum of
    what it is made of only at the mix it was calibrated at; elsewhere `ratio`
    is the physical quantity a unit of composite stands for, and every use of
    the composite in the region is credited with it in `physical` columns.
    Value is exact on the composite side: `composite` times the price users
    pay is what was spent on the varieties, delivered
    """
    composites::DataFrame
    """
    What each process draws of each input: `region`, `process`, `product`,
    `composite` (in the units the balances and prices use) and `physical` (the
    same in the physical units of the flows)
    """
    inputs::DataFrame
    """
    Each input nest of each process against its recipe: `region`, `process`,
    `inputs` (the products of the nest), `recipe` (the nest's requirement per
    unit of activity times the activity), `physical` (the physical quantity of
    its inputs, summed) and `gap` = `physical` − `recipe`. A recipe is met in
    composite units; the gap is the physical quantity the aggregation of
    imperfect substitutes uses up (positive) or makes up (negative) on the
    way, whether in the region's Armington composites or in a CES nest
    itself. It closes a physical account of the economy. The sum means
    something only where the products of a nest share a physical unit
    """
    recipes::DataFrame
    """
    Every solved quantity, unrounded and unfiltered, as dictionaries keyed as in
    the model: `supply` and `consumption` by `(region, product)`, `activity` by
    `(region, process)`, `flows` by `(product, from, to)` — a region's own
    deliveries of an Armington product included — `drawn` by
    `(region, process, nest, product)`, and `composite` and `ratio` by
    `(region, product)`: each region's composite of an Armington product, and
    the physical content of a unit of it. The tables above are built from these
    with small quantities left out; accounts that must close exactly, such as
    [`material_accounts`](@ref), start from them instead
    """
    solution::NamedTuple
    "Bilateral trade flows, positive ones only: `product`, `from`, `to`, `quantity`"
    trade::DataFrame
    "Trade by region: `region`, `product`, `exports`, `imports`, `net` (= exports − imports)"
    net_trade::DataFrame
    """
    Equilibrium prices, the duals of the balances: `region`, `product`,
    `price` (what local users pay) and `producer_price` (what local producers
    get). The two differ only for an [`Armington`](@ref) product, where users
    buy a composite of the regional varieties
    """
    prices::DataFrame
    """
    Process activity levels: `region`, `process`, `level`, the `capacity`
    bounding it and the quasi-`rent` a binding capacity earns per unit of
    activity
    """
    activity::DataFrame
    """
    Every [`Limit`](@ref) and [`TradeLimit`](@ref): `limit`, its `lower` and
    `upper` bounds (`-Inf` and `Inf` where it has none), the `value` its
    expression takes, and its `rent` — the gain in the objective from relaxing
    the bound that binds by one unit: positive for an upper bound, like a
    congestion rent or a quota price, negative for a lower one, the subsidy
    that would get agents to meet it unaided. Empty when there are none
    """
    limits::DataFrame
    """
    What each [`Tax`](@ref) raises: `tax`, what it falls `on`, and its
    `revenue` (negative for a subsidy). Adding the total back to the objective
    recovers welfare, which a taxed objective otherwise understates by exactly
    this much
    """
    taxes::DataFrame
end

# Build the result tables from the solved model.  Internal: users get a
# `Results` from `solve_market`.
"""
$(TYPEDSIGNATURES)

Zero when `x` is within solver noise of it, and `x` otherwise.  Used on rents
so that numerical dust reads as no rent while a genuinely wrong sign survives
to be seen.
"""
denoise(x, atol = 1e-9) = abs(x) < atol ? zero(x) : x

"""
$(TYPEDSIGNATURES)

What the solver reports of `attr`, or `fallback` where it reports nothing:
not every solver supports every attribute, and asking one that does not
throws.
"""
function solver_attribute(m::Model, attr, fallback)
    try
        return MOI.get(m, attr)
    catch
        return fallback
    end
end

function build_results(d, m, D, S, z, T, X, balance, origin_balance, limits_built,
                       tariff, tol, wall_time;
                       A = Dict{Tuple{String,String},VariableRef}(),
                       drawn = Dict{Tuple{String,String,Int,String},Any}())
    val(x) = value(x)
    # bilateral flows: the trade variables of the homogeneous products and the
    # cross-border varieties of the Armington ones
    flows = Dict{Tuple{String,String,String},Any}(T)
    for (key, v) in X
        key[2] == key[3] && continue      # domestic variety, not trade
        flows[key] = v
    end

    # production = primary supply + summed process outputs
    prod = Dict{Tuple{String,String},Float64}()
    for (k, v) in S
        prod[k] = get(prod, k, 0.0) + val(v)
    end
    for proc in d.processes, r in regions_of(proc, d)
        zk = val(z[(r, proc.name)])
        for (p, yield) in proc.outputs
            prod[(r, p)] = get(prod, (r, p), 0.0) + yield * zk
        end
    end
    production = DataFrame(region = String[], product = String[], quantity = Float64[])
    for r in d.regions, p in d.products
        q = get(prod, (r, p), 0.0)
        q > tol && push!(production, (r, p, q))
    end

    # primary supply on its own, with the ceiling and what a binding one earns
    supply = DataFrame(region = String[], product = String[], quantity = Float64[],
                       capacity = Float64[], rent = Float64[])
    for spec in d.supply
        v = S[(spec.region, spec.product)]
        q = val(v)
        # the rent is the gap the ceiling opens between the producer price and
        # the taxed marginal cost on the curve, i.e. the multiplier of the
        # ceiling itself. Taking it from the bound rather than recomputing it
        # from the curve keeps out the tax wedge and, for an Armington product,
        # the gap between what users pay and what producers get. (For a `<=`
        # constraint JuMP's dual is non-positive whatever the sense.)
        rent = has_duals(m) && isfinite(spec.capacity) ?
               denoise(-dual(UpperBoundRef(v))) : 0.0
        push!(supply, (spec.region, spec.product, q, spec.capacity, rent))
    end
    sort!(supply, [:region, :product])

    # physical accounts: what a unit of each region's composite stands for, in
    # the physical units of the flows that make it up
    composites = DataFrame(region = String[], product = String[], received = Float64[],
                           composite = Float64[], ratio = Float64[])
    ratio = Dict{Tuple{String,String},Float64}()
    for ((p, r), a) in A
        composite = val(a)
        received = sum(val(v) for ((q, _, d), v) in X if q == p && d == r; init = 0.0)
        composite > tol || continue
        ratio[(r, p)] = received / composite
        push!(composites, (r, p, received, composite, ratio[(r, p)]))
    end
    sort!(composites, [:region, :product])
    physical(r, p, q) = q * get(ratio, (r, p), 1.0)

    consumption = DataFrame(region = String[], product = String[], quantity = Float64[],
                            physical = Float64[])
    for ((r, p), v) in D
        q = val(v)
        q > tol && push!(consumption, (r, p, q, physical(r, p, q)))   # the same rule as the others
    end
    sort!(consumption, [:region, :product])

    inputs = DataFrame(region = String[], process = String[], product = String[],
                       composite = Float64[], physical = Float64[])
    taken = Dict{Tuple{String,String,String},Float64}()
    for ((r, proc, _, p), e) in drawn
        taken[(r, proc, p)] = get(taken, (r, proc, p), 0.0) + val(e)
    end
    for ((r, proc, p), q) in taken
        q > tol && push!(inputs, (r, proc, p, q, physical(r, p, q)))
    end
    sort!(inputs, [:region, :process, :product])

    recipes = DataFrame(region = String[], process = String[], inputs = String[],
                        recipe = Float64[], physical = Float64[], gap = Float64[])
    for proc in d.processes, r in regions_of(proc, d)
        level = val(z[(r, proc.name)])
        level > tol || continue
        for (k, nest) in enumerate(proc.inputs)
            recipe = nest.composite * level
            phys = sum(physical(r, p, val(drawn[(r, proc.name, k, p)])) for p in nest.products)
            push!(recipes, (r, proc.name, join(nest.products, ", "), recipe, phys, phys - recipe))
        end
    end
    sort!(recipes, [:region, :process])

    trade = DataFrame(product = String[], from = String[], to = String[], quantity = Float64[])
    for ((p, from, to), v) in flows
        q = val(v)
        q > tol && push!(trade, (p, from, to, q))
    end
    sort!(trade, [:product, :from, :to])

    net = DataFrame(region = String[], product = String[],
                    exports = Float64[], imports = Float64[], net = Float64[])
    for r in d.regions, p in d.products
        ex = sum((val(flows[(p, r, o)]) for o in d.regions if haskey(flows, (p, r, o))); init = 0.0)
        im = sum((val(flows[(p, o, r)]) for o in d.regions if haskey(flows, (p, o, r))); init = 0.0)
        (ex > tol || im > tol) && push!(net, (r, p, ex, im, ex - im))
    end

    # Prices are the duals of the balances, with their own sign: a by-product
    # nobody wants, which the balance forces someone to take, has a negative
    # price.  (`shadow_price` must not be used here: for an equality constraint
    # it reports the gain from relaxing it in whichever direction helps, which
    # discards that sign.)
    prices = DataFrame(region = String[], product = String[])
    user = Any[]
    producer = Any[]
    if has_duals(m)
        arm = Set(a.product for a in d.armington if isfinite(a.sigma))
        for r in d.regions, p in d.products
            push!(prices, (r, p))
            # what local users pay, `missing` where the region uses none of it
            push!(user, haskey(balance, (r, p)) ? dual(balance[(r, p)]) : missing)
            # what local producers get: the same, unless the product is an
            # Armington one, and `missing` where the region has no variety of
            # it to sell
            push!(producer, haskey(origin_balance, (r, p)) ? dual(origin_balance[(r, p)]) :
                            p in arm ? missing : user[end])
        end
        prices.price          = identity.(user)       # Float64 unless any is missing
        prices.producer_price = identity.(producer)
        sort!(prices, [:region, :product])
    else
        prices.price          = Float64[]
        prices.producer_price = Float64[]
    end

    activity = DataFrame(region = String[], process = String[], level = Float64[],
                         capacity = Float64[], rent = Float64[])
    for proc in d.processes, r in regions_of(proc, d)
        lvl = val(z[(r, proc.name)])
        lvl > tol || continue
        cap  = capacity_of(proc, r)
        # a capacity that binds earns a quasi-rent per unit of activity: what
        # the process would pay for one more unit of itself. The value is taken
        # as the solver gives it — only numerical dust around zero is swept
        # away, so that a wrong sign would show rather than be clipped to 0
        rent = has_duals(m) && isfinite(cap) ? denoise(reduced_cost(z[(r, proc.name)])) : 0.0
        push!(activity, (r, proc.name, lvl, cap, rent))
    end

    # the value of each limit's expression, and what relaxing the bound that
    # binds would be worth: in a maximisation, minus the dual of the constraint
    limits = DataFrame(limit = String[], lower = Float64[], upper = Float64[],
                       value = Float64[], rent = Float64[])
    for L in limits_built
        value = val(L.expression)
        # complementary slackness: a bound that is not reached earns nothing,
        # whatever numerical dust its dual carries
        reached(bound) = abs(value - bound) <= 1e-6 * max(1.0, abs(bound))
        rent = has_duals(m) ?
            sum((-dual(c) for c in L.constraints if reached(normalized_rhs(c))); init = 0.0) : 0.0
        push!(limits, (L.limit.name, L.limit.lower, L.limit.upper, value, denoise(rent)))
    end

    # what each tax raises, so that the transfer the objective treats as a
    # cost can be added back
    taxes = DataFrame(tax = String[], on = String[], revenue = Float64[])
    price_at(r, p) = (rows = prices[(prices.region .== r) .& (prices.product .== p), :];
                      isempty(rows) ? 0.0 : coalesce(only(rows.price), 0.0))
    producer_at(r, p) = (rows = prices[(prices.region .== r) .& (prices.product .== p), :];
                         isempty(rows) ? 0.0 : coalesce(only(rows.producer_price), 0.0))
    for t in d.taxes
        revenue = 0.0
        if t.on == "supply"
            for spec in d.supply
                selects(t.product, spec.product) && selects(t.region, spec.region) || continue
                q  = val(S[(spec.region, spec.product)])
                mc = spec.p0 * (spec.shift * spec.q0)^(-1 / spec.elasticity) *
                     q^(1 / spec.elasticity)        # the untaxed marginal cost
                revenue += (t.relative ? t.rate * mc : t.rate) * q
            end
        elseif t.on == "demand"
            for spec in d.demand, p in nest_products(spec.product)
                selects(t.product, p) && selects(t.region, spec.region) || continue
                q = val(D[(spec.region, p)])
                revenue += (t.relative ? t.rate * price_at(spec.region, p) : t.rate) * q
            end
        elseif t.on == "process"
            for proc in d.processes, r in regions_of(proc, d)
                selects(t.process, proc.name) && selects(t.region, r) || continue
                revenue += (t.relative ? t.rate * vacost_of(proc, r) : t.rate) * val(z[(r, proc.name)])
            end
        else
            for ((p, o, r), v) in flows
                selects(t.product, p) && selects(t.from, o) && selects(t.to, r) || continue
                revenue += (t.relative ? t.rate * producer_at(o, p) : t.rate) * val(v)
            end
        end
        push!(taxes, (t.name, t.on, revenue))
    end

    solution = (supply      = Dict(k => val(v) for (k, v) in S),
                consumption = Dict(k => val(v) for (k, v) in D),
                activity    = Dict(k => val(v) for (k, v) in z),
                flows       = Dict{Tuple{String,String,String},Float64}(
                                  k => val(v) for (k, v) in merge(Dict{Tuple{String,String,String},Any}(T), X)),
                drawn       = Dict(k => val(e) for (k, e) in drawn),
                composite   = Dict((r, p) => val(a) for ((p, r), a) in A),
                ratio       = ratio)

    objective = has_values(m) ? objective_value(m) : NaN
    iterations = something(solver_attribute(m, MOI.BarrierIterations(), nothing),
                           solver_attribute(m, MOI.SimplexIterations(), nothing),
                           missing)
    return Results(data = d, model = m, status = termination_status(m),
                   raw_status = raw_status(m),
                   solve_time = solver_attribute(m, MOI.SolveTimeSec(), wall_time),
                   iterations = iterations,
                   objective = objective, welfare = objective + sum(taxes.revenue; init = 0.0),
                   supply = supply, production = production,
                   consumption = consumption, composites = composites, inputs = inputs,
                   recipes = recipes, solution = solution, trade = trade, net_trade = net,
                   prices = prices, activity = activity, limits = limits, taxes = taxes)
end


"""
$(TYPEDSIGNATURES)

Slice every table of `res` for a single `region` (trade keeps both incoming and
outgoing flows).

Returns a `NamedTuple`
`(; production, supply, consumption, composites, inputs, recipes, prices,
net_trade, activity, trade)`.
The `limits` and `taxes` tables are not among them: a trade ceiling spans
routes and a tax may span regions, so neither belongs to one region.
"""
function for_region(res::Results, region::Name)
    region = string(region)
    f(df, col) = df[df[!, col] .== region, :]
    (production  = f(res.production, :region),
     supply      = f(res.supply, :region),
     consumption = f(res.consumption, :region),
     composites  = f(res.composites, :region),
     inputs      = f(res.inputs, :region),
     recipes     = f(res.recipes, :region),
     prices      = f(res.prices, :region),
     net_trade   = f(res.net_trade, :region),
     activity    = f(res.activity, :region),
     trade       = res.trade[(res.trade.from .== region) .| (res.trade.to .== region), :])
end
