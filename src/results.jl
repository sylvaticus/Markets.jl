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
    "Final demand: `region`, `product`, `quantity`"
    consumption::DataFrame
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
    Ceilings on trade: `limit`, its `capacity`, how much of it is `used`, and
    the congestion `rent` an extra unit of it would be worth. Empty when the
    economy has no [`TradeLimit`](@ref)
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

function build_results(d, m, D, S, z, T, X, balance, origin_balance, trade_limit,
                       tariff, tol)
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
        q = val(S[(spec.region, spec.product)])
        # the rent is the gap the ceiling opens between the price and what the
        # curve says the last unit cost to produce
        marginal = spec.p0 * (q / (spec.shift * spec.q0))^(1 / spec.elasticity)
        rent = has_duals(m) && isfinite(spec.capacity) ?
               denoise(dual(balance[(spec.region, spec.product)]) - marginal) : 0.0
        push!(supply, (spec.region, spec.product, q, spec.capacity, rent))
    end
    sort!(supply, [:region, :product])

    consumption = DataFrame(region = String[], product = String[], quantity = Float64[])
    for ((r, p), v) in D
        q = val(v)
        q > tol && push!(consumption, (r, p, q))   # the same rule as the others
    end
    sort!(consumption, [:region, :product])

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

    limits = DataFrame(limit = String[], capacity = Float64[], used = Float64[],
                       rent = Float64[])
    for limit in d.trade_limits
        c = trade_limit[limit.name]
        used = sum(get(limit.weights, k[1], 1.0) * val(v)
                   for (k, v) in flows if covers(limit, k...); init = 0.0)
        # `shadow_price` is the gain from relaxing the constraint by one unit,
        # which is exactly the congestion rent. (For an *equality* it is the
        # wrong function, as it discards the sign — see the prices below.)
        rent = has_duals(m) ? denoise(shadow_price(c)) : 0.0
        push!(limits, (limit.name, limit.capacity, used, rent))
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
                revenue += (t.relative ? t.rate * proc.vacost : t.rate) * val(z[(r, proc.name)])
            end
        else
            for ((p, o, r), v) in flows
                selects(t.product, p) && selects(t.from, o) && selects(t.to, r) || continue
                revenue += (t.relative ? t.rate * producer_at(o, p) : t.rate) * val(v)
            end
        end
        push!(taxes, (t.name, t.on, revenue))
    end

    objective = has_values(m) ? objective_value(m) : NaN
    return Results(data = d, model = m, status = termination_status(m),
                   objective = objective, welfare = objective + sum(taxes.revenue; init = 0.0),
                   supply = supply, production = production,
                   consumption = consumption, trade = trade, net_trade = net,
                   prices = prices, activity = activity, limits = limits, taxes = taxes)
end


"""
$(TYPEDSIGNATURES)

Slice every table of `res` for a single `region` (trade keeps both incoming and
outgoing flows).

Returns a `NamedTuple`
`(; production, supply, consumption, prices, net_trade, activity, trade)`.
The `limits` and `taxes` tables are not among them: a trade ceiling spans
routes and a tax may span regions, so neither belongs to one region.
"""
function for_region(res::Results, region::Name)
    region = string(region)
    f(df, col) = df[df[!, col] .== region, :]
    (production  = f(res.production, :region),
     supply      = f(res.supply, :region),
     consumption = f(res.consumption, :region),
     prices      = f(res.prices, :region),
     net_trade   = f(res.net_trade, :region),
     activity    = f(res.activity, :region),
     trade       = res.trade[(res.trade.from .== region) .| (res.trade.to .== region), :])
end
