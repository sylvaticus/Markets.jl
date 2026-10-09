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
    "Primary supply + manufactured output: `region`, `product`, `quantity`"
    production::DataFrame
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
    "Process activity levels: `region`, `process`, `level`"
    activity::DataFrame
end

# Build the result tables from the solved model.  Internal: users get a
# `Results` from `solve_market`.
function build_results(d, m, D, S, z, T, X, balance, origin_balance, tol)
    val(x) = value(x)
    # bilateral flows: the trade variables of the homogeneous products and the
    # cross-border varieties of the Armington ones
    flows = Dict{Tuple{Symbol,Symbol,Symbol},Any}(T)
    for (key, v) in X
        key[2] == key[3] && continue      # domestic variety, not trade
        flows[key] = v
    end

    # production = primary supply + summed process outputs
    prod = Dict{Tuple{Symbol,Symbol},Float64}()
    for (k, v) in S
        prod[k] = get(prod, k, 0.0) + val(v)
    end
    for proc in d.processes, r in regions_of(proc, d)
        zk = val(z[(r, proc.name)])
        for (p, yield) in proc.outputs
            prod[(r, p)] = get(prod, (r, p), 0.0) + yield * zk
        end
    end
    production = DataFrame(region = Symbol[], product = Symbol[], quantity = Float64[])
    for r in d.regions, p in d.products
        q = get(prod, (r, p), 0.0)
        q > tol && push!(production, (r, p, q))
    end

    consumption = DataFrame(region = Symbol[], product = Symbol[], quantity = Float64[])
    for ((r, p), v) in D
        push!(consumption, (r, p, val(v)))
    end
    sort!(consumption, [:region, :product])

    trade = DataFrame(product = Symbol[], from = Symbol[], to = Symbol[], quantity = Float64[])
    for ((p, from, to), v) in flows
        q = val(v)
        q > tol && push!(trade, (p, from, to, q))
    end
    sort!(trade, [:product, :from, :to])

    net = DataFrame(region = Symbol[], product = Symbol[],
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
    prices = DataFrame(region = Symbol[], product = Symbol[])
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

    activity = DataFrame(region = Symbol[], process = Symbol[], level = Float64[])
    for proc in d.processes, r in regions_of(proc, d)
        lvl = val(z[(r, proc.name)])
        lvl > tol && push!(activity, (r, proc.name, lvl))
    end

    return Results(data = d, model = m, production = production,
                   consumption = consumption, trade = trade, net_trade = net,
                   prices = prices, activity = activity)
end

"""
$(TYPEDSIGNATURES)

Slice every table of `res` for a single `region` (trade keeps both incoming and
outgoing flows).

Returns a `NamedTuple` `(; production, consumption, prices, net_trade, activity, trade)`.
"""
function for_region(res::Results, region::Symbol)
    f(df, col) = df[df[!, col] .== region, :]
    (production  = f(res.production, :region),
     consumption = f(res.consumption, :region),
     prices      = f(res.prices, :region),
     net_trade   = f(res.net_trade, :region),
     activity    = f(res.activity, :region),
     trade       = res.trade[(res.trade.from .== region) .| (res.trade.to .== region), :])
end
