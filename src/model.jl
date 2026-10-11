# ----------------------------------------------------------------------------
#  Model construction
# ----------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

The scale of the economy: the smallest reference quantity in its demand and
supply curves, or 1 if it has none.  The quantity floor and the reporting
threshold are taken relative to it, so that they mean the same thing whether
the model is written in m³ or in million m³.
"""
function quantity_scale(d::MarketData)
    # a curve over a group of substitutes describes their composite, so its
    # reference quantity is divided among them for the purpose of scale
    qs = Float64[s.q0 / length(nest_products(s.product)) for s in d.demand]
    append!(qs, Float64[s.q0 for s in d.supply])
    return isempty(qs) ? 1.0 : minimum(qs)
end

"""
$(TYPEDSIGNATURES)

The `(product, region)` pairs something in the region uses: final demand, or an
input of a process operating there.  A region that uses none of a product needs
no composite of it, though it may still produce it for others.
"""
function consumed(d::MarketData)
    uses = Set{Tuple{Symbol,Symbol}}()
    for s in d.demand, p in nest_products(s.product)
        push!(uses, (p, s.region))
    end
    for proc in d.processes, r in regions_of(proc, d), nest in proc.inputs, p in nest.products
        push!(uses, (p, r))
    end
    return uses
end

"""
$(TYPEDSIGNATURES)

The origins that sit directly under a CES node of `node` and therefore enter an
aggregator as ``X^\\rho``.  Only those need a positive lower bound; a variety
pooled by an infinitely substitutable group enters linearly and may go to zero.
"""
function floored_origins(node, acc = Set{Symbol}())
    node isa Symbol && return acc
    for c in node.children
        if c isa Symbol
            isfinite(node.sigma) && push!(acc, c)
        else
            floored_origins(c, acc)
        end
    end
    return acc
end

"""
$(TYPEDSIGNATURES)

Resolve the [`Armington`](@ref) specifications of `d` into a
`(product, destination) => `[`Armington`](@ref) dictionary, covering every
destination of every product whose varieties are imperfect substitutes.
Products absent from it are homogeneous.

A product is Armington in every destination or in none, so a specification with
a finite elasticity for some destinations must be matched either by one for
each of the others or by an `:all` fallback.
"""
function armington_of(d::MarketData)
    specs    = Dict{Tuple{Symbol,Symbol},Armington}()   # (product, destination)
    fallback = Dict{Symbol,Armington}()                 # product, from `:all`
    products = Symbol[]
    for a in d.armington
        a.product in d.products ||
            throw(ArgumentError("Armington specification for unknown product $(a.product)"))
        isfinite(a.sigma) && push!(products, a.product)
        if a.destination === :all
            haskey(fallback, a.product) && throw(ArgumentError(
                "duplicate Armington specification for $(a.product)"))
            fallback[a.product] = a
        else
            for r in (a.destination isa Symbol ? [a.destination] : a.destination)
                r in d.regions || throw(ArgumentError(
                    "Armington specification of $(a.product) for unknown region $r"))
                haskey(specs, (a.product, r)) && throw(ArgumentError(
                    "duplicate Armington specification for $(a.product) in region $r"))
                specs[(a.product, r)] = a
            end
        end
    end

    arm = Dict{Tuple{Symbol,Symbol},Armington}()
    for p in unique(products), r in d.regions
        a = get(specs, (p, r), get(fallback, p, nothing))
        a === nothing && throw(ArgumentError(
            "$p is an Armington product but region $r has no specification for " *
            "it: give one for every region, or one with destination = :all"))
        isfinite(a.sigma) || throw(ArgumentError(
            "$p is an Armington product in some regions but homogeneous in $r: " *
            "a product is either homogeneous everywhere or imperfectly " *
            "substitutable everywhere"))
        arm[(p, r)] = a
    end
    return arm
end

"""
$(TYPEDSIGNATURES)

The `(product, region)` pairs the economy can produce at all: those with a
primary supply curve, and those a process operating in the region yields.  A
region that cannot produce a product has no variety of it, so it is nobody's
origin for it — not even its own.
"""
function producible(d::MarketData)
    can = Set{Tuple{Symbol,Symbol}}()
    for s in d.supply
        push!(can, (s.product, s.region))
    end
    for proc in d.processes, r in regions_of(proc, d), (p, _) in proc.outputs
        push!(can, (p, r))
    end
    return can
end

"""
$(TYPEDSIGNATURES)

The origins that destination `r` can buy product `a.product` from — itself,
plus every origin with a transport route into `r`, keeping only those that can
produce it — together with their normalised value shares.  Origins without a
positive share are left out.
"""
function armington_origins(d::MarketData, a::Armington, r::Symbol,
                           can::Set{Tuple{Symbol,Symbol}} = producible(d))
    p = a.product
    available = [o for o in [r; [o for o in d.regions
                                 if o != r && haskey(d.transport, (p, o, r))]]
                 if (p, o) in can]
    isempty(available) && throw(ArgumentError(
        "no origin can supply $p to region $r: it is produced nowhere that $r " *
        "can buy from"))
    if isempty(a.shares)
        return available, fill(1 / length(available), length(available))
    end
    origins = Symbol[]
    shares  = Float64[]
    for o in available
        s = get(a.shares, (o, r), 0.0)
        s > 0 || continue
        push!(origins, o)
        push!(shares, s)
    end
    isempty(origins) && throw(ArgumentError(
        "no origin with a positive Armington share supplies $p in region $r"))
    shares ./= sum(shares)
    return origins, shares
end

# A node of the CES tree a destination aggregates its origins with: `children`
# are origins (`Symbol`) or deeper nodes, `shares` their value shares among
# themselves, `sigma` the elasticity they substitute at.
struct ArmNode
    sigma::Float64
    children::Vector{Any}
    shares::Vector{Float64}
end

"""
$(TYPEDSIGNATURES)

Lay out the origins available to destination `r` as the CES tree described by
the groups of `a`: origins placed in a group are aggregated by it first, and
every other origin is a direct member of the root.  Each node's shares are
those of the origins beneath it, normalised among its siblings.  Groups that no
available origin belongs to are dropped, and a group left with a single member
collapses into it.
"""
function armington_tree(d::MarketData, a::Armington, r::Symbol,
                        can::Set{Tuple{Symbol,Symbol}} = producible(d))
    origins, shares = armington_origins(d, a, r, can)
    weight = Dict(zip(origins, shares))
    # weight of everything below a child, origins being the leaves
    subtotal(c) = c isa Symbol ? weight[c] : sum(subtotal, c.children)

    # a group becomes a node over the available origins it holds, or `nothing`
    function group(n::OriginNest)
        members = Any[o for o in n.origins if haskey(weight, o)]
        for sub in n.nests
            s = group(sub)
            s === nothing || push!(members, s)
        end
        isempty(members) && return nothing
        length(members) == 1 && return members[1]    # a lone member is no group
        return ArmNode(n.sigma, members, Float64[])  # shares filled in below
    end

    placed = collect_origins!(Symbol[], a.nests)
    children = Any[]
    for n in a.nests
        s = group(n)
        s === nothing || push!(children, s)
    end
    for o in origins                     # origins in no group sit at the root
        o in placed || push!(children, o)
    end

    # give every node the shares of its children among themselves
    function normalise(c)
        c isa ArmNode || return c
        total = sum(subtotal, c.children)
        return ArmNode(c.sigma, Any[normalise(k) for k in c.children],
                       Float64[subtotal(k) / total for k in c.children])
    end
    root = length(children) == 1 ? children[1] :
           normalise(ArmNode(a.sigma, children, Float64[]))
    return root, origins, shares
end

"""
$(TYPEDSIGNATURES)

The quantity of the composite that `node` stands for, as a variable of `m`: the
variety of an origin when `node` is a `Symbol`, otherwise a fresh variable tied
to its children by the CES aggregator of the group.  Called on the root of an
[`armington_tree`](@ref), it adds one variable and one constraint per group.
"""
function armington_quantity(m, node, p::Symbol, r::Symbol, X, qfloor, in_power::Bool = false)
    node isa Symbol && return X[(p, node, r)]
    qs = [armington_quantity(m, c, p, r, X, qfloor, isfinite(node.sigma))
          for c in node.children]
    # the composite only needs a positive floor where it is itself raised to a
    # power, i.e. where the group it belongs to substitutes imperfectly
    composite = @variable(m, lower_bound = in_power ? qfloor : 0.0,
                             start = sum(start_value, qs))
    if isinf(node.sigma)
        # perfect substitutes: the group is one good, pooled by plain addition
        # (the limit of the CES below, where δ^(1/σ) → 1 and ρ → 1)
        @constraint(m, sum(qs) >= composite)
    else
        ρ = (node.sigma - 1) / node.sigma
        aggr = @expression(m, sum(node.shares[i]^(1 / node.sigma) * qs[i]^ρ
                                  for i in eachindex(qs)))
        @constraint(m, aggr^(1 / ρ) >= composite)
    end
    return composite
end

"""
$(TYPEDSIGNATURES)

The quantity of the demand composite that `node` stands for, as a variable of
`m`: the quantity consumed of a product when `node` is a `Symbol`, otherwise a
fresh variable tied to its members by the CES aggregator of the group. Called
on the `product` field of a [`DemandSpec`](@ref), it creates one consumption
variable per product of the group and one composite per nest.
"""
function demand_quantity(m, node, r::Symbol, D, qfloor, start, in_power::Bool = true)
    if node isa Symbol
        haskey(D, (r, node)) && throw(ArgumentError(
            "$(node) has more than one final demand in region $r"))
        return D[(r, node)] = @variable(m, lower_bound = in_power ? qfloor : 0.0,
                                           start = max(start, qfloor))
    end
    weights = Float64[node.shares; [sub.share for sub in node.nests]]
    weights ./= sum(weights)
    members = Any[node.products; node.nests]
    qs = [demand_quantity(m, member, r, D, qfloor, weights[i] * start,
                          isfinite(node.sigma))
          for (i, member) in enumerate(members)]
    length(qs) == 1 && return qs[1]
    composite = @variable(m, lower_bound = in_power ? qfloor : 0.0,
                             start = sum(start_value, qs))
    ρ = (node.sigma - 1) / node.sigma
    aggr = @expression(m, sum(weights[i]^(1 / node.sigma) * qs[i]^ρ for i in eachindex(qs)))
    @constraint(m, aggr^(1 / ρ) >= composite)
    return composite
end

"""
$(TYPEDSIGNATURES)

The tax falling on a final-demand curve.  A curve over a group of substitutes
describes their composite, so a tax on it must be one the whole group shares:
if the products of the group are taxed differently, that is an error rather
than something to average away.
"""
function demand_tax(d::MarketData, s::DemandSpec)
    products = nest_products(s.product)
    rates = [tax_rates(d.taxes, :demand; product = p, region = s.region) for p in products]
    allequal(rates) || throw(ArgumentError(
        "the products $(join(products, ", ")) share one demand curve in region " *
        "$(s.region) but are taxed differently; tax them alike, or give them " *
        "separate curves"))
    return first(rates)
end

"""
$(TYPEDSIGNATURES)

Build and solve the equilibrium for the economy `d` and return a
[`Results`](@ref).

A warning is emitted if the solver does not report an optimal (or locally
optimal) solution.

* `optimizer` is any JuMP solver able to handle nonlinear constraints and to
  return duals, and `silent` suppresses its output.
* `qfloor` is the lower bound given to the quantities that are raised to a
  power below 1, which an interior-point solver cannot differentiate at zero.
* `tol` is the quantity below which a flow, an output or an activity is treated
  as zero and left out of the result tables.

Both default to a multiple of the smallest reference quantity in `d`, so that
they scale with the units the economy is written in.
"""
function solve_once(d::MarketData, tariff::Dict{Tuple{Symbol,Symbol,Symbol},Float64};
                    optimizer = Ipopt.Optimizer, silent::Bool = true,
                    qfloor::Real = 1e-6 * quantity_scale(d),
                    tol::Real = 1e-3 * quantity_scale(d))
    m = Model(optimizer)
    silent && set_silent(m)
    # Ipopt prints its licence banner even when silent unless "sb" is set
    silent && optimizer === Ipopt.Optimizer && set_attribute(m, "sb", "yes")

    supply_of = Dict((s.region, s.product) => s for s in d.supply)
    # the reference quantity of whatever curve a consumed product belongs to,
    # used only to start the solver off in a sensible place
    demand_of = Dict{Tuple{Symbol,Symbol},DemandSpec}()
    for spec in d.demand, p in nest_products(spec.product)
        demand_of[(spec.region, p)] = spec
    end

    # --- decision variables -------------------------------------------------
    # final demand quantities: floored because the marginal benefit a·D^(-1/η)
    # is unbounded at zero
    # Final demand. A plain product gets one variable; a group of substitutes
    # gets one per product plus the composite the curve actually describes, so
    # the balances and the result tables stay per product throughout.
    D = Dict{Tuple{Symbol,Symbol},VariableRef}()
    composites = Tuple{DemandSpec,Any}[]
    for s in d.demand
        push!(composites,
              (s, demand_quantity(m, s.product, s.region, D, qfloor, s.shift * s.q0)))
    end
    # primary supply quantities: the marginal cost b·S^(1/ε) is finite at zero,
    # so no floor is needed
    S = Dict{Tuple{Symbol,Symbol},VariableRef}()
    for s in d.supply
        v = @variable(m, lower_bound = 0.0, start = min(s.shift * s.q0, s.capacity))
        isfinite(s.capacity) && set_upper_bound(v, s.capacity)
        S[(s.region, s.product)] = v
    end
    # process activity levels
    z = Dict{Tuple{Symbol,Symbol},VariableRef}()
    for proc in d.processes, r in regions_of(proc, d)
        cap = capacity_of(proc, r)
        v = @variable(m, lower_bound = 0.0, start = min(1.0, cap))
        isfinite(cap) && set_upper_bound(v, cap)
        z[(r, proc.name)] = v
    end
    # trade flows of the homogeneous products  T[(product, from, to)]
    arm = armington_of(d)
    arm_products = Set(p for (p, _) in keys(arm))
    T = Dict{Tuple{Symbol,Symbol,Symbol},VariableRef}()
    for ((p, from, to), _) in d.transport
        p in arm_products && continue    # Armington products use X instead
        T[(p, from, to)] = @variable(m, lower_bound = 0.0, start = 0.0)
    end
    # Armington products: quantity of the variety of each origin used in each
    # destination, X[(product, origin, destination)], and the composite that
    # the destination actually uses, A[(product, destination)]
    X = Dict{Tuple{Symbol,Symbol,Symbol},VariableRef}()
    A = Dict{Tuple{Symbol,Symbol},VariableRef}()
    tree_of = Dict{Tuple{Symbol,Symbol},Any}()
    can  = producible(d)
    uses = consumed(d)
    for ((p, r), a) in arm
        # a region that uses none of the product needs no composite of it; it
        # may still produce it for others, which the origin balance covers
        (p, r) in uses || continue
        root, origins, shares = armington_tree(d, a, r, can)
        tree_of[(p, r)] = root
        floored = floored_origins(root)
        # a rough but useful starting point: the local reference consumption
        a0 = get(demand_of, (r, p), nothing)
        start = a0 === nothing ? 1.0 : a0.q0
        A[(p, r)] = @variable(m, lower_bound = 0.0, start = start)
        for (i, o) in pairs(origins)
            X[(p, o, r)] = @variable(m, lower_bound = o in floored ? qfloor : 0.0,
                                        start = shares[i] * start)
        end
    end

    # --- ceilings on what may be carried between regions --------------------
    # one constraint per limit, over the weighted sum of the flows it covers,
    # whether those are homogeneous trade variables or Armington varieties
    flow_vars = merge(Dict{Tuple{Symbol,Symbol,Symbol},VariableRef}(T),
                      Dict(k => v for (k, v) in X if k[2] != k[3]))
    trade_limit = Dict{Symbol,ConstraintRef}()
    for limit in d.trade_limits
        haskey(trade_limit, limit.name) && throw(ArgumentError(
            "duplicate trade limit named $(limit.name)"))
        covered = [k for k in keys(flow_vars) if covers(limit, k...)]
        isempty(covered) && @warn "trade limit $(limit.name) covers no flow of this economy"
        trade_limit[limit.name] = @constraint(m,
            sum(get(limit.weights, k[1], 1.0) * flow_vars[k] for k in covered;
                init = 0.0) <= limit.capacity)
    end

    # --- assemble per-(region,product) sources & uses -----------------------
    produced = Dict{Tuple{Symbol,Symbol},Any}()   # process outputs + primary supply
    used     = Dict{Tuple{Symbol,Symbol},Any}()   # process inputs + final demand
    addto!(dict, key, term) = (dict[key] = haskey(dict, key) ? dict[key] + term : term)

    for proc in d.processes, r in regions_of(proc, d)
        zk = z[(r, proc.name)]
        for (p, yield) in proc.outputs
            addto!(produced, (r, p), yield * zk)
        end
        for nest in proc.inputs
            if length(nest.products) == 1            # Leontief
                addto!(used, (r, nest.products[1]), nest.composite * zk)
            else                                     # CES
                qs = VariableRef[]
                for p in nest.products
                    q = @variable(m, lower_bound = qfloor, start = nest.composite)
                    push!(qs, q)
                    addto!(used, (r, p), q)
                end
                # Calibrated CES: shares are *value* shares (Σδ = 1) and enter
                # the quantity aggregator as δ^(1/σ).  This makes the implied
                # composite-input price equal the input price when prices are
                # uniform (and total physical input = composite·z), so the dual
                # prices reproduce the genuine marginal cost of production.
                ρ = (nest.sigma - 1) / nest.sigma
                aggr = @expression(m, sum(nest.shares[i]^(1 / nest.sigma) * qs[i]^ρ
                                          for i in eachindex(qs)))
                @constraint(m, nest.phi * aggr^(1 / ρ) >= nest.composite * zk)
            end
        end
    end

    # --- balances per (region, product)  → duals are the prices -------------
    # `balance` is what a region's users face (its dual is the price they pay),
    # `origin_balance` is what its producers sell (its dual is the price they
    # get).  For a homogeneous product the two coincide and only `balance` is
    # built; for an Armington product they differ by the composition of the
    # CES bundle, and both are built.
    balance        = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    origin_balance = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    for r in d.regions, p in d.products
        src = get(produced, (r, p), 0.0)
        use = get(used, (r, p), 0.0)
        haskey(S, (r, p)) && (src += S[(r, p)])
        haskey(D, (r, p)) && (use += D[(r, p)])

        if p in arm_products
            # everything region r produces goes to one of the destinations that
            # buy its variety (itself included)
            # a region that cannot produce the product has no variety of it,
            # hence nothing to ship and no producer price
            if (p, r) in can
                ships = @expression(m, sum(X[(p, r, dst)] for dst in d.regions
                                           if haskey(X, (p, r, dst)); init = 0.0))
                origin_balance[(r, p)] = @constraint(m, src - ships == 0.0)
            end
            # where the region uses the product, the composite covers those
            # uses and is the (possibly nested) CES aggregate of the varieties
            # bought: one composite variable and one constraint per group
            if haskey(A, (p, r))
                balance[(r, p)] = @constraint(m, A[(p, r)] - use == 0.0)
                @constraint(m, armington_quantity(m, tree_of[(p, r)], p, r, X, qfloor) >=
                               A[(p, r)])
            end
        else
            # homogeneous product: imports into r minus exports from r
            imp = @expression(m, sum(T[(p, o, r)] for o in d.regions if haskey(T, (p, o, r)); init = 0.0))
            exp = @expression(m, sum(T[(p, r, o)] for o in d.regions if haskey(T, (p, r, o)); init = 0.0))
            balance[(r, p)] = @constraint(m, src + imp - use - exp == 0.0)
        end
    end

    # --- objective: net social surplus --------------------------------------
    #   + gross consumer benefit  Σ ∫₀ᴰ P(d)dd
    #   − resource (supply) cost   Σ ∫₀ˢ P(s)ds
    #   − conversion (value-added) cost
    #   − transport cost
    # The benefit term is an antiderivative of the inverse demand, not the
    # integral from zero: for η ≤ 1 that integral diverges, but the constant it
    # differs by changes neither the optimum nor the prices.  It is concave and
    # increasing for every η > 0 — for η < 1 both a/e and e are negative — and
    # at η = 1 the antiderivative is a·log(D).
    # A tax on consumption sits between what the consumer pays, which is what
    # the curve describes, and the market price, which the objective's gradient
    # has to equal: hence the division by (1 + relative) and the term in
    # `absolute` below.
    cons_benefit = @expression(m, sum(
        let t = demand_tax(d, s), w = 1 + t.relative,
            a = s.p0 * (s.shift * s.q0)^(1 / s.elasticity) / w,
            e = 1 - 1 / s.elasticity
            (abs(e) < 1e-8 ? a * log(C) : (a / e) * C^e) - (t.absolute / w) * C
        end for (s, C) in composites))
    # A tax on production sits between the market price and what the producer
    # keeps, so it raises the marginal cost the market has to cover.
    supply_cost = @expression(m, sum(
        let s = supply_of[k],
            t = tax_rates(d.taxes, :supply; product = s.product, region = s.region),
            b = (1 + t.relative) * s.p0 * (s.shift * s.q0)^(-1 / s.elasticity),
            e = 1 + 1 / s.elasticity
            (b / e) * S[k]^e + t.absolute * S[k]
        end for k in keys(S)))
    conv_cost = @expression(m, sum(
        let t = tax_rates(d.taxes, :process; process = proc.name, region = r)
            ((1 + t.relative) * proc.vacost + t.absolute) * z[(r, proc.name)]
        end for proc in d.processes for r in regions_of(proc, d)))
    route_cost(key) = d.transport[key] + get(tariff, key, 0.0) +
        tax_rates(d.taxes, :trade; product = key[1], from = key[2], to = key[3]).absolute
    trans_cost = @expression(m,
        sum(route_cost(key) * T[key] for key in keys(T); init = 0.0) +
        sum(route_cost(key) * X[key] for key in keys(X) if key[2] != key[3]; init = 0.0))

    @objective(m, Max, cons_benefit - supply_cost - conv_cost - trans_cost)

    optimize!(m)
    st = termination_status(m)
    (st == MOI.LOCALLY_SOLVED || st == MOI.OPTIMAL) ||
        @warn "solver returned status $st — results may be unreliable"

    return build_results(d, m, D, S, z, T, X, balance, origin_balance, trade_limit,
                         tariff, tol)
end

"""
$(TYPEDSIGNATURES)

Build and solve the equilibrium for the economy `d` and return a
[`Results`](@ref).

A warning is emitted if the solver does not report an optimal (or locally
optimal) solution.

* `optimizer` is any JuMP solver able to handle nonlinear constraints and to
  return duals, and `silent` suppresses its output.
* `qfloor` is the lower bound given to the quantities that are raised to a
  power below 1, which an interior-point solver cannot differentiate at zero.
* `tol` is the quantity below which a flow, an output or an activity is treated
  as zero and left out of the result tables.

Both default to a multiple of the smallest reference quantity in `d`, so that
they scale with the units the economy is written in.

An **ad valorem tax on trade** is a fraction of the price at the origin, which
is a dual of the problem and so cannot be written into its objective. When `d`
has one, the equilibrium is found by repeated solves: each holds the tariff at
the money equivalent implied by the previous solution's producer prices, and
the loop stops when the two agree to `tax_tolerance` (relative), or warns after
`max_tax_iterations`. Everything else is solved exactly, once.
"""
function solve_market(d::MarketData; max_tax_iterations::Int = 50,
                      tax_tolerance::Real = 1e-8, kwargs...)
    advalorem = [t for t in d.taxes if t.on === :trade && t.relative]
    tariff = Dict{Tuple{Symbol,Symbol,Symbol},Float64}()
    res = solve_once(d, tariff; kwargs...)
    isempty(advalorem) && return res

    routes = collect(keys(d.transport))
    for _ in 1:max_tax_iterations
        # the money a rate of `relative` comes to, at the prices just found
        origin_price(p, o) = begin
            rows = res.prices[(res.prices.region .== o) .& (res.prices.product .== p), :]
            isempty(rows) ? 0.0 : coalesce(only(rows.producer_price), 0.0)
        end
        updated = Dict(key => tax_rates(d.taxes, :trade; product = key[1], from = key[2],
                                        to = key[3]).relative * origin_price(key[1], key[2])
                       for key in routes)
        gap = maximum(abs(updated[k] - get(tariff, k, 0.0)) / max(abs(updated[k]), 1.0)
                      for k in routes; init = 0.0)
        tariff = updated
        res = solve_once(d, tariff; kwargs...)
        gap <= tax_tolerance && return res
    end
    @warn "ad valorem trade taxes did not settle in $max_tax_iterations solves"
    return res
end
