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
    qs = Float64[s.q0 for s in d.demand]
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
    for s in d.demand
        push!(uses, (s.product, s.region))
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
function solve_market(d::MarketData; optimizer = Ipopt.Optimizer, silent::Bool = true,
                      qfloor::Real = 1e-6 * quantity_scale(d),
                      tol::Real = 1e-3 * quantity_scale(d))
    m = Model(optimizer)
    silent && set_silent(m)
    # Ipopt prints its licence banner even when silent unless "sb" is set
    silent && optimizer === Ipopt.Optimizer && set_attribute(m, "sb", "yes")

    demand_of = Dict((s.region, s.product) => s for s in d.demand)
    supply_of = Dict((s.region, s.product) => s for s in d.supply)

    # --- decision variables -------------------------------------------------
    # final demand quantities: floored because the marginal benefit a·D^(-1/η)
    # is unbounded at zero
    D = Dict{Tuple{Symbol,Symbol},VariableRef}()
    for s in d.demand
        D[(s.region, s.product)] = @variable(m, lower_bound = qfloor, start = s.q0)
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
        z[(r, proc.name)] = @variable(m, lower_bound = 0.0, start = 1.0)
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
    cons_benefit = @expression(m, sum(
        let s = demand_of[k], a = s.p0 * (s.shift * s.q0)^(1 / s.elasticity),
            e = 1 - 1 / s.elasticity
            abs(e) < 1e-8 ? a * log(D[k]) : (a / e) * D[k]^e
        end for k in keys(D)))
    supply_cost = @expression(m, sum(
        let s = supply_of[k], b = s.p0 * (s.shift * s.q0)^(-1 / s.elasticity),
            e = 1 + 1 / s.elasticity
            (b / e) * S[k]^e
        end for k in keys(S)))
    conv_cost = @expression(m, sum(proc.vacost * z[(r, proc.name)]
                                   for proc in d.processes for r in regions_of(proc, d)))
    trans_cost = @expression(m,
        sum(d.transport[key] * T[key] for key in keys(T); init = 0.0) +
        sum(d.transport[key] * X[key] for key in keys(X) if key[2] != key[3]; init = 0.0))

    @objective(m, Max, cons_benefit - supply_cost - conv_cost - trans_cost)

    optimize!(m)
    st = termination_status(m)
    (st == MOI.LOCALLY_SOLVED || st == MOI.OPTIMAL) ||
        @warn "solver returned status $st — results may be unreliable"

    return build_results(d, m, D, S, z, T, X, balance, origin_balance, tol)
end
