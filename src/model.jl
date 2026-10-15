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
    uses = Set{Tuple{String,String}}()
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
function floored_origins(node, acc = Set{String}())
    node isa AbstractString && return acc
    for c in node.children
        if c isa AbstractString
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
each of the others or by an `"all"` fallback.
"""
function armington_of(d::MarketData)
    specs    = Dict{Tuple{String,String},Armington}()   # (product, destination)
    fallback = Dict{String,Armington}()                 # product, from `"all"`
    products = String[]
    for a in d.armington
        a.product in d.products ||
            throw(ArgumentError("Armington specification for unknown product $(a.product)"))
        isfinite(a.sigma) && push!(products, a.product)
        if a.destination == "all"
            haskey(fallback, a.product) && throw(ArgumentError(
                "duplicate Armington specification for $(a.product)"))
            fallback[a.product] = a
        else
            for r in (a.destination isa AbstractString ? [a.destination] : a.destination)
                r in d.regions || throw(ArgumentError(
                    "Armington specification of $(a.product) for unknown region $r"))
                haskey(specs, (a.product, r)) && throw(ArgumentError(
                    "duplicate Armington specification for $(a.product) in region $r"))
                specs[(a.product, r)] = a
            end
        end
    end

    arm = Dict{Tuple{String,String},Armington}()
    for p in unique(products), r in d.regions
        a = get(specs, (p, r), get(fallback, p, nothing))
        a === nothing && throw(ArgumentError(
            "$p is an Armington product but region $r has no specification for " *
            "it: give one for every region, or one with destination = \"all\""))
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
    can = Set{Tuple{String,String}}()
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
function armington_origins(d::MarketData, a::Armington, r::String,
                           can::Set{Tuple{String,String}} = producible(d))
    p = a.product
    available = [o for o in [r; [o for o in d.regions
                                 if o != r && haskey(d.transport_costs, (p, o, r))]]
                 if (p, o) in can]
    isempty(available) && throw(ArgumentError(
        "no origin can supply $p to region $r: it is produced nowhere that $r " *
        "can buy from"))
    if isempty(a.shares)
        return available, fill(1 / length(available), length(available))
    end
    origins = String[]
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
# are origins (`String`) or deeper nodes, `shares` their value shares among
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
function armington_tree(d::MarketData, a::Armington, r::String,
                        can::Set{Tuple{String,String}} = producible(d))
    origins, shares = armington_origins(d, a, r, can)
    weight = Dict(zip(origins, shares))
    # weight of everything below a child, origins being the leaves
    subtotal(c) = c isa AbstractString ? weight[c] : sum(subtotal, c.children)

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

    placed = collect_origins!(String[], a.nests)
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
variety of an origin when `node` is a `String`, otherwise a fresh variable tied
to its children by the CES aggregator of the group.  Called on the root of an
[`armington_tree`](@ref), it adds one variable and one constraint per group.
"""
function armington_quantity(m, node, p::String, r::String, X, floor_at, in_power::Bool = false)
    node isa AbstractString && return X[(p, node, r)]
    qs = [armington_quantity(m, c, p, r, X, floor_at, isfinite(node.sigma))
          for c in node.children]
    # the composite only needs a positive floor where it is itself raised to a
    # power, i.e. where the group it belongs to substitutes imperfectly
    composite = @variable(m, lower_bound = in_power ? floor_at(sum(start_value, qs)) : 0.0,
                             start = sum(start_value, qs), base_name = "origin_composite[$p,$r]")
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
`m`: the quantity consumed of a product when `node` is a `String`, otherwise a
fresh variable tied to its members by the CES aggregator of the group. Called
on the `product` field of a [`DemandSpec`](@ref), it creates one consumption
variable per product of the group and one composite per nest.
"""
function demand_quantity(m, node, r::String, D, floor_at, start, in_power::Bool = true)
    if node isa AbstractString
        haskey(D, (r, node)) && throw(ArgumentError(
            "$(node) has more than one final demand in region $r"))
        f = floor_at(start)
        return D[(r, node)] = @variable(m, lower_bound = in_power ? f : 0.0,
                                           start = max(start, f), base_name = "D[$r,$node]")
    end
    weights = Float64[node.shares; [sub.share for sub in node.nests]]
    weights ./= sum(weights)
    members = Any[node.products; node.nests]
    qs = [demand_quantity(m, member, r, D, floor_at, weights[i] * start,
                          isfinite(node.sigma))
          for (i, member) in enumerate(members)]
    length(qs) == 1 && return qs[1]
    composite = @variable(m, lower_bound = in_power ? floor_at(sum(start_value, qs)) : 0.0,
                             start = sum(start_value, qs), base_name = "demand_composite[$r]")
    if isinf(node.sigma)
        # perfect substitutes: one good to the consumer, so the composite is
        # the plain sum, which is the limit of the CES below
        @constraint(m, sum(qs) >= composite)
    else
        ρ = (node.sigma - 1) / node.sigma
        aggr = @expression(m, sum(weights[i]^(1 / node.sigma) * qs[i]^ρ
                                  for i in eachindex(qs)))
        @constraint(m, aggr^(1 / ρ) >= composite)
    end
    return composite
end

"""
$(TYPEDSIGNATURES)

The *relative* rate of tax on a final-demand curve.  An ad valorem rate scales
the price of every member of a group alike, and so scales the group's price
index; a rate that differed between members would not, and there is no single
curve left to put it on.  Per-unit rates have no such problem — they are
charged on each product's own quantity — so only this one has to be shared.
"""
function demand_relative_tax(d::MarketData, s::DemandSpec)
    products = nest_products(s.product)
    rates = [tax_rates(d.taxes, "demand"; product = p, region = s.region).relative
             for p in products]
    allequal(rates) || throw(ArgumentError(
        "the products $(join(products, ", ")) share one demand curve in region " *
        "$(s.region) but carry different ad valorem rates, which cannot be " *
        "expressed on their composite; tax them at one rate, use per-unit rates, " *
        "or give them separate curves"))
    return first(rates)
end

"""
$(TYPEDSIGNATURES)

Build the JuMP model for `d` and solve it once, with the ad valorem trade
tariffs held at the money equivalents given in `tariff`.

Internal: [`solve_market`](@ref) is the entry point, and is what validates the
economy, acts on the solver's verdict and iterates this function when an ad
valorem rate on trade has to be found by fixed point.
"""
function solve_once(d::MarketData, tariff::Dict{Tuple{String,String,String},Float64};
                    solver = Ipopt.Optimizer, solver_options = (;),
                    silent::Bool = true, qfloor::Union{Nothing,Real} = nothing,
                    tol::Real = 1e-3 * quantity_scale(d))
    # The floor exists so that a variable raised to a power below 1 is never
    # differentiated at zero. What counts as "near zero" depends on the
    # variable: a millionth of what the data says that variable should be is
    # negligible for it, whatever the other quantities in the economy are. A
    # single absolute floor would instead be set by the smallest product in the
    # model and applied to the largest. Passing `qfloor` explicitly overrides
    # this with one absolute floor for everything.
    floor_at(reference) = qfloor === nothing ? max(1e-12, 1e-6 * reference) :
                                               float(qfloor)
    smallest_floor = floor_at(quantity_scale(d))
    m = Model(solver)
    silent && set_silent(m)
    # Ipopt, plain or wrapped in `optimizer_with_attributes`. The settings
    # below are defaults: whatever the caller passed in the wrapper is set
    # already, and is left alone.
    wrapped = solver isa MOI.OptimizerWithAttributes
    if (wrapped ? solver.optimizer_constructor : solver) === Ipopt.Optimizer
        chosen = wrapped ? Set(first(p) for p in solver.params) : Set()
        default(name, value) = MOI.RawOptimizerAttribute(name) in chosen ||
                               set_attribute(m, name, value)
        # Ipopt prints its licence banner even when silent unless "sb" is set
        silent && default("sb", "yes")
        # Its feasibility tolerance is absolute, and its default (1e-4) can be
        # coarser than the quantities of an economy written in small units: an
        # infeasible model then comes back reported as solved. Tie it to the
        # scale of the data, as `qfloor` and `tol` already are.
        # It has to be tighter than `qfloor`, or a quantity that the model
        # holds at its floor can violate a balance by more than the solver
        # minds and an infeasible economy comes back reported as solved. Two
        # orders below the floor is enough to see that, without asking for a
        # precision the problem cannot deliver.
        default("constr_viol_tol", max(1e-12, 1e-2 * smallest_floor))
        default("acceptable_constr_viol_tol", max(1e-10, smallest_floor))
    end
    # set last, so that the caller's choice overrides any of the above
    for (option, value) in (solver_options isa NamedTuple ? pairs(solver_options) :
                            solver_options)
        set_attribute(m, option isa Symbol ? String(option) : option, value)
    end

    supply_of = Dict((s.region, s.product) => s for s in d.supply)
    # the reference quantity of whatever curve a consumed product belongs to,
    # used only to start the solver off in a sensible place
    demand_of = Dict{Tuple{String,String},DemandSpec}()
    for spec in d.demand, p in nest_products(spec.product)
        demand_of[(spec.region, p)] = spec
    end

    # --- decision variables -------------------------------------------------
    # final demand quantities: floored because the marginal benefit a·D^(-1/η)
    # is unbounded at zero
    # Final demand. A plain product gets one variable; a group of substitutes
    # gets one per product plus the composite the curve actually describes, so
    # the balances and the result tables stay per product throughout.
    D = Dict{Tuple{String,String},VariableRef}()
    composites = Tuple{DemandSpec,Any}[]
    for s in d.demand
        push!(composites,
              (s, demand_quantity(m, s.product, s.region, D, floor_at, s.shift * s.q0)))
    end
    # primary supply quantities: the marginal cost b·S^(1/ε) is finite at zero,
    # so no floor is needed
    S = Dict{Tuple{String,String},VariableRef}()
    for s in d.supply
        v = @variable(m, lower_bound = 0.0, start = min(s.shift * s.q0, s.capacity), base_name = "S[$(s.region),$(s.product)]")
        isfinite(s.capacity) && set_upper_bound(v, s.capacity)
        S[(s.region, s.product)] = v
    end
    # process activity levels
    z = Dict{Tuple{String,String},VariableRef}()
    for proc in d.processes, r in regions_of(proc, d)
        cap = capacity_of(proc, r)
        v = @variable(m, lower_bound = 0.0, start = min(1.0, cap), base_name = "z[$r,$(proc.name)]")
        isfinite(cap) && set_upper_bound(v, cap)
        z[(r, proc.name)] = v
    end
    # trade flows of the homogeneous products  T[(product, from, to)]
    arm = armington_of(d)
    arm_products = Set(p for (p, _) in keys(arm))
    T = Dict{Tuple{String,String,String},VariableRef}()
    for ((p, from, to), _) in d.transport_costs
        p in arm_products && continue    # Armington products use X instead
        T[(p, from, to)] = @variable(m, lower_bound = 0.0, start = 0.0, base_name = "T[$p,$from,$to]")
    end
    # Armington products: quantity of the variety of each origin used in each
    # destination, X[(product, origin, destination)], and the composite that
    # the destination actually uses, A[(product, destination)]
    X = Dict{Tuple{String,String,String},VariableRef}()
    A = Dict{Tuple{String,String},VariableRef}()
    tree_of = Dict{Tuple{String,String},Any}()
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
        A[(p, r)] = @variable(m, lower_bound = 0.0, start = start, base_name = "A[$p,$r]")
        for (i, o) in pairs(origins)
            X[(p, o, r)] = @variable(m, lower_bound = o in floored ?
                                                      floor_at(shares[i] * start) : 0.0,
                                        start = shares[i] * start, base_name = "X[$p,$o,$r]")
        end
    end

    # --- ceilings on what may be carried between regions --------------------
    # one constraint per limit, over the weighted sum of the flows it covers,
    # whether those are homogeneous trade variables or Armington varieties
    flow_vars = merge(Dict{Tuple{String,String,String},VariableRef}(T),
                      Dict(k => v for (k, v) in X if k[2] != k[3]))
    trade_limit = Dict{String,ConstraintRef}()
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
    produced = Dict{Tuple{String,String},Any}()   # process outputs + primary supply
    used     = Dict{Tuple{String,String},Any}()   # process inputs + final demand
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
                    q = @variable(m, lower_bound = floor_at(nest.composite),
                                     start = nest.composite,
                                     base_name = "x[$r,$(proc.name),$p]")
                    push!(qs, q)
                    addto!(used, (r, p), q)
                end
                # Calibrated CES: shares are *value* shares (Σδ = 1) and enter
                # the quantity aggregator as δ^(1/σ).  This makes the implied
                # composite-input price equal the input price when prices are
                # uniform (and total physical input = composite·z), so the dual
                # prices reproduce the genuine marginal cost of production.
                if isinf(nest.sigma)
                    @constraint(m, nest.phi * sum(qs) >= nest.composite * zk)
                else
                    ρ = (nest.sigma - 1) / nest.sigma
                    aggr = @expression(m, sum(nest.shares[i]^(1 / nest.sigma) * qs[i]^ρ
                                              for i in eachindex(qs)))
                    @constraint(m, nest.phi * aggr^(1 / ρ) >= nest.composite * zk)
                end
            end
        end
    end

    # --- balances per (region, product)  → duals are the prices -------------
    # `balance` is what a region's users face (its dual is the price they pay),
    # `origin_balance` is what its producers sell (its dual is the price they
    # get).  For a homogeneous product the two coincide and only `balance` is
    # built; for an Armington product they differ by the composition of the
    # CES bundle, and both are built.
    balance        = Dict{Tuple{String,String},ConstraintRef}()
    origin_balance = Dict{Tuple{String,String},ConstraintRef}()
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
                @constraint(m, armington_quantity(m, tree_of[(p, r)], p, r, X, floor_at) >=
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
        let w = 1 + demand_relative_tax(d, s),
            a = s.p0 * (s.shift * s.q0)^(1 / s.elasticity) / w,
            e = 1 - 1 / s.elasticity
            (abs(e) < 1e-8 ? a * log(C) : (a / e) * C^e) -
            # per-unit rates are charged on each product's own quantity: the
            # composite is an index, and a tax per m³ is not a tax per index
            # point unless every member happens to cost the same
            sum(tax_rates(d.taxes, "demand"; product = p, region = s.region).absolute /
                w * D[(s.region, p)] for p in nest_products(s.product))
        end for (s, C) in composites))
    # A tax on production sits between the market price and what the producer
    # keeps, so it raises the marginal cost the market has to cover.
    supply_cost = @expression(m, sum(
        let s = supply_of[k],
            t = tax_rates(d.taxes, "supply"; product = s.product, region = s.region),
            b = (1 + t.relative) * s.p0 * (s.shift * s.q0)^(-1 / s.elasticity),
            e = 1 + 1 / s.elasticity
            (b / e) * S[k]^e + t.absolute * S[k]
        end for k in keys(S)))
    conv_cost = @expression(m, sum(
        let t = tax_rates(d.taxes, "process"; process = proc.name, region = r)
            ((1 + t.relative) * proc.vacost + t.absolute) * z[(r, proc.name)]
        end for proc in d.processes for r in regions_of(proc, d)))
    route_cost(key) = d.transport_costs[key] + get(tariff, key, 0.0) +
        tax_rates(d.taxes, "trade"; product = key[1], from = key[2], to = key[3]).absolute
    trans_cost = @expression(m,
        sum(route_cost(key) * T[key] for key in keys(T); init = 0.0) +
        sum(route_cost(key) * X[key] for key in keys(X) if key[2] != key[3]; init = 0.0))

    @objective(m, Max, cons_benefit - supply_cost - conv_cost - trans_cost)

    wall_time = @elapsed optimize!(m)

    return build_results(d, m, D, S, z, T, X, balance, origin_balance, trade_limit,
                         tariff, tol, wall_time)
end

"""
$(TYPEDSIGNATURES)

Build and solve the equilibrium for the economy `d` and return a
[`Results`](@ref).

If the solver does not report an optimal (or locally optimal) solution,
`on_failure` decides what happens: `"error"` (the default) throws, `"warn"`
warns and `"ignore"` says nothing. The status, the solver's own message, the
time spent and the iterations taken are on the [`Results`](@ref).

* `solver` is any JuMP solver able to handle nonlinear constraints and to
  return duals (Ipopt by default), and `silent` suppresses its output.
* `solver_options` are passed on to it, as a `NamedTuple`, a `Dict` or a vector
  of pairs: `(max_iter = 500, tol = 1e-9)`, or
  `Dict("max_iter" => 500, MOI.TimeLimitSec() => 60.0)`. A name is the
  solver's own name for the option, and an MOI attribute is passed as it is.
  They are set after Markets' own settings, and so override them.
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

That iteration is successive substitution, which converges only where the map
contracts — quickly in every case tried here, including rates of several
hundred per cent, but without a guarantee. The step is therefore halved
whenever an iteration increases the error rather than reducing it, which
rescues an oscillating case and leaves a well-behaved one at full speed.
"""
function solve_market(d::MarketData; max_tax_iterations::Int = 50,
                      tax_tolerance::Real = 1e-8, on_failure::Name = "error",
                      kwargs...)
    on_failure = nm(on_failure)
    on_failure in ("error", "warn", "ignore") || throw(ArgumentError(
        "on_failure is \"error\", \"warn\" or \"ignore\", got $on_failure"))
    validate(d)
    advalorem = [t for t in d.taxes if t.on == "trade" && t.relative]
    tariff = Dict{Tuple{String,String,String},Float64}()
    res = solve_once(d, tariff; kwargs...)
    isempty(advalorem) && return report(res, on_failure)

    routes = collect(keys(d.transport_costs))
    previous_gap = Inf
    step = 1.0
    for _ in 1:max_tax_iterations
        # the money a rate of `relative` comes to, at the prices just found
        origin_price(p, o) = begin
            rows = res.prices[(res.prices.region .== o) .&
                              (res.prices.product .== p), :]
            isempty(rows) ? 0.0 : coalesce(only(rows.producer_price), 0.0)
        end
        updated = Dict(key => tax_rates(d.taxes, "trade"; product = key[1], from = key[2],
                                        to = key[3]).relative * origin_price(key[1], key[2])
                       for key in routes)
        gap = maximum(abs(updated[k] - get(tariff, k, 0.0)) / max(abs(updated[k]), 1.0)
                      for k in routes; init = 0.0)
        # Successive substitution has no convergence guarantee: if the error
        # grows the map is not contracting at this step length, so halve it.
        # A full step is kept while the error falls, which is the usual case
        # and costs nothing.
        gap > previous_gap && (step = max(step / 2, 0.05))
        previous_gap = gap
        tariff = Dict(k => step * updated[k] + (1 - step) * get(tariff, k, 0.0)
                      for k in routes)
        res = tally(solve_once(d, tariff; kwargs...), res)
        gap <= tax_tolerance && return report(res, on_failure)
    end
    @warn "ad valorem trade taxes did not settle in $max_tax_iterations solves"
    return report(res, on_failure)
end

"""
$(TYPEDSIGNATURES)

`res`, with the solver's time and iterations summed over it and the solves
in `earlier`: when the ad valorem tariffs are found by iteration, the caller
waited for all of them, not just the last.
"""
function tally(res::Results, earlier::Results)
    fields = (; (f => getfield(res, f) for f in fieldnames(Results))...)
    return Results(; merge(fields, (solve_time = earlier.solve_time + res.solve_time,
                                    iterations = earlier.iterations + res.iterations,
                                    solves = earlier.solves + res.solves))...)
end

"""
$(TYPEDSIGNATURES)

Act on the solver's verdict: a solution that is neither optimal nor locally
optimal is an error by default, since the tables that come with it look exactly
like those of a solved model.  `on_failure` chooses between throwing, warning
and saying nothing.
"""
function report(res::Results, on_failure::AbstractString)
    res.status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) && return res
    if res.status in (MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED)
        # the solver stopped at its looser, "acceptable" tolerance: the
        # solution is an equilibrium, just a less exact one
        on_failure == "ignore" || @warn "the solver returned $(res.status): the " *
            "solution satisfies its acceptable tolerance but not its target one, " *
            "so treat the last digits with suspicion"
        return res
    end
    message = "the solver returned $(res.status): the economy may be infeasible, " *
              "unbounded or too hard to solve, and the tables returned with this " *
              "result do not describe an equilibrium"
    on_failure == "error" && throw(ErrorException(message))
    on_failure == "warn"  && @warn message
    return res
end
