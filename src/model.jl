# Part of Markets.jl package [https://github.com/sylvaticus/Markets.jl] - MIT Licence
# 20261030 verified Antonello Lobianco. Still to fully grasp ces_trees! and armington_of functions TODO 

# ----------------------------------------------------------------------------
#  Model construction
# ----------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

The scale of the economy: the smallest reference quantity in its demand and
supply curves, or 1 if it has none.  The quantity floor and the reporting
threshold are taken relative to it, so that they mean the same thing whether
the model is written in tonnes or in million tonnes.
"""
function quantity_scale(d::MarketData)
    # a curve over a group of substitutes describes their composite, so its
    # reference quantity is divided among them for the purpose of scale; a
    # flat curve has no quantity of its own, whatever its q0 says
    qs = Float64[reference_quantity(s) / length(leaves(s.product))
                 for s in d.demand if isfinite(s.elasticity)]
    append!(qs, Float64[reference_quantity(s) for s in d.supply if isfinite(s.elasticity)])
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
    for s in d.demand, p in leaves(s.product)
        push!(uses, (p, s.region))
    end
    for proc in d.processes, r in regions_of(proc, d), nest in proc.inputs, p in nest.products
        push!(uses, (p, r))
    end
    return uses
end

# Whether the members of a CES group of elasticity `σ` enter its aggregator
# raised to a power, and so need a positive floor: all but the two linear
# cases, Leontief (σ = 0) and perfect substitutes (σ = ∞).
powered(σ) = 0 < σ < Inf

"""
$(TYPEDSIGNATURES)

Require the CES aggregate of the quantities `qs`, with value shares `δ` and
elasticity of substitution `σ`, times `scale`, to be at least `rhs`. The
aggregate is

    ( Σᵢ δᵢ^(1/σ) · qᵢ^ρ )^(1/ρ),     ρ = (σ-1)/σ

and, at the three values of `σ` where that formula cannot be evaluated, its
limit there:

* `σ = 0`, Leontief: `minᵢ qᵢ/δᵢ`, written as one linear constraint per member;
* `σ = 1`, Cobb–Douglas: `Πᵢ (qᵢ/δᵢ)^δᵢ`;
* `σ = ∞`, perfect substitutes: `Σᵢ qᵢ`.

Every one of them is concave in `qs`, so the set where it is above `rhs` is
convex. The price of the composite is the CES index `(Σᵢ δᵢ pᵢ^(1-σ))^(1/(1-σ))`
of the members' prices, and its limits: `Σᵢ δᵢ pᵢ` at 0, `Πᵢ pᵢ^δᵢ` at 1.
"""
function ces_constraint!(m, qs, δ, σ, rhs, name; scale = 1.0)
    scaled(x) = isone(scale) ? x : scale * x
    if isinf(σ)
        # perfect substitutes: the group is one good, pooled by plain addition
        @constraint(m, scaled(sum(qs)) >= rhs, base_name = name)
    elseif iszero(σ)
        # fixed proportions: each member is needed in its share of the composite
        for i in eachindex(qs)
            @constraint(m, scaled(qs[i]) >= δ[i] * rhs, base_name = "$name[$i]")
        end
    elseif isone(σ)
        # the exponents must sum to 1, both for this to be the limit of the
        # CES and for it to be concave
        w = δ ./ sum(δ)
        @constraint(m, scaled(prod((qs[i] / w[i])^w[i] for i in eachindex(qs))) >= rhs,
                    base_name = name)
    else
        ρ = (σ - 1) / σ
        aggr = @expression(m, sum(δ[i]^(1 / σ) * qs[i]^ρ for i in eachindex(qs)))
        @constraint(m, scaled(aggr^(1 / ρ)) >= rhs, base_name = name)
    end
end

# --- one CES tree for groups of origins, of substitutes and of inputs ---------

"""
$(TYPEDEF)

A node of a CES tree as the model builds it: `children` are leaves (`String`s
— origins, or products) or deeper nodes, `shares` their value shares among
themselves, summing to 1, and `sigma` the elasticity they substitute at. It is
the normalised form an [`OriginNest`](@ref), a [`DemandNest`](@ref) and a
[`Nest`](@ref) of process inputs are read into, so that one builder,
[`ces_tree!`](@ref), serves all three.
"""
struct CESNode
    sigma::Float64
    children::Vector{Any}
    shares::Vector{Float64}
end

"""
$(TYPEDSIGNATURES)

`node` with every group of one member replaced by that member, at every
level: a group of one is no group, and dropping it keeps the aggregator, the
floors and the printed model free of trivial nodes.
"""
collapse(leaf::AbstractString) = leaf
function collapse(node::CESNode)
    length(node.children) == 1 && return collapse(only(node.children))
    return CESNode(node.sigma, Any[collapse(c) for c in node.children], node.shares)
end

"""
$(TYPEDSIGNATURES)

The [`CESNode`](@ref) tree of a [`DemandNest`](@ref): its products and
sub-groups, weighted by their `shares` and `share` normalised over the node,
with groups of one collapsed.
"""
function tree(n::DemandNest)
    weights = Float64[n.shares; [sub.share for sub in n.nests]]
    return collapse(CESNode(n.sigma, Any[n.products; [tree(sub) for sub in n.nests]],
                            weights ./ sum(weights)))
end

"""
$(TYPEDSIGNATURES)

The leaves that sit directly under a node of `node` that raises them to a
power. Only those need a positive lower bound; a leaf under a Leontief node or
pooled by an infinitely substitutable group enters linearly and may go to
zero.
"""
function floored_leaves(node, acc = Set{String}())
    node isa AbstractString && return acc
    for c in node.children
        if c isa AbstractString
            powered(node.sigma) && push!(acc, c)
        else
            floored_leaves(c, acc)
        end
    end
    return acc
end

"""
$(TYPEDSIGNATURES)

The quantity of the composite that `node` stands for, as a variable of `m`:
`leaf(name, share, in_power)` for a leaf — the variable of that origin or
product, which the caller creates or looks up, told its share of the start
value and whether it is raised to a power — and otherwise a fresh variable
tied to its children by the CES aggregator of the node. `start` is the
expected quantity of the node, divided among its children by their shares,
and `base` names the constraints, numbered in the order the nodes are built.
A composite needs a positive floor only where it is itself raised to a power,
i.e. where the node it belongs to has a finite, positive elasticity.

Given `rhs`, the root gets no variable of its own: its aggregate, times
`scale`, is required to cover `rhs` directly, under the name `base`. That is a
nest of process inputs, whose composite is its coefficient times the activity.
"""
function ces_tree!(m, node, leaf, floor_at, start, base::String, in_power::Bool = false,
                   counter = Ref(0); rhs = nothing, scale = 1.0)
    node isa AbstractString && return leaf(node, start, in_power)
    qs = [ces_tree!(m, c, leaf, floor_at, s * start, base, powered(node.sigma), counter)
          for (c, s) in zip(node.children, node.shares)]
    if rhs !== nothing
        ces_constraint!(m, qs, node.shares, node.sigma, rhs, base; scale = scale)
        return rhs
    end
    at = sum(start_value, qs)
    composite = @variable(m, lower_bound = in_power ? floor_at(at) : 0.0, start = at,
                          base_name = "$(base)_composite[$(counter[] += 1)]")
    ces_constraint!(m, qs, node.shares, node.sigma, composite, "$base[$(counter[])]")
    return composite
end

"""
$(TYPEDSIGNATURES)

Resolve the [`Armington`](@ref) specifications of `d` into a
`(product, destination) => `[`Armington`](@ref) dictionary, covering every
destination of every product whose varieties are imperfect substitutes.
Products absent from it are homogeneous.

A product is Armington in every destination or in none, so a specification with
a finite elasticity for some destinations must be matched either by one for
each of the others or by a `:` fallback. [`validate`](@ref) runs this, so a
specification that breaks the rule is caught before anything is built.
"""
function armington_of(d::MarketData)
    specs    = Dict{Tuple{String,String},Armington}()   # (product, destination)
    fallback = Dict{String,Armington}()                 # product, from `:`
    products = String[]
    for a in d.armington
        isfinite(a.sigma) && push!(products, a.product)
        if a.destination isa Colon
            haskey(fallback, a.product) && throw(ArgumentError(
                "duplicate Armington specification for $(a.product)"))
            fallback[a.product] = a
        else
            for r in each(a.destination)
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
            "it: give one for every region, or one with destination = :"))
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

"""
$(TYPEDSIGNATURES)

Lay out the origins available to destination `r` as the [`CESNode`](@ref) tree
described by the groups of `a`: origins placed in a group are aggregated by it
first, and every other origin is a direct member of the root.  Each node's
shares are those of the origins beneath it, normalised among its siblings.
Groups that no available origin belongs to are dropped, and a group left with
a single member collapses into it. Returns the tree with the origins and
their shares.
"""
function armington_tree(d::MarketData, a::Armington, r::String,
                        can::Set{Tuple{String,String}} = producible(d))
    origins, shares = armington_origins(d, a, r, can)
    weight = Dict(zip(origins, shares))
    # weight of everything below a child, origins being the leaves
    subtotal(c) = c isa AbstractString ? weight[c] : sum(subtotal, c.children)

    # a group becomes a node over the available origins it holds, or `nothing`;
    # the shares are filled in once the whole tree stands
    function group(n::OriginNest)
        members = Any[o for o in n.origins if haskey(weight, o)]
        for sub in n.nests
            s = group(sub)
            s === nothing || push!(members, s)
        end
        return isempty(members) ? nothing : CESNode(n.sigma, members, Float64[])
    end
    function normalise(c)
        c isa CESNode || return c
        total = sum(subtotal, c.children)
        return CESNode(c.sigma, Any[normalise(k) for k in c.children],
                       Float64[subtotal(k) / total for k in c.children])
    end

    placed = leaves(a.nests)
    children = Any[]
    for n in a.nests
        s = group(n)
        s === nothing || push!(children, s)
    end
    for o in origins                     # origins in no group sit at the root
        o in placed || push!(children, o)
    end
    return collapse(normalise(CESNode(a.sigma, children, Float64[]))), origins, shares
end

"""
$(TYPEDSIGNATURES)

The quantities of the economy `d` that a [`Limit`](@ref) can name, as the model
will have them: the `(region, product)` of every supply curve and of every
product consumed, the `(region, process)` of every activity, the
`(region, process, nest, product)` of every input, and the
`(product, from, to)` of every flow between two regions — the routes of a
homogeneous product, and the origins of an Armington one that a destination
buys from.
"""
function model_catalogue(d::MarketData)
    arm = armington_of(d)
    arm_products = Set(p for (p, _) in keys(arm))
    can, uses = producible(d), consumed(d)
    flows = [k for k in keys(d.transport_costs) if !(k[1] in arm_products)]
    for ((p, r), a) in arm
        (p, r) in uses || continue
        _, origins, _ = armington_tree(d, a, r, can)
        append!(flows, [(p, o, r) for o in origins if o != r])
    end
    supply      = [(s.region, s.product) for s in d.supply]
    consumption = [(s.region, p) for s in d.demand for p in leaves(s.product)]
    processes   = [(r, proc) for proc in d.processes for r in regions_of(proc, d)]
    # the same, indexed by region, so that a limit naming one region does not
    # scan the whole economy
    by(f, xs) = (index = Dict{String,Vector{eltype(xs)}}();
                 foreach(x -> push!(get!(index, f(x), eltype(xs)[]), x), xs); index)
    return (supply = supply, consumption = consumption, processes = processes, flows = flows,
            supply_in = by(first, supply), consumption_in = by(first, consumption),
            processes_in = by(first, processes),
            flows_from = by(f -> f[2], flows), flows_to = by(f -> f[3], flows),
            order = Dict(r => i for (i, r) in enumerate(d.regions)))
end

# The entries of `all` worth looking at for the selector `sel`: those of its
# region only, through `index`, when it names a single region.
candidates(all, index, sel) =
    sel isa AbstractString ? get(index, sel, eltype(all)[]) : all

"""
$(TYPEDSIGNATURES)

The model variables that the quantity reference `q` selects, each with its
coefficient, keyed as `(:supply, region, product)`,
`(:consumption, region, product)`, `(:activity, region, process)`,
`(:input, region, process, nest, product)` or `(:flow, product, from, to)`.
Production is supply plus each process's activity times its yield; what a
fixed-coefficient (Leontief) input draws is its coefficient times the
activity, while a CES input has a variable of its own.
"""
function quantity_keys(q::QuantityRef, d::MarketData, cat)
    sel = q.selection
    is(f, x) = matches(getfield(sel, f), x)
    keys = Pair{Tuple,Float64}[]
    if q.kind == :supply || q.kind == :production
        for (r, p) in candidates(cat.supply, cat.supply_in, sel.region)
            is(:region, r) && is(:product, p) && push!(keys, (:supply, r, p) => 1.0)
        end
        q.kind == :production &&
            for (r, proc) in candidates(cat.processes, cat.processes_in, sel.region),
                (p, y) in proc.outputs
                is(:region, r) && is(:product, p) && push!(keys, (:activity, r, proc.name) => y)
            end
    elseif q.kind == :consumption
        for (r, p) in candidates(cat.consumption, cat.consumption_in, sel.region)
            is(:region, r) && is(:product, p) && push!(keys, (:consumption, r, p) => 1.0)
        end
    elseif q.kind == :activity
        for (r, proc) in candidates(cat.processes, cat.processes_in, sel.region)
            is(:region, r) && is(:process, proc.name) &&
                push!(keys, (:activity, r, proc.name) => 1.0)
        end
    elseif q.kind == :inputs
        for (r, proc) in candidates(cat.processes, cat.processes_in, sel.region),
            (k, nest) in enumerate(proc.inputs), p in nest.products
            is(:region, r) && is(:process, proc.name) && is(:product, p) || continue
            push!(keys, length(nest.products) == 1 ?
                        (:activity, r, proc.name) => nest.composite :
                        (:input, r, proc.name, k, p) => 1.0)
        end
    else
        pool = q.kind == :flow    ? candidates(cat.flows, cat.flows_from, sel.from) :
               q.kind == :exports ? candidates(cat.flows, cat.flows_from, sel.region) :
                                    candidates(cat.flows, cat.flows_to, sel.region)
        q.kind == :flow && pool === cat.flows && (pool = candidates(cat.flows, cat.flows_to, sel.to))
        for (p, o, r) in pool
            is(:product, p) || continue
            # for a group of regions only what crosses its border counts; with
            # no region named, the group is the world and every flow crosses
            every = q.kind != :flow && sel.region isa Colon
            keep = q.kind == :flow    ? is(:from, o) && is(:to, r) :
                   q.kind == :exports ? is(:region, o) && (every || !is(:region, r)) :
                                        is(:region, r) && (every || !is(:region, o))  # imports
            keep && push!(keys, (:flow, p, o, r) => 1.0)
        end
    end
    return keys
end

"""
$(TYPEDSIGNATURES)

The model variables of `limit`, each with its coefficient in the limit's
expression, the terms that select the same variable added up. A term that
selects nothing in the economy is an error.
"""
function limit_terms(limit::Limit, d::MarketData, cat = model_catalogue(d))
    terms = Dict{Tuple,Float64}()
    for (q, w) in limit.expression.terms
        keys = quantity_keys(q, d, cat)
        isempty(keys) && throw(ArgumentError(
            "limit $(limit.name): $q matches no quantity of this economy"))
        for (key, c) in keys
            terms[key] = get(terms, key, 0.0) + w * c
        end
    end
    filter!(t -> !iszero(last(t)), terms)
    return terms
end

describe_key(key::Tuple) =
    key[1] == :input ? "inputs[$(key[2]), $(key[3]), $(key[5])]" :
                       string(key[1], "[", join(key[2:end], ", "), "]")
regions_of_key(key::Tuple) = key[1] == :flow ? [key[3], key[4]] : [key[2]]

# The regions of the variables `keys`, once each, in the order the economy
# declares them.
function regions_in_order(keys, d::MarketData, order = Dict(r => i for (i, r) in enumerate(d.regions)))
    found = unique(r for k in keys for r in regions_of_key(k))
    return sort(found; by = r -> order[r])
end

# The products of the variables `keys`, once each, in the order the economy
# declares them. The activity of a process stands for what it makes, or, for
# one that makes nothing, such as a disposal, for what it takes in.
function products_in_order(keys, d::MarketData)
    process = Dict(p.name => p for p in d.processes)
    function of(key)
        key[1] == :flow  && return [key[2]]
        key[1] == :input && return [key[5]]
        key[1] == :activity || return [key[3]]
        proc = process[key[3]]
        return isempty(proc.outputs) ? [p for n in proc.inputs for p in n.products] :
                                       first.(proc.outputs)
    end
    order = Dict(p => i for (i, p) in enumerate(d.products))
    return sort(unique(p for k in keys for p in of(k)); by = p -> order[p])
end

# The kinds of the variables `keys`, once each: `supply`, `consumption`,
# `activity`, `inputs`, `flow`.
kinds_of(keys) = unique(k[1] == :input ? "inputs" : string(k[1]) for k in keys)

"""
    Markets.scope(limit, d)

What `limit` applies to in the economy `d`: the `variables` of the model its
expression involves, one by one — `supply[region, product]`,
`consumption[region, product]`, `activity[region, process]`,
`inputs[region, process, product]` and `flow[product, from, to]`, production
and fixed-coefficient inputs being expressed through the activities — and the
`regions` and `products` those variables belong to, in the order the economy
declares them, the activity of a process standing for what it makes.
`Results.limits` reports the same regions and products, and only the kinds of
the variables.
"""
function scope(limit::Limit, d::MarketData)
    terms = limit_terms(limit, d)
    keys = sort(collect(Base.keys(terms)); by = describe_key)
    return (regions = regions_in_order(keys, d), products = products_in_order(keys, d),
            variables = describe_key.(keys))
end

"""
$(TYPEDSIGNATURES)

The *relative* rate of tax on a final-demand curve.  An ad valorem rate scales
the price of every member of a group alike, and so scales the group's price
index; a rate that differed between members would not, and there is no single
curve left to put it on.  Per-unit rates have no such problem — they are
charged on each product's own quantity — so only this one has to be shared.
"""
function demand_relative_tax(d::MarketData, s::DemandCurve)
    products = leaves(s.product)
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

Build the JuMP model for `d`, with the ad valorem trade tariffs held at the
money equivalents given in `tariff`, and return it with the variables and
constraints the result tables are read from.

Internal: [`jump_model`](@ref) is what a user calls to see the model, and
[`solve_once`](@ref) what solves it.
"""
function build_model(d::MarketData, tariff::Dict{Tuple{String,String,String},Float64};
                     solver = Ipopt.Optimizer, solver_options = (;),
                     silent::Bool = true, qfloor::Union{Nothing,Real} = nothing)
    # The floor exists so that a variable raised to a power below 1 is never
    # differentiated at zero. What counts as "near zero" depends on the
    # variable: a millionth of what the data say that variable should be is
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
    demand_of = Dict{Tuple{String,String},DemandCurve}()
    for spec in d.demand, p in leaves(spec.product)
        demand_of[(spec.region, p)] = spec
    end

    # --- decision variables -------------------------------------------------
    # final demand quantities: floored because the marginal benefit a·D^(-1/η)
    # is unbounded at zero
    # Final demand. A plain product gets one variable; a group of substitutes
    # gets one per product plus the composite the curve actually describes, so
    # the balances and the result tables stay per product throughout.
    D = Dict{Tuple{String,String},VariableRef}()
    composites = Tuple{DemandCurve,Any}[]
    # a flat curve (η = Inf) or a straight one has a finite marginal benefit
    # at zero: no floor on what it is taken over
    for s in d.demand
        consumed_product(p, start, in_power) = begin
            f = floor_at(start)
            D[(s.region, p)] = @variable(m, lower_bound = in_power ? f : 0.0,
                                         start = max(start, f), base_name = "D[$(s.region),$p]")
        end
        node = s.product isa DemandNest ? tree(s.product) : s.product
        push!(composites,
              (s, ces_tree!(m, node, consumed_product, floor_at,
                            s.shift * reference_quantity(s), "demand_group[$(s.region)]",
                            s.form === :isoelastic && isfinite(s.elasticity))))
    end
    # primary supply quantities: the marginal cost b·S^(1/ε) is finite at zero,
    # so no floor is needed
    S = Dict{Tuple{String,String},VariableRef}()
    for s in d.supply
        S[(s.region, s.product)] = @variable(m, lower_bound = 0.0, start = s.shift * reference_quantity(s),
                                             base_name = "S[$(s.region),$(s.product)]")
    end
    # process activity levels
    z = Dict{Tuple{String,String},VariableRef}()
    for proc in d.processes, r in regions_of(proc, d)
        z[(r, proc.name)] = @variable(m, lower_bound = 0.0, start = 1.0,
                                      base_name = "z[$r,$(proc.name)]")
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
        floored = floored_leaves(root)
        # a rough but useful starting point: the local reference consumption
        a0 = get(demand_of, (r, p), nothing)
        start = a0 === nothing ? 1.0 : reference_quantity(a0)
        A[(p, r)] = @variable(m, lower_bound = 0.0, start = start, base_name = "A[$p,$r]")
        for (i, o) in pairs(origins)
            X[(p, o, r)] = @variable(m, lower_bound = o in floored ?
                                                      floor_at(shares[i] * start) : 0.0,
                                        start = shares[i] * start, base_name = "X[$p,$o,$r]")
        end
    end

    # --- assemble per-(region,product) sources & uses -----------------------
    produced = Dict{Tuple{String,String},Any}()   # process outputs + primary supply
    used     = Dict{Tuple{String,String},Any}()   # process inputs + final demand
    # what each nest of each process draws of each product, kept for the
    # physical accounts of the results: (region, process, nest, product)
    drawn    = Dict{Tuple{String,String,Int,String},Any}()
    addto!(dict, key, term) = (dict[key] = haskey(dict, key) ? dict[key] + term : term)

    for proc in d.processes, r in regions_of(proc, d)
        zk = z[(r, proc.name)]
        for (p, yield) in proc.outputs
            addto!(produced, (r, p), yield * zk)
        end
        for (k, nest) in enumerate(proc.inputs)
            if length(nest.products) == 1            # Leontief
                addto!(used, (r, nest.products[1]), nest.composite * zk)
                drawn[(r, proc.name, k, nest.products[1])] = nest.composite * zk
            else                                     # CES
                # one input variable per product, floored only where the
                # aggregator raises it to a power: a Leontief nest (σ = 0) or
                # perfect substitutes (σ = ∞) use their inputs linearly
                function input(p, start, in_power)
                    q = @variable(m, lower_bound = in_power ? floor_at(nest.composite) : 0.0,
                                     start = start, base_name = "x[$r,$(proc.name),$p]")
                    addto!(used, (r, p), q)
                    drawn[(r, proc.name, k, p)] = q
                    return q
                end
                # Calibrated CES: shares are *value* shares (Σδ = 1) and enter
                # the quantity aggregator as δ^(1/σ).  This makes the implied
                # composite-input price equal the input price when prices are
                # uniform (and total physical input = composite·z), so the dual
                # prices reproduce the genuine marginal cost of production.
                # The nest is a tree of one level, whose root is no variable
                # but the coefficient times the activity.
                ces_tree!(m, CESNode(nest.sigma, Any[nest.products...], nest.shares), input,
                          floor_at, nest.composite, "input_nest[$r,$(proc.name),$k]";
                          rhs = nest.composite * zk, scale = nest.phi)
            end
        end
    end

    # an expression that holds no variable, whether a number or an empty sum
    constant(x) = x isa Number || (x isa AbstractJuMPScalar && isempty(x.terms) && iszero(x.constant))

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
                origin_balance[(r, p)] = @constraint(m, src - ships == 0.0,
                                                        base_name = "origin_balance[$r,$p]")
            end
            # where the region uses the product, the composite covers those
            # uses and is the (possibly nested) CES aggregate of the varieties
            # bought: one composite variable and one constraint per group
            if haskey(A, (p, r))
                balance[(r, p)] = @constraint(m, A[(p, r)] - use == 0.0,
                                              base_name = "balance[$r,$p]")
                φ = scale_of(arm[(p, r)], r)
                # the leaves are the varieties, whose variables exist already
                aggregate = ces_tree!(m, tree_of[(p, r)], (o, _, _) -> X[(p, o, r)], floor_at,
                                      start_value(A[(p, r)]), "origin_group[$p,$r]")
                @constraint(m, (isone(φ) ? aggregate : φ * aggregate) >= A[(p, r)],
                            base_name = "armington[$p,$r]")
            end
        else
            # homogeneous product: imports into r minus exports from r
            imp = @expression(m, sum(T[(p, o, r)] for o in d.regions if haskey(T, (p, o, r)); init = 0.0))
            exp = @expression(m, sum(T[(p, r, o)] for o in d.regions if haskey(T, (p, r, o)); init = 0.0))
            # a product nothing in the region makes, uses or carries has no
            # balance there: `0 == 0` would only be a constraint on nothing,
            # which Ipopt counts against the degrees of freedom
            all(constant, (src, imp, use, exp)) && continue
            balance[(r, p)] = @constraint(m, src + imp - use - exp == 0.0,
                                          base_name = "balance[$r,$p]")
        end
    end

    # --- limits: linear constraints on the quantities ------------------------
    # A limit that is only an upper bound on one variable becomes a bound on it,
    # which the solver handles at no cost; every other one is a constraint, one
    # per finite bound, or an equality. Each keeps what is needed to report it:
    # the value of its expression, and for every bound the dual giving its rent.
    cat = model_catalogue(d)
    flow_vars = merge(Dict{Tuple{String,String,String},Any}(T),
                      Dict{Tuple{String,String,String},Any}(k => v for (k, v) in X if k[2] != k[3]))
    variable(key) = key[1] == :supply      ? S[(key[2], key[3])] :
                    key[1] == :consumption ? D[(key[2], key[3])] :
                    key[1] == :activity    ? z[(key[2], key[3])] :
                    key[1] == :input       ? drawn[(key[2], key[3], key[4], key[5])] :
                                             flow_vars[(key[2], key[3], key[4])]
    bounded = Set{Tuple}()
    limits = NamedTuple[]
    penalised = AffExpr(0.0)      # what crossing the soft limits costs
    for limit in d.limits
        terms = limit_terms(limit, d, cat)
        e = AffExpr(0.0)
        for (key, w) in terms
            add_to_expression!(e, w, variable(key))
        end
        base = "limit[$(limit.name)]"
        # each dual is (bound, reference, scale, side): side says which way the
        # bound faces, so that a soft limit crossed still counts as reached
        duals = Tuple{Float64,Any,Float64,Symbol}[]
        slacks = VariableRef[]
        soft = isfinite(limit.penalty)
        if !soft && length(terms) == 1 && only(values(terms)) > 0 && isinf(limit.lower) &&
           isfinite(limit.upper) && !(only(keys(terms)) in bounded)
            key, w = only(terms)
            limit.upper / w >= lower_bound(variable(key)) || throw(ArgumentError(
                "limit $(limit.name) caps $(describe_key(key)) at $(limit.upper / w), " *
                "below the floor of $(lower_bound(variable(key))) it must keep: a " *
                "quantity raised to a power cannot be held at zero"))
            set_upper_bound(variable(key), limit.upper / w)
            push!(bounded, key)
            push!(duals, (limit.upper, UpperBoundRef(variable(key)), 1 / w, :upper))
        elseif !soft && limit.lower == limit.upper
            push!(duals, (limit.upper, @constraint(m, e == limit.upper, base_name = base), 1.0, :equal))
        else
            # a soft bound is crossed by a slack that costs the penalty per unit
            function slack(side)
                soft || return 0.0
                x = @variable(m, lower_bound = 0.0, base_name = "excess[$(limit.name),$side]")
                push!(slacks, x)
                add_to_expression!(penalised, limit.penalty, x)
                return x
            end
            isfinite(limit.upper) &&
                push!(duals, (limit.upper, @constraint(m, e - slack(:upper) <= limit.upper,
                                                       base_name = base), 1.0, :upper))
            isfinite(limit.lower) &&
                push!(duals, (limit.lower, @constraint(m, e + slack(:lower) >= limit.lower,
                                                       base_name = base), 1.0, :lower))
        end
        push!(limits, (limit = limit, expression = e, duals = duals, terms = terms,
                       slacks = slacks))
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
    # at η = 1 the antiderivative is a·log(D). Along a straight line
    # D = shift·(α + β·p) the inverse demand is (D/shift - α)/β, whose integral
    # from zero is quadratic and concave, since β < 0.
    # A tax on consumption sits between what the consumer pays, which is what
    # the curve describes, and the market price, which the objective's gradient
    # has to equal: hence the division by (1 + relative) and the term in
    # `absolute` below.
    cons_benefit = @expression(m, sum(
        let w = 1 + demand_relative_tax(d, s),
            # per-unit rates are charged on each product's own quantity: the
            # composite is an index, and a tax per m³ is not a tax per index
            # point unless every member happens to cost the same
            taxed = [(tax_rates(d.taxes, "demand"; product = p, region = s.region).absolute, p)
                     for p in leaves(s.product)],
            benefit = curve_integral(s, C) / w
            # untaxed curves are left without the term, which would only print
            # as a string of zeros in `jump_model`
            all(iszero(first(t)) for t in taxed) ? benefit :
                benefit - sum(t / w * D[(s.region, p)] for (t, p) in taxed)
        end for (s, C) in composites))
    # A tax on production sits between the market price and what the producer
    # keeps, so it raises the marginal cost the market has to cover.
    supply_cost = @expression(m, sum(
        let s = supply_of[k],
            t = tax_rates(d.taxes, "supply"; product = s.product, region = s.region)
            (1 + t.relative) * curve_integral(s, S[k]) + t.absolute * S[k]
        end for k in keys(S)))
    conv_cost = @expression(m, sum(
        let t = tax_rates(d.taxes, "process"; process = proc.name, region = r)
            ((1 + t.relative) * vacost_of(proc, r) + t.absolute) * z[(r, proc.name)]
        end for proc in d.processes for r in regions_of(proc, d)))
    route_cost(key) = d.transport_costs[key] + get(tariff, key, 0.0) +
        tax_rates(d.taxes, "trade"; product = key[1], from = key[2], to = key[3]).absolute
    trans_cost = @expression(m,
        sum(route_cost(key) * T[key] for key in keys(T); init = 0.0) +
        sum(route_cost(key) * X[key] for key in keys(X) if key[2] != key[3]; init = 0.0))

    # crossing a soft limit costs its penalty: a device of the model, which
    # the results report apart from welfare
    @objective(m, Max, cons_benefit - supply_cost - conv_cost - trans_cost - penalised)

    # Welfare is the same surplus with every tax left out, a tax being a
    # transfer, and no penalty. Adding the revenue back to the objective would
    # recover it for per-unit rates only: a relative rate on supply or demand
    # scales a whole curve integral, not what the tax raises.
    surplus = @expression(m,
        sum(curve_integral(s, C) for (s, C) in composites; init = 0.0) -
        sum(curve_integral(supply_of[k], S[k]) for k in keys(S); init = 0.0) -
        sum(vacost_of(proc, r) * z[(r, proc.name)]
            for proc in d.processes for r in regions_of(proc, d); init = 0.0) -
        sum(d.transport_costs[key] * T[key] for key in keys(T); init = 0.0) -
        sum(d.transport_costs[key] * X[key] for key in keys(X) if key[2] != key[3]; init = 0.0))

    return m, (D = D, S = S, z = z, T = T, X = X, A = A, drawn = drawn, balance = balance,
               origin_balance = origin_balance, limits = limits, surplus = surplus)
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
                    tol::Real = 1e-3 * quantity_scale(d), model_options...)
    m, v = build_model(d, tariff; model_options...)
    wall_time = @elapsed optimize!(m)
    return build_results(d, m, v.D, v.S, v.z, v.T, v.X, v.balance, v.origin_balance,
                         v.limits, tariff, tol, wall_time; A = v.A, drawn = v.drawn,
                         surplus = v.surplus)
end

"""
$(TYPEDSIGNATURES)

The JuMP model that [`solve_market`](@ref) would solve for the economy `d`,
built but not solved: `print(jump_model(d))` shows its variables, objective
and constraints in full, and `latex_formulation(jump_model(d))` the same in
LaTeX.

!!! warning "Signs in the LaTeX"
    JuMP 1.32 drops the minus sign of a negative coefficient written in
    scientific notation when it typesets a model
    ([JuMP.jl#4268](https://github.com/jump-dev/JuMP.jl/issues/4268)), and
    the benefit terms of an inelastic demand are such coefficients. Read the
    signs from `print` until that is fixed.

`d` is validated first, as `solve_market` does, and the keywords are those of
`solve_market` that shape the model — `solver`, `solver_options`, `silent`
and `qfloor` — so the model comes back with its solver attached and can be
modified and solved with JuMP directly.

An **ad valorem tax on trade** is not in the objective: it is a fraction of a
price, which is a dual of the problem, and `solve_market` finds its money
equivalent by solving repeatedly. The model here is the first of those solves,
with those tariffs at zero. Every other tax is in it as it will be solved.
"""
function jump_model(d::MarketData; kwargs...)
    unknown = setdiff(keys(kwargs), (:solver, :solver_options, :silent, :qfloor))
    isempty(unknown) || throw(ArgumentError("jump_model does not take " *
        "$(join(unknown, ", ")): the model has no solution yet, so only the keywords " *
        "that shape it apply — solver, solver_options, silent and qfloor"))
    validate(d)
    return first(build_model(d, Dict{Tuple{String,String,String},Float64}(); kwargs...))
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
  By default each such quantity is floored at a millionth of its own reference
  magnitude; a number here sets one absolute floor for all of them instead.
* `tol` is the quantity below which a flow, an output or an activity is treated
  as zero and left out of the result tables. It defaults to a thousandth of the
  smallest reference quantity in `d`, so that it scales with the units the
  economy is written in.

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
        # the money a rate of `relative` comes to, at the prices just found:
        # the producer price at the origin, 0 where it has none (it then has
        # no variety to ship, so the route is empty whatever the tariff)
        at_origin = Dict((r, p) => coalesce(price, 0.0) for (r, p, price) in
                         zip(res.prices.region, res.prices.product, res.prices.producer_price))
        updated = Dict(key => tax_rates(d.taxes, "trade"; product = key[1], from = key[2],
                                        to = key[3]).relative * get(at_origin, (key[2], key[1]), 0.0)
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
