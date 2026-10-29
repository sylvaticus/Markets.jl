# Part of Markets.jl package [https://github.com/sylvaticus/Markets.jl] - MIT Licence

# ----------------------------------------------------------------------------
#  Material accounts: what the goods of a solved economy are made of
# ----------------------------------------------------------------------------

"""
$(TYPEDEF)

A material carried by the goods of an economy — carbon, dry matter, a metal,
water — whose flows [`material_accounts`](@ref) traces through a solved
economy, from the primary supply it enters with to the goods it ends up in.

The material enters with primary supply, at a given `content` per unit. It
travels with the goods: a region's own production carries the content of what
made it, and what a region uses carries the average content of what it
received. In a process, part of what enters may be lost — burnt, emitted,
discarded — and the rest leaves in the outputs, divided among them by
`allocation`. Nothing else creates or destroys it, so what was extracted is
what was consumed plus what was lost.

How much of the material a unit of a manufactured good carries is not given
but follows from the solution: a mill that changes its mix of inputs changes
the content of its output. This is the material-intensity accounting of
[Lorang & Gerlagh (2026)](@ref lorang2026), applied after the solve: the
equilibrium does not depend on it.

# Fields
$(TYPEDFIELDS)

# Example
```julia
carbon = Material(name = "carbon",
                  content = Dict("logs" => 0.25),            # t C per m³ harvested
                  losses  = Dict("pulpmill" => 0.5),         # half burnt for energy
                  allocation = Dict("sawmill" => Dict("lumber" => 1.0, "chips" => 1.0)))
```
"""
Base.@kwdef struct Material
    "Name of the material, used to report it"
    name::String
    """
    Content per unit of primary supply, keyed by product — the same in every
    region — or by `(region, product)` where it differs, the more specific key
    taking precedence. A supplied product not listed carries none
    """
    content::Dict{Union{String,Tuple{String,String}},Float64}
    """
    Share of the material entering each process that leaves it other than in
    its outputs, keyed by process name; 0 for a process not listed. A process
    with no outputs, such as a [`disposal`](@ref), loses all of it
    """
    losses::Dict{String,Float64} = Dict{String,Float64}()
    """
    How a process with several outputs divides the material leaving it, keyed
    by process name: a weight per unit of each output, as `product => weight`.
    By default every unit of every output gets the same content
    """
    allocation::Dict{String,Dict{String,Float64}} = Dict{String,Dict{String,Float64}}()

    function Material(name, content, losses, allocation)
        key(k) = k isa Tuple ? (nm(k[1]), nm(k[2])) : nm(k)
        c = Dict{Union{String,Tuple{String,String}},Float64}(key(k) => float(v) for (k, v) in content)
        all(>=(0), values(c)) || throw(ArgumentError(
            "the content of $name must not be negative"))
        l = namekeys(losses)
        all(x -> 0 <= x <= 1, values(l)) || throw(ArgumentError(
            "the losses of $name are shares, between 0 and 1"))
        a = Dict{String,Dict{String,Float64}}(nm(k) => namekeys(v) for (k, v) in allocation)
        all(w >= 0 for v in values(a) for w in values(v)) || throw(ArgumentError(
            "the allocation weights of $name must not be negative"))
        new(nm(name), c, l, a)
    end
end

content_of(m::Material, r::String, p::String) =
    get(m.content, (r, p), get(m.content, p, 0.0))

"""
$(TYPEDSIGNATURES)

Trace `material` through the solved economy `res`, and return three tables:

* `intensities`: `region`, `product`, the content of a unit of what the region
  `produced`, and of a unit of what its users `used` — the average of what it
  received, for an [`Armington`](@ref) product, or of its own production and
  imports, for a homogeneous one. `missing` where there is none.
* `balances`: by `region` and `product`, the material `extracted` with primary
  supply, `made` into the product by processes, `imported`, `exported`,
  `consumed` in final demand and `processed` as an input. For every row,
  `extracted + made + imported − exported = consumed + processed`.
* `processes`: by `region` and `process`, the material `input`, `lost` and
  `output`.

The `consumed` column is the material embodied in each region's final
consumption: its footprint. Summed over the world, what was extracted is what
was consumed plus what was lost.

The content of a manufactured good follows from its inputs: a process's
output carries what entered it, less its losses, divided among its outputs by
the allocation weights. These intensities form a linear system, solved by
iteration from zero, which converges as long as the material cannot circulate
for ever without loss; `tol` is the relative precision asked of them.

The accounts are built from `res.solution`, so they close exactly, and do not
change the equilibrium: they describe it.
"""
function material_accounts(res::Results, material::Material; tol::Real = 1e-12,
                           maxiter::Integer = 100_000)
    d, sol = res.data, res.solution
    for (k, _) in material.content
        r, p = k isa Tuple ? k : (nothing, k)
        r === nothing || r in d.regions || throw(ArgumentError(
            "the content of $(material.name) names the region $r, which the economy does not declare"))
        p in d.products || throw(ArgumentError(
            "the content of $(material.name) names the product $p, which the economy does not declare"))
    end
    processes = Dict(proc.name => proc for proc in d.processes)
    for k in keys(material.losses)
        haskey(processes, k) || throw(ArgumentError(
            "the losses of $(material.name) name the process $k, which the economy does not declare"))
    end
    for (k, w) in material.allocation
        haskey(processes, k) || throw(ArgumentError(
            "the allocation of $(material.name) names the process $k, which the economy does not declare"))
        made = Set(first.(processes[k].outputs))
        for p in keys(w)
            p in made || throw(ArgumentError(
                "the allocation of $(material.name) gives a weight to $p, which process $k does not make"))
        end
    end

    arm = Set(p for (p, _) in keys(armington_of(d)))
    supplied(r, p) = get(sol.supply, (r, p), 0.0)
    level(r, k) = get(sol.activity, (r, k), 0.0)
    production = sol.production
    # what each region receives of each product, physically: the varieties of
    # an Armington product, its own included, or the imports of a homogeneous one
    into = Dict{Tuple{String,String},Vector{Tuple{String,Float64}}}()
    for ((p, o, r), x) in sol.flows
        x > 0 && push!(get!(into, (r, p), Tuple{String,Float64}[]), (o, x))
    end
    sent(r, p) = sum(x for ((q, o, s), x) in sol.flows if q == p && o == r && s != r; init = 0.0)
    received(r, p) = sum(last, get(into, (r, p), Tuple{String,Float64}[]); init = 0.0)

    # the unknowns: the content of a unit of each region's production, θ, and
    # of a unit of what its users draw on (a unit of composite, for an
    # Armington product), μ
    θ = Dict{Tuple{String,String},Float64}()
    μ = Dict{Tuple{String,String},Float64}()
    weight(proc, p) = get(get(material.allocation, proc.name, Dict{String,Float64}()), p, 1.0)
    weights(proc) = sum((weight(proc, p) * y for (p, y) in proc.outputs); init = 0.0)
    # with no output to carry it, everything that enters a process is lost
    loss(proc) = weights(proc) > 0 ? get(material.losses, proc.name, 0.0) : 1.0
    material_in(r, proc) = sum(x * get(μ, (r, p), 0.0)
                               for ((r2, k, _, p), x) in sol.drawn if r2 == r && k == proc.name;
                               init = 0.0)
    function sweep!()
        change = 0.0
        update!(dict, key, value) = (change = max(change, abs(value - get(dict, key, 0.0)));
                                     dict[key] = value)
        for r in d.regions, p in d.products
            if p in arm
                A = get(sol.composite, (r, p), 0.0)
                A > 0 || continue
                update!(μ, (r, p), sum(x * get(θ, (o, p), 0.0)
                                       for (o, x) in get(into, (r, p), Tuple{String,Float64}[]);
                                       init = 0.0) / A)
            else
                pool = get(production, (r, p), 0.0) + received(r, p)
                pool > 0 || continue
                update!(μ, (r, p), (get(θ, (r, p), 0.0) * get(production, (r, p), 0.0) +
                                    sum(x * get(μ, (o, p), 0.0)
                                        for (o, x) in get(into, (r, p), Tuple{String,Float64}[]);
                                        init = 0.0)) / pool)
            end
        end
        made = Dict{Tuple{String,String},Float64}()
        for proc in d.processes, r in regions_of(proc, d)
            total = weights(proc)
            total > 0 || continue
            out = (1 - loss(proc)) * material_in(r, proc)
            for (p, y) in proc.outputs
                made[(r, p)] = get(made, (r, p), 0.0) + out * weight(proc, p) * y / total
            end
        end
        for ((r, p), q) in production
            q > 0 || continue
            update!(θ, (r, p), (content_of(material, r, p) * supplied(r, p) +
                                get(made, (r, p), 0.0)) / q)
        end
        return change
    end
    scale_of_values() = max(1.0, maximum(values(θ); init = 0.0), maximum(values(μ); init = 0.0))
    settled = false
    for _ in 1:maxiter
        if sweep!() <= tol * scale_of_values()
            settled = true
            break
        end
    end
    settled || throw(ErrorException(
        "the content of $(material.name) did not settle in $maxiter rounds: it circulates " *
        "through a loop of processes or trade that loses none of it"))

    # --- the tables -------------------------------------------------------------
    intensities = DataFrame(region = String[], product = String[],
                            produced = Union{Float64,Missing}[], used = Union{Float64,Missing}[])
    for r in d.regions, p in d.products
        prod = get(production, (r, p), 0.0) > 0 ? θ[(r, p)] : missing
        use = if p in arm
            got = received(r, p)
            got > 0 && haskey(μ, (r, p)) ? μ[(r, p)] * sol.composite[(r, p)] / got : missing
        else
            get(μ, (r, p), missing)
        end
        (ismissing(prod) && ismissing(use)) || push!(intensities, (r, p, prod, use))
    end

    balances = DataFrame(region = String[], product = String[], extracted = Float64[],
                         made = Float64[], imported = Float64[], exported = Float64[],
                         consumed = Float64[], processed = Float64[])
    for r in d.regions, p in d.products
        extracted = content_of(material, r, p) * supplied(r, p)
        made = get(production, (r, p), 0.0) > 0 ?
               θ[(r, p)] * get(production, (r, p), 0.0) - extracted : 0.0
        carried(o) = p in arm ? get(θ, (o, p), 0.0) : get(μ, (o, p), 0.0)
        imported = sum(x * carried(o) for (o, x) in get(into, (r, p), Tuple{String,Float64}[])
                       if o != r; init = 0.0)
        exported = sent(r, p) * carried(r)
        consumed = get(sol.consumption, (r, p), 0.0) * get(μ, (r, p), 0.0)
        processed = sum(x for ((r2, _, _, q), x) in sol.drawn if r2 == r && q == p;
                        init = 0.0) * get(μ, (r, p), 0.0)
        row = (r, p, extracted, made, imported, exported, consumed, processed)
        any(!iszero, row[3:end]) && push!(balances, row)
    end

    flows_through = DataFrame(region = String[], process = String[], input = Float64[],
                              lost = Float64[], output = Float64[])
    for proc in d.processes, r in regions_of(proc, d)
        level(r, proc.name) > 0 || continue
        input = material_in(r, proc)
        lost = loss(proc) * input
        push!(flows_through, (r, proc.name, input, lost, input - lost))
    end
    return (intensities = intensities, balances = balances, processes = flows_through)
end
