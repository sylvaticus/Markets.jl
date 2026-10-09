# ----------------------------------------------------------------------------
#  Data schema  (this is all a user touches to describe an economy)
#
#  All the types and the builders below take keyword arguments only.
# ----------------------------------------------------------------------------

"""
$(TYPEDEF)

Constant-elasticity *final* demand for one product in one region, calibrated to
pass through the reference point `(q0, p0)`.

Inverse demand is `P = a · D^(-1/η)` with `a = p0 · q0^(1/η)`.

# Fields
$(TYPEDFIELDS)

# Example
```julia
DemandSpec(product = :sawn_sw, region = :EU, p0 = 250, q0 = 90, elasticity = 1.3)
```
"""
Base.@kwdef struct DemandSpec
    "Consumed product"
    product::Symbol
    "Region where the product is consumed"
    region::Symbol
    "Reference price, i.e. the price at which the quantity demanded is `q0`"
    p0::Float64
    "Reference quantity, i.e. the quantity demanded at price `p0`"
    q0::Float64
    """
    Own-price demand elasticity ``\\eta`` (positive). **Use ``\\eta > 1``**, so
    that the consumer-surplus integral in the objective is finite
    """
    elasticity::Float64
end

"""
$(TYPEDEF)

Constant-elasticity *primary* supply (e.g. roundwood from the forest, crude oil
from wells, crops from land) for one product in one region, calibrated to pass
through the reference point `(q0, p0)`.

Inverse supply is `P = b · S^(1/ε)` with `b = p0 · q0^(-1/ε)`.

# Fields
$(TYPEDFIELDS)

# Example
```julia
SupplySpec(product = :swr, region = :NA, p0 = 65, q0 = 380, elasticity = 0.6)
```
"""
Base.@kwdef struct SupplySpec
    "Supplied product"
    product::Symbol
    "Region where the product is supplied"
    region::Symbol
    "Reference price, i.e. the price at which the quantity supplied is `q0`"
    p0::Float64
    "Reference quantity, i.e. the quantity supplied at price `p0`"
    q0::Float64
    "Own-price supply elasticity ``\\varepsilon`` (positive)"
    elasticity::Float64
end

"""
$(TYPEDEF)

One *input requirement* of a [`Process`](@ref). Inputs in the same nest
substitute for each other; different nests of the same process are needed in
fixed proportion to the process activity.

* A **Leontief** (fixed-coefficient) requirement is a nest with a single
  product: `composite` units of it are needed per unit of process activity.
  Build it with [`leontief`](@ref).
* A **CES** (smoothly substitutable) requirement bundles several `products`
  into a composite via a CES aggregator with value `shares` (δ, Σδ = 1),
  elasticity of substitution `sigma` (σ) and scale `phi` (φ):

      composite · z  ≤  φ · ( Σᵢ δᵢ^(1/σ) · qᵢ^ρ )^(1/ρ),     ρ = (σ-1)/σ

  As an input's price rises the optimiser smoothly shifts the mix toward
  cheaper inputs. Build it with [`ces`](@ref), which requires `sigma > 1`.

Use the [`leontief`](@ref) and [`ces`](@ref) builders rather than this
constructor: they normalise the shares and check the arguments.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct Nest
    "Units of the (composite) input needed per unit of process activity"
    composite::Float64
    "Products of the nest. A single product makes the nest a Leontief one"
    products::Vector{Symbol}
    """
    Value shares ``\\delta`` of each product, summing to 1 (equal shares by
    default). They set the input mix when all input prices are equal. Ignored
    for a single-product (Leontief) nest
    """
    shares::Vector{Float64} = fill(1 / length(products), length(products))
    """
    Elasticity of substitution ``\\sigma`` between the products of the nest.
    It must be `> 1`; it is ignored (and conventionally set to 0) for a
    single-product (Leontief) nest
    """
    sigma::Float64
    "Scale parameter ``\\phi`` of the CES aggregator"
    phi::Float64 = 1.0
end

"""
$(TYPEDSIGNATURES)

Build a fixed-coefficient (Leontief) input requirement: `coeff` units of
`product` per unit of process activity.

# Example
```julia
leontief(product = :pulp, coeff = 1.10)    # 1.1 t of pulp per t of paper
```
"""
leontief(; product::Symbol, coeff::Real) =
    Nest(composite = coeff, products = [product], shares = [1.0], sigma = 0.0, phi = 1.0)

"""
$(TYPEDSIGNATURES)

Build a smoothly substitutable (CES) input requirement of `composite` units of
a bundle of `products` per unit of process activity.

`shares` are the value shares of each input when all input prices are equal
(equal shares by default); they need not sum to 1, as they are normalised
internally. `sigma` is the elasticity of substitution and must be `> 1`.

# Example
```julia
ces(composite = 1.30, products = [:chips, :swr, :hwr],
    shares = [0.55, 0.30, 0.15], sigma = 2.5)
```
"""
function ces(; composite::Real, products::Vector{Symbol},
               shares::Vector{<:Real} = fill(1.0, length(products)),
               sigma::Real, phi::Real = 1.0)
    @assert length(products) == length(shares) "products and shares must align"
    @assert sigma > 1 "use σ > 1 for a well-behaved (smooth) CES nest"
    s = collect(float.(shares)); s ./= sum(s)
    Nest(composite = composite, products = products, shares = s, sigma = sigma, phi = phi)
end

"""
$(TYPEDEF)

A transformation activity, operated at a level of activity chosen by the model.

# Fields
$(TYPEDFIELDS)

# Example
```julia
Process(name    = :sawmill_sw,
        inputs  = [leontief(product = :swr, coeff = 1.0)],
        outputs = [:sawn_sw => 0.50, :chips => 0.35],
        vacost  = 40)
```
"""
Base.@kwdef struct Process
    "Name of the process"
    name::Symbol
    """
    Input requirements, as a vector of [`Nest`](@ref)s built with
    [`leontief`](@ref) or [`ces`](@ref). Each nest is needed in fixed
    proportion to the process activity
    """
    inputs::Vector{Nest}
    """
    Outputs, as `product => yield` pairs per unit of activity. Several pairs
    make a joint-production process (e.g. sawnwood *and* chips)
    """
    outputs::Vector{Pair{Symbol,Float64}}
    "Regions where the process may operate, or `:all` for every region"
    regions::Union{Symbol,Vector{Symbol}} = :all
    """
    Value-added cost per unit of activity: the cost of everything not modelled
    as an explicit input (labour, energy, capital, ...)
    """
    vacost::Float64 = 0.0
end

"""
$(TYPEDEF)

A group of origins that substitute for each other more closely than they do for
the rest, inside an [`Armington`](@ref) composite.

Regions whose products are technically interchangeable — because they share
grades, standards or certification, or because the same non-tariff measures
let them through — belong in one group. Within it buyers substitute with the
group's own `sigma`; against anything outside it they substitute with the
`sigma` of the enclosing nest, which is lower.

Groups may be nested further through their own `nests` field, so the structure
is a tree whose leaves are origins. This is the general shape an elasticity
specification may take: the substitution between any two origins is the
elasticity of the smallest group containing both, and a group must be at least
as substitutable inside as it is with the outside.

# Fields
$(TYPEDFIELDS)

# Example
```julia
# in this market, EU and NA products are near-interchangeable
OriginNest(sigma = 12, origins = [:EU, :NA])
```
"""
Base.@kwdef struct OriginNest
    """
    Elasticity of substitution ``\\sigma`` between the members of the group. It
    must be `> 1`, and at least the `sigma` of the nest that encloses it.
    `Inf` makes the members perfect substitutes: the group is pooled into a
    single good, as in the homogeneous case, while still being a distinct
    variety to the outside
    """
    sigma::Float64
    "Origins belonging to the group"
    origins::Vector{Symbol} = Symbol[]
    "Sub-groups of the group, for a deeper tree"
    nests::Vector{OriginNest} = OriginNest[]

    function OriginNest(sigma, origins, nests)
        sigma > 1 || throw(ArgumentError(
            "an origin group's elasticity must be > 1, got $sigma"))
        isempty(origins) && isempty(nests) && throw(ArgumentError(
            "an origin group must contain at least one origin or sub-group"))
        for n in nests
            n.sigma ≥ sigma || throw(ArgumentError(
                "the sub-group of elasticity $(n.sigma) is less substitutable " *
                "inside than with the outside ($sigma): a group's members must " *
                "substitute for each other at least as easily as for non-members"))
        end
        new(sigma, origins, nests)
    end
end

"""
$(TYPEDEF)

Imperfect substitution between the regional varieties of one product
(Armington, 1969).

By default a product is **homogeneous**: it is the same good wherever it comes
from, regional prices differ by at most the transport cost, and a region never
both imports and exports it. Listing a product here instead makes the varieties
from the different origins imperfect substitutes: each region uses a CES
composite of them, so a region can import and export the same product at once
and its price responds to supply and demand everywhere, not only to the cost of
the cheapest source.

Each destination `r` draws on the varieties of the origins it can buy from: its
own, plus every origin `o` with a transport route `(product, o, r)`. The
composite is

    A ≤ ( Σₒ δₒ^(1/σ) · Xₒ^ρ )^(1/ρ),     ρ = (σ-1)/σ

with `Xₒ` the quantity bought from origin `o` and `δₒ` its value share. The
price a region pays for the composite is then the CES price index of the
delivered prices,

    Pᶜ = ( Σₒ δₒ · (Pₒ + τₒ)^(1-σ) )^(1/(1-σ))

which tends to the cheapest delivered price as `σ → ∞`: `sigma = Inf` is
exactly the homogeneous (Samuelson spatial price equilibrium) case, and a large
finite `sigma` approaches it.

A single `sigma` makes every origin substitute equally well for every other,
which is the defining restriction of a flat CES. Two ways out, both available
here:

* **Groups of origins** ([`OriginNest`](@ref)): origins whose products are
  mutually interchangeable are nested together with a higher `sigma`, and
  substitute with the rest at the lower `sigma` of this specification.
* **One specification per destination** (`destination`): the structure may
  differ from one importing market to the next, which is what technical
  standards and non-tariff measures actually do, since the importer sets them.

# Fields
$(TYPEDFIELDS)

# Example
```julia
# moderate substitution, equal shares, same in every market
Armington(product = :paper, sigma = 4)

# calibrated shares: 70% domestic, the rest split between the two other regions
Armington(product = :paper, sigma = 4,
          shares = merge(Dict((r, r) => 0.7 for r in regions),
                         Dict((o, r) => 0.15 for o in regions, r in regions if o != r)))

# in the EU market, EU and NA sawnwood are near-interchangeable (shared grading
# rules) while AS sawnwood is a distant substitute; elsewhere everything
# substitutes freely
Armington(product = :sawn_sw, destination = :EU, sigma = 2.5,
          nests = [OriginNest(sigma = 12, origins = [:EU, :NA])])
Armington(product = :sawn_sw, destination = [:NA, :AS], sigma = 8)
```
"""
Base.@kwdef struct Armington
    "Product whose regional varieties are imperfect substitutes"
    product::Symbol
    """
    Elasticity of substitution ``\\sigma`` between the varieties of the
    different origins. It must be `> 1`; `Inf` (the default for every product
    not listed) means perfect substitutes, i.e. a homogeneous product
    """
    sigma::Float64 = Inf
    """
    Value shares ``\\delta`` of each origin in each destination's composite,
    keyed `(origin, destination)` — the shares that would be observed if the
    delivered prices of all origins were equal. They are normalised per
    destination. An origin with no share given (or a share of 0) is left out of
    that destination's composite. When this is empty, every origin available to
    a destination gets an equal share
    """
    shares::Dict{Tuple{Symbol,Symbol},Float64} = Dict{Tuple{Symbol,Symbol},Float64}()
    """
    Importing regions this specification applies to, or `:all` (the default)
    for every region without one of their own. A product may have one
    specification per destination plus an `:all` fallback, which is how the
    structure, the elasticity and the shares can differ between importing
    markets
    """
    destination::Union{Symbol,Vector{Symbol}} = :all
    """
    Groups of origins that substitute for each other more closely than for the
    rest, as a vector of [`OriginNest`](@ref)s. Every available origin not
    placed in a group stays a direct member of the composite. Empty (the
    default) gives the flat CES in which all origins substitute equally
    """
    nests::Vector{OriginNest} = OriginNest[]

    function Armington(product, sigma, shares, destination, nests)
        sigma > 1 || throw(ArgumentError(
            "the Armington elasticity of $product must be > 1 (or Inf for a " *
            "homogeneous product), got $sigma"))
        all(≥(0), values(shares)) || throw(ArgumentError(
            "the Armington shares of $product must not be negative"))
        for n in nests
            n.sigma ≥ sigma || throw(ArgumentError(
                "a group of origins of $product substitutes less inside " *
                "($(n.sigma)) than with the outside ($sigma): a group's members " *
                "must substitute for each other at least as easily as for non-members"))
        end
        grouped = Symbol[]
        collect_origins!(grouped, nests)
        allunique(grouped) || throw(ArgumentError(
            "an origin appears in more than one group of $product"))
        new(product, sigma, shares, destination, nests)
    end
end

"Append every origin appearing anywhere in `nests` to `acc`."
function collect_origins!(acc::Vector{Symbol}, nests::Vector{OriginNest})
    for n in nests
        append!(acc, n.origins)
        collect_origins!(acc, n.nests)
    end
    return acc
end

"""
$(TYPEDEF)

The complete description of an economy, to be passed to [`solve_market`](@ref).

# Fields
$(TYPEDFIELDS)

# Example
```julia
MarketData(regions   = [:EU, :NA],
           products  = [:swr, :sawn_sw],
           demand    = demand,      # vector of DemandSpec
           supply    = supply,      # vector of SupplySpec
           processes = processes,   # vector of Process
           tradable  = [:swr, :sawn_sw],
           transport = Dict((:swr, :EU, :NA) => 12.0))
```
"""
Base.@kwdef struct MarketData
    "Regions of the model"
    regions::Vector{Symbol}
    "Products of the model"
    products::Vector{Symbol}
    "Final demand curves, one [`DemandSpec`](@ref) per consumed product and region"
    demand::Vector{DemandSpec} = DemandSpec[]
    "Primary supply curves, one [`SupplySpec`](@ref) per supplied product and region"
    supply::Vector{SupplySpec} = SupplySpec[]
    "Transformation activities, a vector of [`Process`](@ref)es"
    processes::Vector{Process} = Process[]
    "Products that may be traded between regions"
    tradable::Vector{Symbol} = Symbol[]
    """
    Products whose regional varieties are imperfect substitutes, as a vector of
    [`Armington`](@ref) specifications. A product that is not listed here (the
    default) is homogeneous: the same good wherever it comes from
    """
    armington::Vector{Armington} = Armington[]
    """
    Unit transport costs, mapping `(product, from_region, to_region)` to the
    cost per unit of quantity shipped. A route that is missing from this
    dictionary cannot be used
    """
    transport::Dict{Tuple{Symbol,Symbol,Symbol},Float64} =
        Dict{Tuple{Symbol,Symbol,Symbol},Float64}()
end

regions_of(p::Process, d::MarketData) = p.regions === :all ? d.regions : p.regions
