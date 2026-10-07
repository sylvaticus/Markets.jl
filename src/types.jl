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
    Unit transport costs, mapping `(product, from_region, to_region)` to the
    cost per unit of quantity shipped. A route that is missing from this
    dictionary cannot be used
    """
    transport::Dict{Tuple{Symbol,Symbol,Symbol},Float64} =
        Dict{Tuple{Symbol,Symbol,Symbol},Float64}()
end

regions_of(p::Process, d::MarketData) = p.regions === :all ? d.regions : p.regions
