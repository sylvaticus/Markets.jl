# ----------------------------------------------------------------------------
#  Data schema  (this is all a user touches to describe an economy)
# ----------------------------------------------------------------------------

"""
    DemandSpec(product, region; p0, q0, elasticity)

Constant-elasticity *final* demand for `product` in `region`, calibrated to pass
through the reference point `(q0, p0)`.  Inverse demand is `P = a · D^(-1/η)`
with `a = p0 · q0^(1/η)`.  `elasticity` is the (positive) own-price demand
elasticity `η`; **use η > 1** so the consumer-surplus integral is finite.
"""
struct DemandSpec
    product::Symbol
    region::Symbol
    p0::Float64
    q0::Float64
    elasticity::Float64
end
DemandSpec(p, r; p0, q0, elasticity) = DemandSpec(p, r, p0, q0, elasticity)

"""
    SupplySpec(product, region; p0, q0, elasticity)

Constant-elasticity *primary* supply (e.g. roundwood from the forest, crude
oil from wells) for `product` in `region`.  Inverse supply is `P = b · S^(1/ε)`
with `b = p0 · q0^(-1/ε)`.  `elasticity` is the (positive) supply elasticity `ε`.
"""
struct SupplySpec
    product::Symbol
    region::Symbol
    p0::Float64
    q0::Float64
    elasticity::Float64
end
SupplySpec(p, r; p0, q0, elasticity) = SupplySpec(p, r, p0, q0, elasticity)

"""
    Nest(composite, products, shares, sigma, phi)

One *input requirement* of a process.  Inputs in the same nest substitute for
each other; different nests of the same process are needed in fixed proportion
to the process activity.

* A **Leontief** (fixed-coefficient) requirement is a nest with a single
  product: `composite` units of it are needed per unit of process activity.
  Build it with [`leontief`](@ref).
* A **CES** (smoothly substitutable) requirement bundles several `products`
  into a composite via a CES aggregator with value `shares` (δ, Σδ = 1),
  elasticity of substitution `sigma` (σ) and scale `phi` (φ):

      composite · z  ≤  φ · ( Σᵢ δᵢ^(1/σ) · qᵢ^ρ )^(1/ρ),     ρ = (σ-1)/σ

  As an input's price rises the optimiser smoothly shifts the mix toward
  cheaper inputs.  Build it with [`ces`](@ref), which requires `sigma > 1`.

The engine distinguishes the two cases by the number of products in the nest,
so the `sigma` stored by [`leontief`](@ref) (0) is never used.
"""
struct Nest
    composite::Float64
    products::Vector{Symbol}
    shares::Vector{Float64}
    sigma::Float64
    phi::Float64
end

"""
    leontief(product, coeff)

A fixed-coefficient input: `coeff` units of `product` per unit of activity.
"""
leontief(product::Symbol, coeff::Real) = Nest(float(coeff), [product], [1.0], 0.0, 1.0)

"""
    ces(composite, products, shares; sigma, phi=1.0)

A smoothly substitutable (CES) input requirement of `composite` units per unit
of activity.  `shares` are the value shares of each input at equal input prices;
they need not sum to 1, as they are normalised internally.  `sigma` is the
elasticity of substitution and must be `> 1`.
"""
function ces(composite::Real, products::Vector{Symbol}, shares::Vector{<:Real};
             sigma::Real, phi::Real = 1.0)
    @assert length(products) == length(shares) "products and shares must align"
    @assert sigma > 1 "use σ > 1 for a well-behaved (smooth) CES nest"
    s = collect(float.(shares)); s ./= sum(s)
    Nest(float(composite), products, s, float(sigma), float(phi))
end

"""
    Process(name, inputs, outputs; regions=:all, vacost=0.0)

A transformation activity.  `inputs` is a vector of [`Nest`](@ref)s; `outputs`
is a vector of `product => yield` pairs per unit of activity (joint products are
allowed, e.g. sawmilling yields sawnwood *and* chips).  `vacost` is the
value-added (conversion) cost per unit of activity, i.e. the cost of everything
not modelled as an explicit input (labour, energy, capital, ...).  `regions`
lists where the process may operate (`:all` for every region).
"""
struct Process
    name::Symbol
    inputs::Vector{Nest}
    outputs::Vector{Pair{Symbol,Float64}}
    regions::Union{Symbol,Vector{Symbol}}
    vacost::Float64
end
Process(name, inputs, outputs; regions = :all, vacost = 0.0) =
    Process(name, inputs, Pair{Symbol,Float64}[p => float(v) for (p, v) in outputs], regions, vacost)

"""
    MarketData(; regions, products, demand, supply, processes, tradable, transport)

The complete description of an economy.  `transport` maps
`(product, from_region, to_region) => unit_cost` (cost per unit of quantity
shipped); missing pairs are not tradable on that route.  `tradable` lists which
products may be traded at all.
"""
struct MarketData
    regions::Vector{Symbol}
    products::Vector{Symbol}
    demand::Vector{DemandSpec}
    supply::Vector{SupplySpec}
    processes::Vector{Process}
    tradable::Vector{Symbol}
    transport::Dict{Tuple{Symbol,Symbol,Symbol},Float64}
end
MarketData(; regions, products, demand, supply, processes, tradable, transport) =
    MarketData(regions, products, demand, supply, processes, tradable, transport)

regions_of(p::Process, d::MarketData) = p.regions === :all ? d.regions : p.regions
