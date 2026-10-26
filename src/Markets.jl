"""
    Markets

Solve a spatial partial-equilibrium model of an economic sector made of an
arbitrary set of products, regions and transformation processes.

The package is the *engine* only: it knows nothing about specific products.
You describe an economy with the data types [`MarketData`](@ref),
[`DemandSpec`](@ref), [`SupplySpec`](@ref) and [`Process`](@ref) (built from
[`leontief`](@ref) and [`ces`](@ref) input nests), hand it to
[`solve_market`](@ref), and get back production, consumption, trade and prices
per region in a [`Results`](@ref) object.

The engine implements a **spatial price equilibrium** (Samuelson 1952 /
Takayama & Judge 1971) solved as a **net-social-surplus maximisation**.  Under
perfect competition this single optimisation *is* the market equilibrium:

* the **quantities** (production, consumption, trade) are the primal solution;
* the **prices** are the dual variables (shadow prices) of the per-region,
  per-product material-balance constraints.

Smooth (CES) input substitution, multi-stage transformation with joint products
and residues, and trade between non-identical regional markets (transport
costs) all fall out of the same formulation.

A complete example (the forest-products sector) is in `examples/forest/`.
"""
module Markets

using JuMP, Ipopt, DataFrames
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

export MarketData, DemandSpec, SupplySpec, Process, Nest, Armington, OriginNest,
       TradeLimit, DemandNest, Tax,
       leontief, ces, disposal, armington_shares, armington_calibration, ces_calibration,
       impute_flows, exogenous_shift,
       validate, solve_market, jump_model, Results, for_region

# re-exported from JuMP, so that a model from `jump_model` can be typeset
# without `using JuMP`
export latex_formulation

include("names.jl")
include("types.jl")
include("validate.jl")
include("results.jl")     # defines Results, which model.jl dispatches on
include("model.jl")
include("show.jl")

end # module
