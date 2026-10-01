# =============================================================================
#  Run the example forest-products market and print trade, consumption and
#  production for every world region.
#
#  From the package root:   julia --project=. examples/forest/run_example.jl
#  Or, in any environment where Markets and DataFrames are installed:
#                           include("examples/forest/run_example.jl")
# =============================================================================

using Markets
using DataFrames

include(joinpath(@__DIR__, "example_data.jl"))   # defines `example_market`

res = solve_market(example_market)

section(t) = (println("\n", "="^70); println("  ", t); println("="^70))

section("PRODUCTION  (primary supply + manufactured output, by region)")
show(unstack(res.production, :region, :product, :quantity); allcols = true, allrows = true)
println()

section("CONSUMPTION  (final demand, by region)")
show(unstack(res.consumption, :region, :product, :quantity); allcols = true, allrows = true)
println()

section("PRICES  (equilibrium shadow prices = duals of material balance)")
show(unstack(res.prices, :region, :product, :price); allcols = true, allrows = true)
println()

section("TRADE FLOWS  (positive shipments only)")
show(res.trade; allrows = true)
println()

section("NET TRADE  (exports − imports, by region & product)")
show(res.net_trade; allrows = true)
println()

# --- example of slicing everything for one region ---------------------------
section("EVERYTHING FOR ONE REGION  (Asia-Pacific :AS)")
as = for_region(res, :AS)
println("\nConsumption:");  show(as.consumption; allrows = true)
println("\n\nProduction:");  show(as.production;  allrows = true)
println("\n\nNet trade:");   show(as.net_trade;  allrows = true)
println("\n\nMill activity:"); show(as.activity; allrows = true)
println()
β