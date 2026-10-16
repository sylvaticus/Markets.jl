# Markets.jl

!!! warning "Pre-release: under active development and validation"
    This project is at a _proof-of-concept_ stage and has not been released. A substantial part of it was generated with the assistance of an AI model and is still being reviewed.
    Both the modelling framework (its assumptions and equations) and its software implementation are currently undergoing verification and validation.

    - Results have not yet been validated and should not be used for publications, policy advice, or operational decisions.
    - The API may change without notice and without backward compatibility.

Markets.jl solves a **spatial partial-equilibrium model** of an economic sector
made of an arbitrary set of products, regions and transformation processes.
You describe the economy with data — demand and supply curves, processes that
turn some products into others, transport costs between regions — and a single
optimisation returns **production, consumption, trade and prices** for every
region and product.

The engine knows no specific product or sector. The same data types describe
any sector with primary supply, transformation and trade: an oil refinery
splitting crude into petrol, diesel and fuel oil, a steel chain from iron ore
to cars, a fishery under quota, a dairy turning milk into cheese. Two worked
examples ship with the package: the dairy, as the smallest model that is still
a market, and the forest-products sector (roundwood → sawnwood, panels, pulp,
paper and pellets, with sawmill residues recycled and smallwood either burnt or
milled) at full size.

GitHub: [https://github.com/sylvaticus/Markets.jl](https://github.com/sylvaticus/Markets.jl)

## Features

* **Market equilibrium as an optimisation.** The competitive equilibrium is
  computed as the maximum of net social surplus. Quantities are the primal
  solution; prices are the duals of the material-balance constraints.
* **Multi-stage transformation with joint products.** A process can have several
  outputs (e.g. a main product and a residue), and residues can be inputs to
  other processes.
* **Smooth input substitution.** Inputs grouped in a CES nest substitute for each
  other continuously as relative prices change; fixed-proportion (Leontief)
  inputs are also available.
* **Linked regional markets.** Regions trade at a transport cost, so equilibrium
  prices differ across regions by at most that cost.
* **Imperfect substitution between origins, where you want it.** Give a product
  an Armington elasticity and its regional varieties become imperfect
  substitutes: regions cross-haul and every region's price responds to supply
  and demand everywhere. An infinite elasticity, the default, is the
  homogeneous case above.
* **Bounds where they belong.** Plant capacities, port and pipeline ceilings,
  quotas and licences truncate the model where a real limit exists, and a binding
  one reports its rent.
* **Substitution in final demand.** Products that compete for the same use are
  grouped in a nest, so one product's price moves another's demand — and
  substitutes outside the model enter as curve shifters.
* **Taxes and subsidies.** Wedges on production, consumption, transformation or
  trade, per unit or ad valorem, with the revenue of each reported.
* **Tidy results.** Every result is a `DataFrame`.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/sylvaticus/Markets.jl")   # this model
```

`Pkg.add("Markets")` would install the registered v0.0.1 instead, which is the
older and entirely different formulation described in
[Previous formulation (v0.0.1)](@ref); its documentation is kept under the
`v0.0.2` entry of the version selector.

## Quick start

```julia
using Markets
include(joinpath(pkgdir(Markets), "examples", "minimal", "minimal_market.jl"))
result.prices
```

That script is the [Two regions and a dairy](@ref) page, generated from the
source with Literate.jl, as is the larger
[forest example](@ref "The forest sector: France in the world").

## How this documentation is organised

* [Using the module](@ref) — how to describe an economy, solve it and read the
  results.
* **Examples** — [Two regions and a dairy](@ref) first, the smallest model
  that is still a market; then
  [The forest sector: France in the world](@ref), the same machinery at full
  size, from the data to a storm-salvage scenario.
* [Modelling choices](@ref) — what the model assumes and why: partial
  equilibrium, surplus maximisation, spatial trade, CES substitution, dynamics.
* [Code implementation](@ref) — the exact mathematical program the code builds,
  how prices are recovered, numerical details, and the [API reference](@ref).
* [References](@ref) — two lists: the papers behind each concept, cited as
  *Author (year)*, and open-source code in which the same concept can be read,
  cited as *[C1]*, *[C2]* and so on.
