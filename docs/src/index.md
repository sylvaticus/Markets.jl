# Markets.jl

Markets.jl solves a **spatial partial-equilibrium model** of an economic sector
made of an arbitrary set of products, regions and transformation processes.
You describe the economy with data — demand and supply curves, processes that
turn some products into others, transport costs between regions — and a single
optimisation returns **production, consumption, trade and prices** for every
region and product.

The engine knows no specific product or sector. The forest-products sector
(roundwood → sawnwood, panels, pulp, paper, with sawmill residues recycled) is
the example shipped with the package, but the same data types describe any
sector with primary supply, multi-stage transformation and trade.

GitHub: [https://github.com/sylvaticus/Markets.jl](https://github.com/sylvaticus/Markets.jl)

!!! warning
    Pre-alpha status: the API may change between versions.

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
* **Tidy results.** Every result is a `DataFrame`.

## Installation

```julia
using Pkg
Pkg.add("Markets")                                          # registered version
Pkg.add(url = "https://github.com/sylvaticus/Markets.jl")   # development version
```

## Quick start

```julia
using Markets
include(joinpath(pkgdir(Markets), "examples", "forest", "example_data.jl"))
res = solve_market(example_market)
res.prices
```

## How this documentation is organised

* [Using the module](@ref) — how to describe an economy, solve it and read the
  results, with the forest example worked through.
* [Modelling choices](@ref) — what the model assumes and why: partial
  equilibrium, surplus maximisation, spatial trade, CES substitution, dynamics.
* [Code implementation](@ref) — the exact mathematical program the code builds,
  how prices are recovered, numerical details, and the [API reference](@ref).
