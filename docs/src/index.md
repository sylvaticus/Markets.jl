# Markets.jl

!!! danger "NOT RELEASED — UNDER ACTIVE DEVELOPMENT AND VALIDATION"
    **This model is not finished, not calibrated and not validated. Do not use
    its results for anything that matters yet.**

    * **It is not released.** The current version, 0.1.0, is not registered,
      and the API changes between commits without notice.
    * **`Pkg.add("Markets")` installs a different model.** The registry still
      has v0.0.1, the 2022 equation-system formulation that this version
      replaced entirely — see
      [Previous formulation (v0.0.1)](@ref). Install from the repository
      instead.
    * **Nothing here is calibrated.** Every elasticity, share, cost and
      reference quantity in the forest example is illustrative — plausible in
      order of magnitude, fitted to nothing. No part of the model has been
      confronted with observed data, and the shares of the example's Armington
      nests come from a gravity rule rather than from a trade matrix.
    * **What *is* checked** is the internal logic. The test suite verifies the
      equilibrium conditions themselves: prices equal to inverse demand and to
      marginal cost, zero profit on active processes, no-arbitrage across trade
      routes, the CES demand and price-index conditions, and the limits in
      which the Armington formulation collapses to the homogeneous one. That
      makes the model internally coherent. It does not make it right about the
      world.
    * **Read some results more cautiously than others.** Unless a process is
      given a capacity, it has constant returns and nothing holds it in place,
      so processing location is bang-bang; see
      [What the model does not constrain](@ref). The forest example bounds its
      pellet mills but nothing else, so its prices and final demand are on
      firmer ground than its geography of milling.

Markets.jl solves a **spatial partial-equilibrium model** of an economic sector
made of an arbitrary set of products, regions and transformation processes.
You describe the economy with data — demand and supply curves, processes that
turn some products into others, transport costs between regions — and a single
optimisation returns **production, consumption, trade and prices** for every
region and product.

The engine knows no specific product or sector. The forest-products sector
(roundwood → sawnwood, panels, pulp, paper and pellets, with sawmill residues
recycled and smallwood either burnt or milled) is
the example shipped with the package, but the same data types describe any
sector with primary supply, multi-stage transformation and trade.

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
* **Bounds where they belong.** Mill capacities, port and corridor ceilings and
  allowable cuts truncate the model where a real limit exists, and a binding
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
include(joinpath(pkgdir(Markets), "examples", "forest", "forest_market.jl"))
res.prices
```

That script is the [Forest example](@ref "The forest sector: France in the world")
page, generated from the source with Literate.jl.

## How this documentation is organised

* [Using the module](@ref) — how to describe an economy, solve it and read the
  results.
* [The forest sector: France in the world](@ref) — a complete model worked
  through, from the data to a storm-salvage scenario.
* [Modelling choices](@ref) — what the model assumes and why: partial
  equilibrium, surplus maximisation, spatial trade, CES substitution, dynamics.
* [Code implementation](@ref) — the exact mathematical program the code builds,
  how prices are recovered, numerical details, and the [API reference](@ref).
