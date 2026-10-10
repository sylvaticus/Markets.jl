# Markets.jl

Solve a spatial partial-equilibrium model of an economic sector made of an
arbitrary set of products, regions and transformation processes.

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sylvaticus.github.io/Markets.jl/dev/)
[![Build status (Github Actions)](https://github.com/sylvaticus/Markets.jl/workflows/CI/badge.svg)](https://github.com/sylvaticus/Markets.jl/actions)
[![codecov.io](http://codecov.io/github/sylvaticus/Markets.jl/coverage.svg?branch=main)](http://codecov.io/github/sylvaticus/Markets.jl?branch=main)

> [!WARNING]
> ### NOT RELEASED — UNDER ACTIVE DEVELOPMENT AND VALIDATION
>
> **This model is not finished, not calibrated and not validated. Do not use
> its results for anything that matters yet.**
>
> * **It is not released.** The current version, 0.1.0, is not registered, and
>   the API changes between commits without notice.
> * **`Pkg.add("Markets")` installs a different model.** The registry still has
>   v0.0.1, the 2022 equation-system formulation that this version replaced
>   entirely. Install from this repository instead (see below).
> * **Nothing here is calibrated.** Every elasticity, share, cost and reference
>   quantity in the forest example is illustrative — plausible in order of
>   magnitude, fitted to nothing. No part of the model has been confronted with
>   observed data.
> * **What *is* checked** is the internal logic: the test suite verifies the
>   equilibrium conditions themselves — prices equal to inverse demand and
>   marginal cost, zero profit on active processes, no-arbitrage across trade
>   routes, the CES demand and price-index conditions, and the limits in which
>   the Armington model collapses to the homogeneous one.
> * **Read processing location with particular caution.** Processes have
>   constant returns and no capacity bounds, so where milling happens is
>   bang-bang: a small cost advantage takes the whole industry.

You describe the economy with data — demand and supply curves, processes that
transform products (with substitutable inputs and joint by-products), transport
costs between regions — and a single optimisation returns **production,
consumption, trade and prices** for every region and product. The equilibrium
is computed as the maximum of net social surplus (JuMP + Ipopt); prices are the
duals of the material balances.

Products are homogeneous by default (Samuelson spatial price equilibrium), or
imperfect substitutes by origin where you give them an **Armington** elasticity,
which lets regions cross-haul and links every region's price to supply and
demand everywhere.

The engine is sector-agnostic. The forest-products sector is the example
shipped in [`examples/forest/`](examples/forest/).

## Quick start

```julia
using Pkg; Pkg.add(url = "https://github.com/sylvaticus/Markets.jl")
using Markets
include(joinpath(pkgdir(Markets), "examples", "forest", "forest_market.jl"))
res.production; res.consumption; res.trade; res.prices
```

That example — a forest sector with four French regions and four world blocs —
is also the [Forest example](https://sylvaticus.github.io/Markets.jl/dev/generated/forest_market/)
page of the documentation, generated from the script itself with Literate.jl.

## Documentation

The [documentation](https://sylvaticus.github.io/Markets.jl/dev/) has four
parts: *Using the module*, the worked *Forest example*, *Modelling choices*,
and *Code implementation* (with the API reference).
