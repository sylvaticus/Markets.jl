# Markets.jl

Solve a spatial partial-equilibrium model of an economic sector made of an
arbitrary set of products, regions and transformation processes.

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sylvaticus.github.io/Markets.jl/dev/)
[![Build status (Github Actions)](https://github.com/sylvaticus/Markets.jl/workflows/CI/badge.svg)](https://github.com/sylvaticus/Markets.jl/actions)
[![codecov.io](http://codecov.io/github/sylvaticus/Markets.jl/coverage.svg?branch=main)](http://codecov.io/github/sylvaticus/Markets.jl?branch=main)

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

Pre-alpha status: the API may change between versions.

## Quick start

```julia
using Pkg; Pkg.add("Markets")
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
