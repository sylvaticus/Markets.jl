# Markets.jl

Solve a spatial partial-equilibrium model of an economic sector made of an
arbitrary set of products, regions and transformation processes, with partial substitution in transformation, trade and demand.

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
>   quantity in the minimal and forest examples is illustrative — plausible in order of
>   magnitude, fitted to nothing. No part of the model has been confronted with
>   observed data.
> * **What *is* checked** is the internal logic: the test suite verifies the
>   equilibrium conditions themselves — prices equal to inverse demand and
>   marginal cost, zero profit on active processes, no-arbitrage across trade
>   routes, the CES demand and price-index conditions, and the limits in which
>   the Armington model collapses to the homogeneous one.
> * **Read processing location with particular caution.** Unless a process is
>   given a capacity, it has constant returns and nothing holds it in place, so
>   where milling happens is bang-bang: a small cost advantage takes the whole
>   industry. The forest example does bound its pellet mills, but nothing
>   else.

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

The smallest model that is still a market: two regions, two products, one
transformation, one trade route. The **Pasture** has the cows, the **City** has
the people, and milk is either drunk or turned into cheese.

```julia
using Pkg; Pkg.add(url = "https://github.com/sylvaticus/Markets.jl")
using Markets, DataFrames

economy = MarketData(
    regions   = [:Pasture, :City],
    products  = [:milk, :cheese],
    supply    = [SupplySpec(product = :milk, region = :Pasture, p0 = 400, q0 = 2000, elasticity = 0.6),
                 SupplySpec(product = :milk, region = :City,    p0 = 450, q0 =  300, elasticity = 0.6)],
    # milk is both drunk and made into cheese: each use has its own demand
    demand    = [DemandSpec(product = :milk,   region = :Pasture, p0 =  400, q0 = 200, elasticity = 0.3),
                 DemandSpec(product = :milk,   region = :City,    p0 =  450, q0 = 600, elasticity = 0.3),
                 DemandSpec(product = :cheese, region = :Pasture, p0 = 5300, q0 =  30, elasticity = 0.5),
                 DemandSpec(product = :cheese, region = :City,    p0 = 5350, q0 = 120, elasticity = 0.5)],
    processes = [Process(name    = :dairy,      # 10 t of milk and EUR 1300 a tonne
                         inputs  = [leontief(product = :milk, coeff = 10.0)], # ces nests are also available
                         outputs = [:cheese => 1.0],  # multi-outputs is also supported
                         vacost  = 1300)],
    transport_costs = Dict(
        (:milk,   :Pasture, :City) => 45.0, (:milk,   :City, :Pasture) => 45.0,
        (:cheese, :Pasture, :City) => 45.0, (:cheese, :City, :Pasture) => 45.0))

res = solve_market(economy)
unstack(res.prices, :region, :product, :price)
unstack(res.production, :region, :product, :quantity)
```

```
Price                                      Quantity
 Row │ region   cheese    milk             Row │ region   milk      cheese
─────┼─────────────────────────────        ────┼────────────────────────────────
   1 │ City      5353.55   445.855           1 │ Pasture  2002.56       149.936
   2 │ Pasture   5308.55   400.855           2 │ City      298.339  missing
```

Nothing above said how much milk to drink and how much to turn into cheese,
where the dairy should stand, or what should cross the border. All three fall
out. The 2301 tonnes of milk split themselves 801 drunk against 1500 to the
dairy, because the two uses bid for the same herds. **The City makes no cheese
at all**: a tonne of cheese costs 45 to carry, the ten tonnes of milk behind it
450. And milk does cross, but only the milk that is drunk. The cheese price
covers exactly its milk plus the making (10 × 400.855 + 1300 = 5308.55), and
each price gap between the regions is exactly the freight (5353.55 − 5308.55 =
45).

Walked through step by step in [Two regions and a dairy](https://sylvaticus.github.io/Markets.jl/dev/generated/minimal_market/);
the same machinery at full size — eight regions, eleven products, a storm — is in [Forest sector economics](https://sylvaticus.github.io/Markets.jl/dev/generated/forest_market/).
Both pages are generated from runnable scripts with Literate.jl, so every number
on them is the solver's own output.

## Documentation

The [documentation](https://sylvaticus.github.io/Markets.jl/dev/) has five
parts: *Using the module*, two worked *Examples*, *Modelling choices*,
*Code implementation* (with the API reference), and *References*.
