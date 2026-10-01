# Using the module

Using Markets.jl takes three steps:

1. **Describe the economy** with the data types [`DemandSpec`](@ref),
   [`SupplySpec`](@ref), [`Process`](@ref) (built from [`leontief`](@ref) and
   [`ces`](@ref) nests), collected in a [`MarketData`](@ref).
2. **Solve** it with [`solve_market`](@ref).
3. **Read the results** from the `DataFrame`s in the returned [`Results`](@ref).

The engine never hard-codes a product, process or region name: everything
specific to a sector lives in the data.

## Describing an economy

### Regions and products

Regions and products are plain `Symbol`s:

```julia
regions  = [:EU, :NA, :AS]
products = [:swr, :hwr, :chips, :pulp, :sawn_sw, :sawn_hw, :panel, :paper]
```

Every product gets a material balance in every region, whether or not it is
produced, consumed or traded there.

### Final demand

A [`DemandSpec`](@ref) gives a constant-elasticity demand curve for one
product in one region, calibrated to pass through a reference point: at price
`p0` the quantity demanded is `q0`. `elasticity` is the (positive) own-price
elasticity ``\eta`` and must be greater than 1 (see
[Modelling choices](@ref "Demand and supply curves")).

```julia
DemandSpec(:sawn_sw, :EU; p0 = 250, q0 = 90, elasticity = 1.3)
```

Only products with a `DemandSpec` in a region are consumed there; other
products are intermediates in that region.

### Primary supply

A [`SupplySpec`](@ref) gives a constant-elasticity supply curve for a product
that enters the economy from outside the modelled processes — roundwood from
the forest, ore from a mine, crops from land:

```julia
SupplySpec(:swr, :NA; p0 = 65, q0 = 380, elasticity = 0.6)
```

### Processes and input nests

A [`Process`](@ref) turns inputs into outputs at a level of activity ``z``
chosen by the model. It has:

* a list of **input nests**. Each nest is one input *requirement* per unit of
  activity. Different nests are needed in fixed proportion to each other;
  products inside the same nest substitute for each other.
* a list of **outputs**, as `product => yield` pairs per unit of activity. Several
  outputs make a joint-production process (a main product and its residues).
* a **value-added cost** `vacost` per unit of activity, covering everything not
  modelled as an explicit input (labour, energy, capital, ...).
* optionally the `regions` where it may operate (default: all).

There are two kinds of nest:

* [`leontief(product, coeff)`](@ref leontief): exactly `coeff` units of a single
  product per unit of activity.
* [`ces(composite, products, shares; sigma)`](@ref ces): `composite` units of a
  CES bundle of several products. `shares` are the value shares of each input
  when all input prices are equal, and `sigma` (> 1) is the elasticity of
  substitution: the higher it is, the more the mix moves away from the
  expensive input when relative prices change.

For example, a sawmill that turns one m³ of softwood logs into 0.5 m³ of
lumber and 0.35 m³ of chips, and a panel mill that uses 1.3 m³ of a
substitutable mix of chips and logs:

```julia
Process(:sawmill_sw,
        [leontief(:swr, 1.0)],
        [:sawn_sw => 0.50, :chips => 0.35];
        vacost = 40)

Process(:panelmill,
        [ces(1.30, [:chips, :swr, :hwr], [0.55, 0.30, 0.15]; sigma = 2.5)],
        [:panel => 1.0];
        vacost = 120)
```

A process with two nests, e.g. `[ces(1.0, [:a, :b], [0.5, 0.5]; sigma = 2), leontief(:c, 0.2)]`,
needs one unit of the ``a``/``b`` bundle *and* 0.2 units of ``c`` per unit of
activity.

### Trade

`tradable` lists the products that may be traded, and `transport` maps each
allowed route `(product, from, to)` to a unit transport cost, per unit of
**quantity** shipped. A route missing from `transport` cannot be used.

```julia
transport = Dict((:swr, :NA, :AS) => 13.2, (:swr, :AS, :NA) => 13.2, ...)
```

### Putting it together

```julia
market = MarketData(; regions, products, demand, supply,
                      processes, tradable, transport)
```

## Solving

```julia
res = solve_market(market)                       # Ipopt, silent
res = solve_market(market; silent = false)       # show the Ipopt log
res = solve_market(market; optimizer = MyOpt)    # any JuMP nonlinear solver
```

`solve_market` warns if the solver does not report an optimal (or locally
optimal) solution.

## Reading the results

A [`Results`](@ref) holds six `DataFrame`s:

| Field | Columns | Content |
|:------|:--------|:--------|
| `production`  | region, product, quantity | primary supply + manufactured output |
| `consumption` | region, product, quantity | final demand |
| `trade`       | product, from, to, quantity | bilateral flows (positive only) |
| `net_trade`   | region, product, exports, imports, net | net = exports − imports |
| `prices`      | region, product, price | equilibrium price of every product in every region |
| `activity`    | region, process, level | process activity levels |

[`for_region`](@ref) slices all tables for one region at once. The input data
and the underlying JuMP model are kept in `res.data` and `res.model`.

## Worked example: the forest-products sector

The package ships an example economy in `examples/forest/example_data.jl`,
with three regions (Europe, North America, Asia-Pacific) and this product
chain:

```
  primary supply:   swr  softwood roundwood       hwr  hardwood roundwood

  sawmill_sw :  swr               → sawn_sw + chips    (fixed proportions)
  sawmill_hw :  hwr               → sawn_hw + chips    (fixed proportions)
  panelmill  :  CES(chips,swr,hwr) → panel              (smooth substitution)
  pulpmill   :  CES(swr,hwr,chips) → pulp               (chips compete here)
  papermill  :  pulp               → paper

  final demand:     sawn_sw, sawn_hw, panel, paper
```

Chips are a by-product of sawmilling and an input to both panel and pulp mills,
so the competition for residues is part of the equilibrium. Quantities are in
million m³ (wood) or million t (pulp, paper), prices in USD per unit; the
numbers are illustrative.

```@example forest
using Markets, DataFrames
include(joinpath(pkgdir(Markets), "examples", "forest", "example_data.jl"))
res = solve_market(example_market)
nothing # hide
```

Production, by region:

```@example forest
unstack(res.production, :region, :product, :quantity)
```

Consumption:

```@example forest
unstack(res.consumption, :region, :product, :quantity)
```

Prices — they differ across regions by at most the transport cost:

```@example forest
unstack(res.prices, :region, :product, :price)
```

Trade flows:

```@example forest
res.trade
```

Everything for one region:

```@example forest
for_region(res, :AS).net_trade
```

To run the full example with printed tables from a terminal, from the package
folder:

```bash
julia --project=. examples/forest/run_example.jl
```

### Running a scenario

Results are plain data, so scenarios are written by changing the data and
solving again. For example, a 20% increase in North American softwood
roundwood supply at every price:

```@example forest
sup = [s.product == :swr && s.region == :NA ?
           SupplySpec(s.product, s.region; p0 = s.p0, q0 = 1.2 * s.q0, elasticity = s.elasticity) : s
       for s in example_market.supply]
d2   = MarketData(example_market.regions, example_market.products, example_market.demand,
                  sup, example_market.processes, example_market.tradable, example_market.transport)
res2 = solve_market(d2)
comp = innerjoin(res.prices, res2.prices, on = [:region, :product], renamecols = "_base" => "_scen")
comp.change_pct = 100 .* (comp.price_scen ./ comp.price_base .- 1)
comp[comp.product .== :swr, :]
```

## How to amend an economy

* **Add a product**: add its symbol to `products`; give it a `DemandSpec` and/or
  `SupplySpec` where it is consumed or harvested; wire it into processes; list it
  in `tradable` and add `transport` entries if it can be shipped.
* **Add or change a process**: add a `Process`. Chains of any length are built by
  making the output of one process an input of another (e.g.
  roundwood → sawnwood → pallets). A product can be both consumed and used as an
  input.
* **Add a region**: add it to `regions` and give it demand, supply and transport
  data.
* **Change substitutability**: change the `sigma` of a CES nest. Values near 1
  give a near Cobb–Douglas mix; large values make the inputs close to perfect
  substitutes.

No change to the package code is needed for any of these.
