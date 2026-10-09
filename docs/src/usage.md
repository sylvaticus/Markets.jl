# Using the module

Using Markets.jl takes three steps:

1. **Describe the economy** with the data types [`DemandSpec`](@ref),
   [`SupplySpec`](@ref), [`Process`](@ref) (built from [`leontief`](@ref) and
   [`ces`](@ref) nests), collected in a [`MarketData`](@ref).
2. **Solve** it with [`solve_market`](@ref).
3. **Read the results** from the `DataFrame`s in the returned [`Results`](@ref).

The engine never hard-codes a product, process or region name: everything
specific to a sector lives in the data.

All the types and builders below take **keyword arguments only**, so every
number in a model description says what it is. Optional arguments have
defaults; the others must be given.

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
DemandSpec(product = :sawn_sw, region = :EU, p0 = 250, q0 = 90, elasticity = 1.3)
```

Only products with a `DemandSpec` in a region are consumed there; other
products are intermediates in that region.

### Primary supply

A [`SupplySpec`](@ref) gives a constant-elasticity supply curve for a product
that enters the economy from outside the modelled processes — roundwood from
the forest, ore from a mine, crops from land:

```julia
SupplySpec(product = :swr, region = :NA, p0 = 65, q0 = 380, elasticity = 0.6)
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

* [`leontief`](@ref): exactly `coeff` units of a single `product` per unit of
  activity.
* [`ces`](@ref): `composite` units of a CES bundle of several `products`.
  `shares` are the value shares of each input when all input prices are equal
  (equal shares if omitted), and `sigma` (> 1) is the elasticity of
  substitution: the higher it is, the more the mix moves away from the
  expensive input when relative prices change.

For example, a sawmill that turns one m³ of softwood logs into 0.5 m³ of
lumber and 0.35 m³ of chips, and a panel mill that uses 1.3 m³ of a
substitutable mix of chips and logs:

```julia
Process(name    = :sawmill_sw,
        inputs  = [leontief(product = :swr, coeff = 1.0)],
        outputs = [:sawn_sw => 0.50, :chips => 0.35],
        vacost  = 40)

Process(name    = :panelmill,
        inputs  = [ces(composite = 1.30,
                       products  = [:chips, :swr, :hwr],
                       shares    = [0.55, 0.30, 0.15],
                       sigma     = 2.5)],
        outputs = [:panel => 1.0],
        vacost  = 120)
```

A process with two nests, e.g.

```julia
inputs = [ces(composite = 1.0, products = [:a, :b], shares = [0.5, 0.5], sigma = 2),
          leontief(product = :c, coeff = 0.2)]
```

needs one unit of the ``a``/``b`` bundle *and* 0.2 units of ``c`` per unit of
activity.

### Trade

`tradable` lists the products that may be traded, and `transport` maps each
allowed route `(product, from, to)` to a unit transport cost, per unit of
**quantity** shipped. A route missing from `transport` cannot be used.

```julia
transport = Dict((:swr, :NA, :AS) => 13.2, (:swr, :AS, :NA) => 13.2, ...)
```

By default a traded product is **homogeneous**: the same good wherever it comes
from. Regional prices then differ by at most the transport cost and no region
imports and exports it at once.

### Imperfect substitution between origins (optional)

Listing a product in [`Armington`](@ref) makes the varieties of the different
origins imperfect substitutes instead: each region uses a CES composite of
them, so it can import and export the same product at once and its price
responds to supply and demand everywhere rather than only to the cheapest
source.

```julia
armington = [Armington(product = :paper, sigma = 4),             # equal shares
             Armington(product = :panel, sigma = 6,
                       shares = Dict((:EU, :EU) => 0.7, (:NA, :EU) => 0.2, (:AS, :EU) => 0.1,
                                     (:NA, :NA) => 0.8, (:EU, :NA) => 0.1, (:AS, :NA) => 0.1,
                                     (:AS, :AS) => 0.9, (:EU, :AS) => 0.05, (:NA, :AS) => 0.05))]
```

* `sigma` is the elasticity of substitution between origins. It must be `> 1`;
  the lower it is, the more buyers stick to their usual origin when prices
  move.
* `shares` are keyed `(origin, destination)` and are the shares that would be
  observed if all delivered prices were equal, so they carry the home bias.
  Calibrate them on a base-year trade matrix. Omitted, every origin available
  to a destination gets an equal share. An origin with no share (or a share of
  0) is left out of that destination's composite.
* `sigma = Inf`, and any product not listed, is the homogeneous case above.

With an Armington product, local producers and local users no longer face the
same price, and both are reported — see [Reading the results](@ref).

### Putting it together

```julia
market = MarketData(regions   = regions,
                    products  = products,
                    demand    = demand,
                    supply    = supply,
                    processes = processes,
                    tradable  = tradable,
                    transport = transport)
```

Only `regions` and `products` are required; `demand`, `supply`, `processes`,
`tradable` and `transport` default to empty. When your variables already carry
the field names, Julia's shorthand saves the repetition:

```julia
market = MarketData(; regions, products, demand, supply, processes, tradable, transport)
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
| `prices`      | region, product, price, producer\_price | what local users pay, and what local producers get |
| `activity`    | region, process, level | process activity levels |

`price` and `producer_price` differ only for an [`Armington`](@ref) product,
where local users buy a composite of all origins' varieties while local
producers sell their own; for every other product the two columns are equal.

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
solving again. A three-line helper that copies an economy with some fields
replaced makes this comfortable:

```@example forest
reconfigure(d; kwargs...) =
    MarketData(; regions = d.regions, products = d.products, demand = d.demand,
                 supply = d.supply, processes = d.processes, tradable = d.tradable,
                 transport = d.transport, armington = d.armington, kwargs...)
nothing # hide
```

For example, a 20% increase in North American softwood roundwood supply at
every price:

```@example forest
sup = [s.product == :swr && s.region == :NA ?
           SupplySpec(product = s.product, region = s.region,
                      p0 = s.p0, q0 = 1.2 * s.q0, elasticity = s.elasticity) : s
       for s in example_market.supply]
res2 = solve_market(reconfigure(example_market, supply = sup))
comp = innerjoin(res.prices, res2.prices, on = [:region, :product], renamecols = "_base" => "_scen")
comp.change_pct = 100 .* (comp.price_scen ./ comp.price_base .- 1)
comp[comp.product .== :swr, [:region, :product, :price_base, :price_scen, :change_pct]]
```

### Imperfect substitution between origins

In the solution above, paper is not traded at all: regional paper prices differ
by less than the cost of shipping it, so no shipment pays for itself.

```@example forest
res.trade[res.trade.product .== :paper, :]
```

That is the homogeneous-product logic. Making paper an [`Armington`](@ref)
product instead — buyers mildly prefer their usual origin, 85% of a region's
paper coming from home at equal delivered prices — gives the two-way trade that
is actually observed:

```@example forest
home   = 0.85
shares = Dict((o, r) => (o == r ? home : (1 - home) / (length(regions) - 1))
              for o in regions, r in regions)
res_a  = solve_market(reconfigure(example_market,
             armington = [Armington(product = :paper, sigma = 4, shares = shares)]))
res_a.trade[res_a.trade.product .== :paper, :]
```

Every region now both imports and exports paper, and local producers no longer
face the price local users pay:

```@example forest
res_a.prices[res_a.prices.product .== :paper, :]
```

Raising `sigma` tightens the varieties together again, and the solution walks
back to the homogeneous one — which is what `sigma = Inf`, the default, builds
exactly:

```@example forest
homog  = res.prices[res.prices.product .== :paper, :price]
ladder = DataFrame(sigma = Float64[], paper_trade = Float64[], max_price_gap_pct = Float64[])
for σ in [4, 20, 150, Inf]
    a = solve_market(reconfigure(example_market,
            armington = [Armington(product = :paper, sigma = σ, shares = shares)]))
    traded = sum(a.trade[a.trade.product .== :paper, :quantity])
    gap    = maximum(abs.(a.prices[a.prices.product .== :paper, :price] ./ homog .- 1))
    push!(ladder, (σ, traded, 100gap))
end
ladder
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
* **Change substitutability between inputs**: change the `sigma` of a CES nest.
  Values near 1 give a near Cobb–Douglas mix; large values make the inputs close
  to perfect substitutes.
* **Change substitutability between origins**: add or amend an `Armington`
  entry. Leave a product out of `armington` (or give it `sigma = Inf`) for the
  homogeneous, spatial-equilibrium treatment.

No change to the package code is needed for any of these.
