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
`p0` the quantity demanded is `q0`. `elasticity` is the own-price elasticity
``\eta``, which need only be positive — **inelastic demand, ``\eta < 1``, is
allowed and is what most forest products call for**:

```julia
DemandSpec(product = :sawn_sw, region = :EU, p0 = 250, q0 = 90, elasticity = 0.45)
```

With ``\eta \le 1`` the objective value is no longer an absolute measure of
welfare, only a basis for comparing scenarios; quantities and prices are
unaffected. See [Demand and supply curves](@ref).

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
  They are **not** the observed shares: calibrate them from a base-year trade
  matrix with [`armington_shares`](@ref), which inverts the CES demand
  condition for you.

  ```julia
  shares = armington_shares(flows = base_year_quantities,   # (origin, dest) => quantity
                            prices = delivered_prices,      # (origin, dest) => price
                            sigma = 4)
  ```

  Omitted, every origin available to a destination gets an equal share. An
  origin with no share (or a share of 0) is left out of that destination's
  composite, and so can never gain trade later.
* `sigma = Inf`, and any product not listed, is the homogeneous case above.
  Inside an [`OriginNest`](@ref) it instead pools that group of origins into a
  single good, which is how a set of regions that form one internal market —
  the regions of a country, say — enters a world model as one variety.

With an Armington product, local producers and local users no longer face the
same price, and both are reported — see [Reading the results](@ref).

!!! tip "Shares or elasticity?"
    Use the **shares** to say how much of an origin a market takes when prices
    are equal — that is where a home bias, or a regulatory penalty against an
    origin, belongs. Use **`sigma`** to say how readily buyers switch when
    prices move. A low `sigma` is not a barrier: it means buyers *cannot* get
    away from that origin, so an expensive one keeps its share rather than
    losing it.

#### Origins that are closer substitutes than others

A single `sigma` makes every origin substitute equally well for every other.
When some origins are near-interchangeable — shared grades, standards or
certification — and others are not, group them in an [`OriginNest`](@ref) with
its own, higher elasticity:

```julia
armington = [Armington(product = :sawn_sw, sigma = 2.5,
                       nests = [OriginNest(sigma = 12, origins = [:EU, :NA])],
                       shares = shares)]
```

Buyers then swap EU for NA sawnwood readily (σ = 12) and replace either with
Asian sawnwood only slowly (σ = 2.5). Origins left out of every group are
direct members of the composite, groups may contain groups, and a group must be
at least as substitutable inside as it is with the outside — the engine rejects
the reverse, which would be inconsistent.

#### Different structures in different markets

Standards and non-tariff measures are set by the importer, so they are
asymmetric. Give a specification a `destination` to apply it to those markets
only, and keep an `:all` entry as the fallback:

```julia
armington = [
    # the EU market keeps Asian sawnwood at arm's length, EU and NA are alike
    Armington(product = :sawn_sw, destination = :EU, sigma = 2.5,
              nests = [OriginNest(sigma = 12, origins = [:EU, :NA])], shares = shares),
    # every other market substitutes freely between all three origins
    Armington(product = :sawn_sw, sigma = 8, shares = shares)]
```

The elasticity, the groups and the shares may all differ by destination. What
may not differ is whether the product is Armington at all: it has one variety
per origin everywhere, or none anywhere.

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

Either price column is `missing` where it has no meaning: a region that cannot
produce an Armington product has no `producer_price` for it, and a region that
uses none of one has no `price`. Prices are signed, so a by-product that must
be disposed of at a cost shows a negative one.

`price` and `producer_price` differ only for an [`Armington`](@ref) product,
where local users buy a composite of all origins' varieties while local
producers sell their own; for every other product the two columns are equal.

[`for_region`](@ref) slices all tables for one region at once. The input data
and the underlying JuMP model are kept in `res.data` and `res.model`.

## A worked example

The [Forest example](@ref "The forest sector: France in the world") page builds
a complete economy with this API and solves it: eight regions (four of them
French), the roundwood-to-paper chain, joint products, a CES input bundle, a
process restricted to the regions that have the capacity, and a three-level
Armington structure in which the French regions are perfect substitutes for
each other, close substitutes for the rest of the EU and more distant ones for
the other blocs. It then runs a storm-salvage scenario and traces how far the
shock travels.

The page is generated from
[`examples/forest/forest_market.jl`](https://github.com/sylvaticus/Markets.jl/blob/main/examples/forest/forest_market.jl),
which is an ordinary script: run it to get the same tables in a terminal.

```bash
julia --project=. examples/forest/forest_market.jl
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
