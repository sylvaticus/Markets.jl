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

Regions and products are named with strings:

```julia
regions  = ["EU", "US", "CN"]
products = ["crude", "gasoline", "diesel", "fuel_oil", "scrap", "pig_iron", "steel"]
```

A name may be written as a `Symbol` instead — `:EU` is accepted wherever `"EU"`
is — but it is converted on the way in, so what is stored, compared and
reported back is always the string. The two spellings may be mixed within one
economy, and a field read off a specification compares to the string that built
it:

```julia
SupplySpec(product = :crude, region = :EU, p0 = 550, q0 = 60, elasticity = 0.4).region == "EU"
```

Every product gets a material balance in every region, whether or not it is
produced, consumed or traded there.

### Final demand

A [`DemandSpec`](@ref) gives a constant-elasticity demand curve for one
product in one region, calibrated to pass through a reference point: at price
`p0` the quantity demanded is `q0`. `elasticity` is the own-price elasticity
``\eta``, which need only be positive — **inelastic demand, ``\eta < 1``, is
allowed and is what most fuels, foods and basic materials call for**:

```julia
DemandSpec(product = "gasoline", region = "EU", p0 = 900, q0 = 80, elasticity = 0.3)
```

With ``\eta \le 1`` the objective value is no longer an absolute measure of
welfare, only a basis for comparing scenarios; quantities and prices are
unaffected. See [Demand and supply curves](@ref).

Only products with a `DemandSpec` in a region are consumed there; other
products are intermediates in that region.

#### Products that substitute for each other

Give `product` a [`DemandNest`](@ref) instead of a single name and the curve
describes a **group of substitutes**: the `p0`, `q0` and `elasticity` then
belong to the composite, while the nest says how its members trade off against
each other.

```julia
DemandSpec(region = "EU", p0 = 30000, q0 = 12, elasticity = 0.5,
           product = DemandNest(sigma = 1.6, products = ["electric_car"], shares = [0.25],
                       nests = [DemandNest(sigma = 5.0, share = 0.75,
                                           products = ["petrol_car", "diesel_car"],
                                           shares   = [0.65, 0.35])]))
```

Petrol and diesel cars substitute readily here, electric cars less so with
either.

The data seldom come in that form: they give a price and a quantity for each
product, not a share for each and a reference point for the bundle.
[`ces_calibration`](@ref) converts the one into the other, so that at the
base-year prices the model buys exactly the base-year quantities:

```julia
fuel = ces_calibration(quantities = [7.0, 1.4], prices = [55.0, 250.0], sigma = 3)
DemandSpec(region = "FR", p0 = fuel.price, q0 = fuel.quantity, elasticity = 0.35,
           product = DemandNest(sigma = 3, products = ["firewood", "pellets"],
                                shares = fuel.shares))
```

`fuel.shares` are the shares at equal prices, `fuel.price` is the CES price
index of the bundle and `fuel.quantity` the base-year spending divided by it.
A nested bundle is calibrated from the inside out: the price and quantity of an
inner group make it one member of the level above. The
[forest example](@ref "Calibrating a bundle") does this for wood in
construction. The structure is a tree for the same reason the Armington one is:
a free matrix of cross-price elasticities is not the gradient of any benefit
function, and the model would stop being the market it claims to be. See
[Substitution in final demand](@ref).

Substitutes **outside** the model — public transport against cars — need no
nest at all. An exogenous price is a constant, so put the cross-price
elasticity into `shift` through [`exogenous_shift`](@ref).

Everything that moves demand other than its own price — income, population, the
price of a substitute outside the model — goes in `shift`, which scales the
whole curve and leaves `p0` and `q0` holding the calibration:

```julia
DemandSpec(product = "electricity", region = "EU", p0 = 90, q0 = 2700, elasticity = 0.3,
           shift = exogenous_shift(income = (level = 112.0, reference = 100.0,
                                             elasticity = 0.8)))
```

### Primary supply

A [`SupplySpec`](@ref) gives a constant-elasticity supply curve for a product
that enters the economy from outside the modelled processes — crude oil from a
field, ore from a mine, fish from the sea, crops from land:

```julia
SupplySpec(product = "crude", region = "US", p0 = 550, q0 = 600, elasticity = 0.4)
```

Two further fields carry the conditions of production:

* `shift` scales the curve, exactly as for demand. The size of the resource, the
  roads and ports, the workers and machines available all belong here, and
  [`exogenous_shift`](@ref) turns driver levels and their elasticities into the
  multiplier.
* `capacity` is a hard ceiling that no price can pass — a fishing quota, a
  mining licence, an exhausted resource. Whereas a shifter scales the curve, this
  truncates it, and the price then rises above the marginal cost on the curve
  by the **scarcity rent** of the resource.

```julia
SupplySpec(product = "cod", region = "EU", p0 = 3000, q0 = 300, elasticity = 0.40,
           shift    = exogenous_shift(stock = (level = 118.0, reference = 100.0,
                                               elasticity = 0.6)),
           capacity = 350)                      # the fishing quota
```

Which of the two to use for which driver, and when to model the factors
explicitly instead, is set out in
[Everything other than the price](@ref).

### Processes and input nests

A [`Process`](@ref) turns inputs into outputs at a level of activity ``z``
chosen by the model. It has:

* a list of **input nests**. Each nest is one input *requirement* per unit of
  activity. Different nests are needed in fixed proportion to each other;
  products inside the same nest substitute for each other.
* a list of **outputs**, as `product => yield` pairs per unit of activity. Several
  outputs make a joint-production process (a main product and its residues).
* a **value-added cost** `vacost` per unit of activity, covering everything not
  modelled as an explicit input (labour, energy, capital, ...) — a single number
  for every region, or a `region => cost` dictionary covering each region the
  process operates in.
* optionally the `regions` where it may operate (default: all), and a
  `capacity` bounding its activity — a single number for every region, or a
  `region => bound` dictionary. Without one, a process with constant returns
  concentrates wherever production is cheapest; with one, the bound binds and
  earns a quasi-rent reported in `Results.activity`.

There are two kinds of nest:

* [`leontief`](@ref): exactly `coeff` units of a single `product` per unit of
  activity.
* [`ces`](@ref): `composite` units of a CES bundle of several `products`.
  `shares` are the value shares of each input when all input prices are equal
  (equal shares if omitted; [`ces_calibration`](@ref) derives them from
  observed input quantities and prices), and `sigma` (≥ 0) is the elasticity of
  substitution: the higher it is, the more the mix moves away from the
  expensive input when relative prices change. `sigma = 0` needs the inputs in
  fixed proportions, `sigma = 1` is Cobb–Douglas, and `sigma = Inf` makes them
  perfect substitutes.

For example, a refinery that splits each tonne of crude into 0.45 t of
petrol, 0.35 t of diesel and 0.15 t of fuel oil, and a steel mill that needs
1.1 t of a substitutable mix of scrap and pig iron per tonne of steel:

```julia
Process(name    = "refinery",
        inputs  = [leontief(product = "crude", coeff = 1.0)],
        outputs = ["gasoline" => 0.45, "diesel" => 0.35, "fuel_oil" => 0.15],
        vacost  = 40)

Process(name    = "steel_mill",
        inputs  = [ces(composite = 1.10,
                       products  = ["scrap", "pig_iron"],
                       shares    = [0.4, 0.6],
                       sigma     = 2.0)],
        outputs = ["steel" => 1.0],
        vacost  = 150)
```

A process with two nests, e.g.

```julia
inputs = [ces(composite = 1.0, products = ["a", "b"], shares = [0.5, 0.5], sigma = 2),
          leontief(product = "c", coeff = 0.2)]
```

needs one unit of the ``a``/``b`` bundle *and* 0.2 units of ``c`` per unit of
activity.

### Trade

`transport_costs` maps each allowed route `(product, from, to)` to a unit
transport cost, per unit of **quantity** shipped. It is also the whole
declaration of trade: a route missing from it cannot be used, and a product
with no route anywhere is simply not traded, so there is no separate list of
tradable products to keep in step with it. Routes are directional, which is
what lets a corridor cost more one way than the other, so a product that moves
both ways needs an entry each way.

```julia
transport_costs = Dict(("crude", "US", "EU") => 25.0, ("crude", "EU", "US") => 25.0, ...)
```

A free route is written `=> 0.0`, not omitted.

A [`TradeLimit`](@ref) caps what may actually move: the capacity of a port or a
corridor, or a quota. It applies to the **sum** of the flows it covers, so one
entry can cap a single product on a single route or everything crossing a
border, and `weights` converts the quantities into whatever unit the ceiling is
written in:

```julia
trade_limits = [TradeLimit(name = "steel_quota", products = "steel", from = "CN", to = "EU",
                           capacity = 8.0),
                TradeLimit(name = "container_port", from = "CN", capacity = 6.0,
                           weights = Dict("cars" => 2.0, "wine" => 0.04))]
```

A binding ceiling shows up in `Results.limits` with its congestion rent, and
the price gap across the route it throttles becomes transport **plus** that
rent.

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
armington = [Armington(product = "wine", sigma = 4),              # equal shares
             Armington(product = "cars", sigma = 6,
                       shares = Dict(("EU", "EU") => 0.7, ("US", "EU") => 0.2, ("CN", "EU") => 0.1,
                                     ("US", "US") => 0.8, ("EU", "US") => 0.1, ("CN", "US") => 0.1,
                                     ("CN", "CN") => 0.9, ("EU", "CN") => 0.05, ("US", "CN") => 0.05))]
```

* `sigma` is the elasticity of substitution between origins, any value `≥ 0`:
  the lower it is, the more buyers stick to their usual origin when prices
  move. Below 1 the varieties are poor substitutes, or complements; `sigma = 1`
  is Cobb–Douglas (each origin keeps a fixed share of the spending) and
  `sigma = 0` fixes the mix of origins whatever the prices.
* `shares` are keyed `(origin, destination)` and are the shares that would be
  observed if all delivered prices were equal, so they carry the home bias.
  They are **not** the observed shares: calibrate them from a base-year trade
  matrix with [`armington_shares`](@ref), which inverts the CES demand
  condition for you — it is [`ces_calibration`](@ref) applied to each
  destination.

  ```julia
  shares = armington_shares(flows = base_year_quantities,   # (origin, dest) => quantity
                            prices = delivered_prices,      # (origin, dest) => price
                            sigma = 4)
  ```

  If omitted, every origin available to a destination gets an equal share. An
  origin with no share (or a share of 0) is left out of that destination's
  composite, and so can never gain trade later.

  With origin groups (below), [`armington_calibration`](@ref) builds the whole
  specification from the same two dictionaries — shares through every level of
  the tree, and the `scale` that keeps each composite in the physical units
  of the flows at the base year, so that base-year data balance exactly.
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
armington = [Armington(product = "cars", sigma = 2.5,
                       nests = [OriginNest(sigma = 12, origins = ["EU", "US"])],
                       shares = shares)]
```

Buyers then swap EU for US cars readily (σ = 12) and replace either with
Chinese cars only slowly (σ = 2.5). Origins left out of every group are
direct members of the composite, groups may contain groups, and a group must be
at least as substitutable inside as it is with the outside — the engine rejects
the reverse, which would be inconsistent.

#### Different structures in different markets

Standards and non-tariff measures are set by the importer, so they are
asymmetric. Give a specification a `destination` to apply it to those markets
only, and keep an `"all"` entry as the fallback:

```julia
armington = [
    # the EU market keeps Chinese cars at arm's length, EU and US are alike
    Armington(product = "cars", destination = "EU", sigma = 2.5,
              nests = [OriginNest(sigma = 12, origins = ["EU", "US"])], shares = shares),
    # every other market substitutes freely between all three origins
    Armington(product = "cars", sigma = 8, shares = shares)]
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
                    transport_costs = transport_costs)
```

Only `regions` and `products` are required; `demand`, `supply`, `processes` and
`transport_costs` all default to empty. When your variables already carry the
field names, Julia's shorthand saves the repetition:

```julia
market = MarketData(; regions, products, demand, supply, processes, transport_costs)
```

### Taxes and subsidies

A [`Tax`](@ref) drives a wedge between what one side pays and what the other
receives, on production, on consumption, on a transformation or on a trade
route. `rate` is an amount per unit, or a fraction of the base when
`relative = true`; a subsidy is a negative rate.

```julia
taxes = [Tax(name = "oil_royalty",  on = "supply",  product = "crude", region = "US", rate = 4.0),
         Tax(name = "fuel_vat",     on = "demand",  product = "gasoline", rate = 0.2,
             relative = true),
         Tax(name = "green_steel",  on = "process", process = "steel_mill", region = "EU",
             rate = -15.0),
         Tax(name = "steel_tariff", on = "trade",   product = "steel", to = "US", rate = 0.25,
             relative = true)]
```

Selectors left at `"all"` match everything, and several taxes may fall on one
base. `Results.taxes` reports what each raises — which matters, because a taxed
objective counts the revenue as a cost rather than a transfer and so is no
longer welfare. See [Taxes, subsidies and other wedges](@ref).

## Looking at an economy

Typing the name of a `MarketData` in the REPL gives one line, such as
`MarketData(3 regions, 7 products, 4 processes)`. Printing it, with
`print(market)`, gives a summary: what is supplied, consumed, made and used of
each product, the processes in the order they chain into one another (stage 1
uses only primary products, stage 2 also the outputs of stage 1, and so on),
the trade routes and their costs, the Armington specifications, the trade
limits and the taxes and subsidies. For a steel chain that ends in cars:

```
Transformations, by stage
  1  blast_furnace  1.6 iron_ore + 0.8 coal ───────→ pig_iron + 0.3 slag  cost 60
  2  steel_mill     1.1 CES(scrap, pig_iron; σ 2) ─→ steel                cost 150
  2  dispose_slag   slag ──────────────────────────→ (disposed of)        cost 5
  3  car_plant      0.9 steel ─────────────────────→ cars                 cost 20000 · cap EU 14, US 13, CN 30
```

The mathematical program the economy becomes is returned, built but not
solved, by [`jump_model`](@ref):

```julia
m = jump_model(market)      # the JuMP model solve_market would solve
print(m)                    # objective, balances, bounds, in full
latex_formulation(m)        # the same, in LaTeX
```

`latex_formulation` is JuMP's, re-exported so that it needs no `using JuMP`.
In a notebook or on a documentation page it is typeset; in a terminal it
prints the LaTeX source.

!!! warning "Check the signs against `print`"
    JuMP 1.32 drops the minus sign of a negative coefficient written in
    scientific notation when it typesets a model
    ([JuMP.jl#4268](https://github.com/jump-dev/JuMP.jl/issues/4268)). The
    benefit term of a demand with ``\eta < 1`` has a negative coefficient,
    usually a large one, so in most models some objective terms are shown with
    the wrong sign. `print(m)` is not affected; read the signs from it until
    the bug is fixed.

Its variables and constraints are named after what they stand for —
`D[EU,steel]`, `S[US,crude]`, `z[EU,steel_mill]`, `balance[EU,steel]`,
`trade_limit[steel_quota]` — and it takes the same `solver`,
`solver_options`, `silent` and `qfloor` keywords as `solve_market`, so it can
be modified and solved with JuMP directly. An ad valorem tax on trade is the
one thing missing from it: `solve_market` finds its money equivalent by
solving repeatedly, and the model returned is the first of those solves, with
the tax at zero.

## Solving

```julia
res = solve_market(market)                       # Ipopt, silent
res = solve_market(market; silent = false)       # show the Ipopt log
res = solve_market(market; solver = MyOpt)       # any JuMP nonlinear solver
res = solve_market(market; solver_options = (max_iter = 500, tol = 1e-9))
```

`solver_options` takes a `NamedTuple`, a `Dict` or a vector of pairs. A name
is the solver's own name for the option, written as a symbol or a string, and
an MOI attribute such as `MOI.TimeLimitSec()` is passed as it is. The options
are set after the ones Markets chooses for Ipopt (see
[Numerical details](@ref)), so they override them.

Before building anything, `solve_market` calls [`validate`](@ref) on the
economy. That check is not a formality: the model has one balance per declared
region and product, so a name outside those lists would get variables and no
balance — a mistyped process input would become a *free* input and the solver
would report success. Every name, duplicate and selector is therefore checked
first, and a problem is an `ArgumentError` saying which entry is at fault.

If the solver then returns anything but an optimal or locally optimal solution,
`solve_market` **throws**. An infeasible or unbounded economy comes back with
tables that look exactly like those of a solved one, so returning them quietly
is the more dangerous default. Pass `on_failure = "warn"` or `"ignore"` to get
them anyway, and read `res.status` to see what happened.

## Reading the results

A [`Results`](@ref) carries the solver's verdict in `status` and its own
message in `raw_status`, the `solve_time` in seconds and the `iterations` it
took (summed over all `solves` when an ad valorem tariff is found by iteration,
and `missing` for a solver that does not count them), the value of the
maximised `objective`, the `welfare` that objective implies once tax revenue is
counted as the transfer it is, and nine `DataFrame`s:

| Field | Columns | Content |
|:------|:--------|:--------|
| `production`  | region, product, quantity | primary supply + manufactured output |
| `supply`      | region, product, quantity, capacity, rent | primary supply alone, with its ceiling and the scarcity rent a binding one earns |
| `consumption` | region, product, quantity | final demand |
| `trade`       | product, from, to, quantity | bilateral flows (positive only) |
| `net_trade`   | region, product, exports, imports, net | net = exports − imports |
| `prices`      | region, product, price, producer\_price | what local users pay, and what local producers get |
| `activity`    | region, process, level, capacity, rent | process activity, its bound, and the quasi-rent a binding bound earns |
| `limits`      | limit, capacity, used, rent | trade ceilings and their congestion rent |
| `taxes`       | tax, on, revenue | what each instrument raises |

Either price column is `missing` where it has no meaning: a region that cannot
produce an Armington product has no `producer_price` for it, and a region that
uses none of one has no `price`. Prices are signed, so a by-product that must
be disposed of at a cost shows a negative one.

`price` and `producer_price` differ only for an [`Armington`](@ref) product,
where local users buy a composite of all origins' varieties while local
producers sell their own; for every other product the two columns are equal.

[`for_region`](@ref) slices all tables for one region at once. The input data
and the underlying JuMP model are kept in `res.data` and `res.model`.

## Worked examples

[Two regions and a dairy](@ref) is the smallest complete model — two products,
one process, one route, and milk that is both drunk and turned into cheese —
and the place to start.

The [forest example](@ref "The forest sector: France in the world") page builds
a complete economy with this API and solves it: eight regions (four of them
French), the roundwood-to-paper-and-pellets chain, a product that is both burnt
directly and milled, joint products, CES input bundles, a
process restricted to the regions that have the capacity, and a three-level
Armington structure in which the French regions are near-perfect substitutes
for each other, close substitutes for the rest of the EU and more distant ones for
the other blocs. It builds all of it from observed base-year data read from
CSV files — calibrating every parameter it does not observe or choose — and
compares the equilibrium with those data before running a storm-salvage
scenario that traces how far the shock travels, and a gas-price one.

The page is generated from
[`examples/forest/forest_market.jl`](https://github.com/sylvaticus/Markets.jl/blob/main/examples/forest/forest_market.jl),
which is an ordinary script: run it to get the same tables in a terminal.

```bash
julia --project=. examples/forest/forest_market.jl
```

## How to amend an economy

* **Add a product**: add its name to `products`; give it a `DemandSpec` and/or
  `SupplySpec` where it is consumed or produced; wire it into processes;
  add `transport_costs` entries if it can be shipped.
* **Add or change a process**: add a `Process`. Chains of any length are built by
  making the output of one process an input of another (e.g.
  iron ore → steel → cars). A product can be primary-supplied,
  manufactured, consumed as final demand and used as an input all at once — its
  balance simply sums every source against every use.
* **Add a region**: add it to `regions` and give it demand, supply and transport
  data.
* **Change substitutability between inputs**: change the `sigma` of a CES nest.
  Values near 1 give a near Cobb–Douglas mix; large values make the inputs close
  to perfect substitutes.
* **Change substitutability between origins**: add or amend an `Armington`
  entry. Leave a product out of `armington` (or give it `sigma = Inf`) for the
  homogeneous, spatial-equilibrium treatment.

No change to the package code is needed for any of these.
