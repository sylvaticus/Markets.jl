# Accounts, calibration and scope

## What the model does not constrain

Two limitations are inherent to the formulation rather than to the data, and
both bear on how results should be read.

**Processing location does not respond smoothly.** Processes have constant
returns to scale and a constant value-added cost, so a process has no supply
response of its own: below a threshold set by freight it does not move at all,
and above it the whole industry moves at once.
[Why processes, and not supply and demand per transformation](@ref) measures
this and explains where it bites — chiefly on products left homogeneous, since
an [`Armington`](@ref) specification on the output restores a smooth response.

Otherwise the remedy is to say what capacity exists. A [`Limit`](@ref) on the
activity of a process bounds it, region by region if need be, and a binding
bound earns a quasi-rent reported in `Results.activity` and `Results.limits` —
what another unit of plant capacity would be worth there. Without such a bound, and without Armington
on its output, read the geography of *processing* far more cautiously than that
of prices or of final demand.

**Prices come from a continuous relaxation.** All quantities are continuous, so
the model knows nothing of indivisible investments, minimum efficient scale or
shutdown decisions.

Linear constraints on any of the model's quantities — a resource shared by
several products, a quota, a minimum share — can be added with
[`@limit`](@ref) without leaving the formulation: they keep the problem convex,
and each comes back with its rent; see [Limits on any quantity](@ref).

## Physical accounts

Every quantity the model solves for is physical — what each region produced,
what it shipped to every other, what its households consumed, what its plants
took in — and every balance between them is an equality. Two kinds of
quantity are not sums, though, but CES aggregates:

* the **composite** a region uses of an [`Armington`](@ref) product, made of
  the varieties of every origin;
* the **bundle** a process uses from a [`ces`](@ref) nest, made of several
  inputs.

A CES aggregate equals the sum of its members only at the mix it was
calibrated at. Away from it, a mixed bundle counts for more, or less, than
the physical quantity in it: a region's households and plants can draw a
composite of 100 out of 104 m³ delivered, or a mill can turn 2.47 m³ of wood
into a tonne of pellets where its recipe says 2.39. That is not an error of
the solver but what a CES is — a production or preference function, with a
gain from mixing built in — and it is the price of the smooth, scale-free
substitution it gives.

The results therefore keep two accounts side by side:

* in **composite units**, the quantities the demand curves, the recipes and
  the prices refer to. This account is exact in *value*: the composite times
  the price its users pay is what was spent on the varieties, delivered.
* in **physical units**, the quantities that moved. The varieties a region
  received are its `received`, in `Results.composites`, and every use of the
  region's composite — final demand and each plant — is credited with the
  same physical content per unit, the `ratio` received / composite. That
  assumes every user in a region gets the same mix of origins, the one
  assumption the account makes. Final consumption has a `physical` column,
  process inputs are in `Results.inputs`, and `Results.recipes` sets the
  physical inputs of each nest against what its recipe requires: the `gap`
  is the quantity the aggregation used up or made up.

With them the physical account closes: production equals shipments, shipments
equal receipts, receipts equal the physical uses, and a plant's physical
inputs equal its recipe plus the gap. In the
[forest examples](@ref "The forest sector: calibrating with partial data") the
gaps come to a few hundredths of a per cent of the world's wood, and up to a
few per cent for one product in one region under a large shock — small, but
worth reporting wherever the model feeds a physical balance, such as carbon or
resource use.

### What it would take to make the composites add up

Where a physical identity has to hold exactly, there are exact aggregators:
`sigma = 0` (fixed proportions) and `sigma = Inf` (perfect substitutes, a
plain sum), and, for inputs, a choice between alternative fixed recipes
written as separate processes. Keeping imperfect substitution *and* an
aggregate that is a plain sum is possible too, but each way of doing it gives
up something else, and there are three properties to choose from:

| | one convex program, prices as duals | substitution depends on price ratios only | the aggregate is a plain sum |
|:---|:---:|:---:|:---:|
| CES, as here | ✓ | ✓ | |
| logit shares, with a variety term in the objective | ✓ | | ✓ |
| additive CES ([van der Mensbrugghe & Peters, 2016](@ref vdmp2016)) | | ✓ | ✓ |

The logit keeps the program and the plain sum, but its shares respond to price
*differences*: doubling every price moves them, and users pay more per unit of
the aggregate than was spent per unit. The additive CES keeps the plain sum and
the ratios — its shares follow the usual CES rule, and the aggregate is priced
at the average unit value — but it chooses its mix as if minimising the cost
of a CES aggregate while counting the result as a sum, which no objective in
the quantities alone describes; it belongs in an equilibrium model posed as a
complementarity problem, the paradigm this package does not take. Its results
also depend on how base-year values are split into prices and volumes, which
those of the CES do not. Van der Mensbrugghe & Peters find the deviations of
the standard form, for land in a global model, of 3.5% in the world and over
8% in some countries — the same order as here.

### Material accounts

What the goods are made of — carbon, dry matter, a metal — can be traced
through the same solution with [`material_accounts`](@ref), following
[Lorang & Gerlagh (2026)](@ref lorang2026). A [`Material`](@ref) enters with
primary supply at a given content per unit; it travels with the goods, each
region's production carrying the content of what made it and its users the
average content of what they received; and in a process it is either lost —
burnt, emitted, discarded, at a given share — or carried into the outputs.
How much of it a unit of a manufactured good carries is therefore not given
but follows from the solution: a mill that changes its mix changes the
content of what it makes, and the gap of the physical account becomes, in
these terms, a change in the material intensity of the output rather than
a quantity left unexplained.

These accounts close exactly — every region's material balance, and the
world's: what was extracted was consumed or lost — and the `consumed` column
is each region's footprint, the material embodied in its final consumption.
They describe the solution and do not change it, as long as nothing in the
model depends on the content of the goods. A tax per tonne of embodied
material, or a cap on it, would make the intensities part of the equilibrium;
that is not implemented.



## Calibration

A reasonable first reaction to all of this is: if I already know what was
produced last year, and at what price, what is the model for?

The answer is that an observation is a **point**, and a point with an
elasticity attached is a **curve**. The model never tells you the base year —
you measured it. It tells you where the curves' intersection moves when
something changes, and it does so for every product, region and process at
once, with material balance, zero profit and no-arbitrage all holding
simultaneously. That joint consistency is the thing no observation can give
you, because it is a property of the system rather than of any one market.

### A consistent benchmark is a fixed point

There is a sharp way to see what calibration does. If the data you anchor on
*are* an equilibrium, the model reproduces them exactly and adds nothing:

```
A  steel  price 182.3283 → 182.3283    relative difference 2.5e-13
B  steel  price 190.3283 → 190.3283    relative difference 1.3e-13
A  ore    price  63.9552 →  63.9552    relative difference 5.1e-13
```

That is *benchmark replication*, and reproducing it is the standard test that
the data and the structure agree with each other.

Real data are never an equilibrium — they carry measurement error, aggregation
over products and periods, margins the model does not represent, and markets
that were not clearing. Anchor one price 5% away from its consistent value and
the model refuses to return it:

```
A  steel  anchored at 191.445   model returns  183.64   gap 4.1%
A  ore    anchored at  63.955   model returns   64.77   gap 1.3%
```

It splits the difference, because your numbers cannot all be true at once. That
gap is not a nuisance: it is the model telling you so, and its size is a
measure of how far your data are from being a coherent description of a market.

### What is actually calibrated

Entering `p0` and `q0` is not calibration, it is data entry. Calibration is
choosing the parameters you do **not** observe so that the benchmark becomes an
equilibrium:

* `vacost` of each process in each region, from the observed output prices
  less the observed cost of its inputs there — which makes zero profit hold
  exactly at base-year prices;
* `transport_costs` on each route, from the observed price gap between its ends on
  routes that were actually used;
* [`armington_calibration`](@ref) for the trade shares and the scale of each
  composite, from the base-year bilateral matrix at base-year delivered prices,
  through every level of the origin groups — not from the raw observed shares,
  for the reason given under
  [Imperfect substitution between origins (Armington)](@ref);
  [`armington_shares`](@ref) does the shares alone for a flat CES;
* [`ces_calibration`](@ref) for the shares of any other CES group — the
  substitutes in a [`DemandNest`](@ref), the inputs of a [`ces`](@ref) nest —
  and for the `p0` and `q0` of a demand curve over a bundle, which are its CES
  price index and the base-year spending divided by that index;
* `shift` on the curves, to absorb what remains.

The forest examples do all of it from base-year data read from files:
[from partial, rounded statistics](@ref "The forest sector: calibrating with partial data"),
where the model lands near the data rather than on them, and
[from complete, consistent data](@ref "The forest sector: complete data, and what calibration changes"),
where it returns them exactly.

### Calibrating with partial data

Data will be missing, and what is there will not quite agree. Neither stops a
calibration, and neither matters as much as it seems, because scenario results
are differences from the model's own reference run, never from the data. The
question is not whether the reference run reproduces every observed number,
but whether it is right about what the differences depend on.

That depends on the mechanics of the model:

* **Constant elasticities give the same percentage response wherever the
  reference run sits.** A curve shifted by 10% moves its price by the same
  percentage whether the price is 90 or 110, so much of a scenario is robust to
  a reference run that sits some way from the data.
* **Nothing protects a scenario from a reference run that is wrong about what
  the shock runs into.** A capacity or a logistic limit that the reference run
  places nearer to binding than it was binds sooner, and changes the result
  in kind, not only in size. Market shares matter the same way — how strongly
  one region's price passes into another's is ``\sigma`` times the share —
  and so do relative sizes, the margins of the processes, and where they run.
* **The detail of flows inside a group of regions matters little**, as long
  as the totals around it are right: the group's trade with the outside, and
  each member's production and use.

The [forest examples](@ref "The forest sector: complete data, and what calibration changes")
measure this on one economy and two scenarios. A reference run calibrated on
partial data — every flow touching the regions of one country imputed from
national totals — gives practically the results of one calibrated on complete
data. One left uncalibrated, with a single price per product, trade shares
from a rule and guessed costs, gets most of the results about right and
doubles the main one, because it places a logistic limit closer to binding
than the data do.

That says where to spend the effort:

| get right | can be imputed or approximated |
|:---|:---|
| production and use by region | flows between regions of one country — [`impute_flows`](@ref) |
| trade at the level customs report it | regional prices: the national price, or it adjusted by transport |
| a price level for each product, in each region where it differs | the recipes of the processes: published conversion factors, one region's input mix |
| the constraints your scenarios will run into: capacities, logistic limits, policies | consumer prices: derived from flows and producer prices |

[`impute_flows`](@ref) fills in a trade matrix with the flows closest to a
gravity pattern that meet every known total. It works best on the part that is
actually unknown. When customs report a country's trade with each partner but
not which of its regions it went through, impute only the country's corner of
the matrix: its regions with their production and use, and each partner
standing for its trade with the country as a whole — what it sent in, as an
origin, and what it took out, as a destination:

```julia
## for one product; `country` are the regions of the country, `partners` the rest
corner = impute_flows(
    production = merge(Dict(r => production[r] for r in country),
                       Dict(b => imports_from[b] for b in partners)),
    use        = merge(Dict(r => use[r] for r in country),
                       Dict(b => exports_to[b] for b in partners)),
    distance = distance, decay = 1.4, home_bias = 3.0,
    totals = Dict((partners, partners) => 0.0))   # their trade with each other is known
```

The rest of the matrix is what customs report between partners, and what each
partner keeps of its own production; the whole goes to
[`armington_calibration`](@ref) as if it had been observed. The
[second forest page](@ref "The forest sector: calibrating with partial data")
does exactly this.

Inconsistent data need no special treatment: the reference run *is* their
reconciliation, the one point where every balance holds, and its distance from
the data says where the data disagree. Only two things are hard errors: a
region exporting more than it produced, which leaves its own deliveries
negative, and a process with a negative margin, which the model would close.

### Never re-anchor on your own output

In a multi-period run the temptation is to re-base each period on the last:
set `q0` and `p0` to the quantities and prices the model just produced, and
carry on. **Do not.** Written out, that makes demand

```math
D_t = q_{t-1}\left(\frac{P_t}{P_{t-1}}\right)^{-\eta}
```

which is precisely the lagged-ratio form this package avoids, arrived at from
another direction. Its results depend on the path taken rather than the state
reached, its errors compound, and the elasticity you estimated is no longer
anchored at the level you estimated it for.

`p0` and `q0` stay at the calibration point for the whole run. Everything that
moves between periods moves through `shift` and the limits, computed from state
and from exogenous drivers. The discipline in one sentence:

> Use model output to update **physical state**, never **behavioural
> parameters**.

Updating a fish stock by `stock + growth − catch` is an accounting
identity, and the catch may perfectly well come from the model. Updating `q0`
to last period's quantity is a silent re-estimation of behaviour from your own
output, and is a different thing entirely.

The one re-anchoring that is legitimate is to **new data**: if observations for
a later year arrive and you want to forecast from there, that is a fresh
calibration against reality.

## Dynamics are somebody else's job

This package solves **one equilibrium, at one point in time**, and that is the
whole of its remit. There is no time index anywhere in it and there will not
be: a model that needs several periods runs `solve_market` once per period and
carries the state between them itself. Keeping the boundary there is what lets
the same engine sit under a fishery model with annual steps, a capacity-planning
study with five-year ones, and a one-shot policy comparison with no time at
all.

What follows is therefore advice to whoever writes that loop, not a description
of something `Markets.jl` does.

The approach worth taking is **recursive dynamics**: solve a full static
equilibrium each period and link periods only through physical state, updated
from the previous period's solution:

* resource stocks and their growth shift the primary supply curves;
* installed capacity bounds process activity;
* exogenous drivers (income, population, ...) shift the demand curves.

Demand and supply remain *level* functions of price in every period,
e.g. ``D_t(P) = \alpha_t P^{-\eta}``, with the shifters of
[Everything other than the price](@ref) carrying ``\alpha_t``. That is the
intended seam: the loop updates `shift` and the limits between periods and
never touches the fitted curves.

This deliberately avoids the form used in some recursive models, where
demand in a period is a function of demand and price in the previous period,
``D_t / D_{t-1} = (P_t / P_{t-1})^{\eta}``. That form is a finite-difference
proxy for an elasticity, so its results depend on the path taken and its errors
compound across periods. (The ratio form is at least dimensionless, so the
objection is about path dependence, not about units.) Anchoring each period on a static
equilibrium with level functions keeps the dynamics interpretable.

## Planned extensions

These are things the engine itself might grow. A time loop is not among them,
for the reason given under [Dynamics are somebody else's job](@ref).

* **Armington refinements**: a separate composite per user (final demand and
  each process) instead of one per region, as GTAP does for the
  domestic-versus-imported margin; and a CET nest on the export side, so that
  redirecting output between destinations is costly for producers too, which is
  the mirror image of the import-side nest built here.
* **Inputs with an exogenous price** inside a nest (e.g. a chemical bought at a world
  price). A `SupplyCurve` with a very high elasticity is an approximation
  available today.
* **Ad valorem taxes without iteration**, which would need the model restated
  as a complementarity problem rather than an optimisation.
* **Estimation** of the parameters from data.
