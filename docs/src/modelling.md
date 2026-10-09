# Modelling choices

This page explains what the model assumes and why. The exact program that the
code builds is on the [Code implementation](@ref) page.

## Partial, not general, equilibrium

The model represents one sector in detail (its products, processes and trade)
and treats the rest of the economy as exogenous: income, population or
construction activity enter only through the position of the demand curves,
and the prices of non-modelled inputs (labour, energy, ...) enter only through
the value-added cost of each process. A general-equilibrium model would add
the accounting of the whole economy without adding much detail to the markets
of the sector.

## The equilibrium as surplus maximisation

Under perfect competition, the competitive spatial equilibrium is the solution
of a single optimisation that maximises **net social surplus** (Samuelson,
1952; Takayama & Judge, 1971):

```math
\max \; \underbrace{\sum \int_0^{D} P^D(x)\,dx}_{\text{gross consumer benefit}}
     \;-\; \underbrace{\sum \int_0^{S} P^S(x)\,dx}_{\text{primary resource cost}}
     \;-\; \underbrace{\sum c_k z_k}_{\text{value-added cost}}
     \;-\; \underbrace{\sum \tau\, T}_{\text{transport cost}}
```

subject to a material balance for every product in every region. This is
equivalent to "clearing the markets": the first-order conditions of the
optimisation are the market-clearing conditions. One solve gives everything:

* the **quantities** (production, consumption, trade, process activity) are the
  primal solution;
* the **prices** are the dual variables of the material balances: the value of
  one more unit of a product available in a region.

There is no need to iterate on prices or quantities to clear markets.

### Where is the producer surplus?

The objective contains the producer surplus even though it is not written as a
separate term. Let ``\pi`` be the equilibrium prices. Adding and subtracting the
values ``\pi D`` and ``\pi S`` splits net social surplus into:

* **consumer surplus**, ``\int_0^D P^D - \pi D``, the area between the demand
  curve and the price;
* **primary producer surplus**, ``\pi S - \int_0^S P^S``, the area between the
  price and the supply curve (the rent of the resource owners, e.g. forest
  owners);
* **processors' profit**, value of outputs − value of inputs − value-added cost;
* **traders' profit**, ``(\pi_{to} - \pi_{from} - \tau)\,T`` on each route.

The material balances make all the ``\pi`` terms cancel, so the sum of the four
is exactly the objective. Processes have constant returns to scale and trade is
competitive, so at the equilibrium the last two are zero: an active process
earns exactly its costs, and goods are shipped only where the price gap equals
the transport cost. All producer surplus is therefore the rent of primary
supply, which is in the objective as the area above the supply curves.

## Demand and supply curves

Final demand and primary supply have constant elasticity, calibrated to pass
through a reference point ``(q_0, p_0)``:

```math
P^D(D) = p_0 \left(\frac{D}{q_0}\right)^{-1/\eta}, \qquad
P^S(S) = p_0 \left(\frac{S}{q_0}\right)^{1/\varepsilon}
```

Both elasticities need only be **positive**. Demand may be inelastic, which
matters because that is how demand for most forest products is estimated:
earlier versions of this package required ``\eta > 1``, on the grounds that for
``\eta \le 1`` the area under the demand curve from zero is infinite. That
argument does not hold. The objective needs an *antiderivative* of the inverse
demand, not the integral from zero, and the additive constant — infinite though
it is — changes neither the optimum nor the first-order conditions, hence
neither quantities nor prices:

```math
\int^D P^D = \frac{a}{1-1/\eta}\, D^{\,1-1/\eta} \quad (\eta \ne 1),
\qquad \int^D P^D = a \log D \quad (\eta = 1)
```

With ``\eta < 1`` the coefficient and the exponent are both negative, so the
term is still increasing and concave and the program is still convex. The one
thing lost is the *level* of welfare: with ``\eta \le 1`` the objective value is
not consumer surplus in absolute terms. Differences between scenarios remain
well defined, which is what a comparative exercise uses.

!!! warning
    Inelastic demand changes results substantially, and in the direction people
    often forget: with ``\eta`` well below 1 a supply shock is absorbed by
    *prices* rather than by quantities. Using ``\eta > 1`` for a product whose
    demand is really inelastic makes a model systematically under-react in
    price and over-react in volume.

Each final demand depends on its own price only; there are no cross-price
effects between final products. Substitution between products happens on the
input side, inside processes. This separability is also what makes the
surplus-maximisation equivalence legitimate: the demand system is integrable by
construction.

Each final demand depends on its own price only; there are no cross-price
effects between final products. Substitution between products happens on the
input side, inside processes.

## Processes, nests and input substitution

A **process** is a transformation with constant returns to scale: per unit of
activity it needs a fixed set of input requirements (its **nests**) and yields
fixed quantities of one or more outputs.

* Inputs in **different nests** are needed in fixed proportion to each other.
  If a process needs 1 unit of nest A and 0.2 units of nest B, it cannot use
  more of A to save on B.
* Inputs **in the same nest** substitute for each other.

### Leontief nests

A nest with a single product is a fixed coefficient (Leontief): exactly
`coeff` units per unit of activity. Fixed proportions between several products
are written as several Leontief nests.

### CES nests

A nest with several products is a constant-elasticity-of-substitution (CES)
bundle. The aim is **smooth** substitution: if one input becomes more
expensive, the process uses proportionally less of it, rather than switching
entirely from one input to another as a linear model with alternative
processes would do. For instance, a mill that uses 30% chips at equal prices
may use 20% chips when chips become more expensive, and the share moves
continuously with relative prices.

The cost-minimising mix satisfies

```math
\frac{q_i}{q_j} = \frac{\delta_i}{\delta_j}\left(\frac{p_j}{p_i}\right)^{\sigma}
```

where:

* the **shares** ``\delta_i`` (normalised to sum to 1) fix the mix *at equal
  input prices*: at equal prices, input ``i`` makes up the fraction
  ``\delta_i`` of the nest, both in quantity and in value. They are the
  calibration of where the substitution starts from. Without them all inputs
  would be used in equal amounts at equal prices.
* the **elasticity of substitution** ``\sigma`` fixes how fast the mix moves
  away from an input whose relative price rises. Values just above 1 give a mix
  close to Cobb–Douglas (value shares stay nearly constant); large values make
  the inputs nearly perfect substitutes.

The shares are written in calibrated form, so that the price of the CES
composite equals the input price when all input prices are equal. As a result
the equilibrium price of an output equals its true marginal cost of
production.

### Range of the elasticity of substitution

[`ces`](@ref) requires ``\sigma > 1``. Mathematically, CES is also defined for
``0 < \sigma < 1`` (inputs that are poor substitutes, or complements). The
code restricts it because the aggregator is evaluated as
``q_i^{\rho}`` with ``\rho = (\sigma - 1)/\sigma``: for ``\sigma < 1``,
``\rho`` is negative and ``q_i^{\rho}`` grows without bound as an input
goes to zero, which is numerically harder for the solver. ``\sigma = 1``
(Cobb–Douglas) is a singular case of this formula (``\rho = 0``).

[`leontief`](@ref) stores ``\sigma = 0``, but that value is never used: a
single-product nest has nothing to substitute, and the code treats any
single-product nest as fixed-coefficient regardless of ``\sigma``.

### Value-added cost

`vacost` is the cost per unit of activity of everything not modelled as an
explicit input: labour, energy, capital, chemicals. It is a constant, so
these inputs have an exogenous price and cannot substitute for the modelled
inputs.

## Joint products and residues

A process can have several outputs. Residues (e.g. sawmill chips) are ordinary
products with their own material balance in each region, so they can be
consumed, traded or used as inputs to other processes, and their price is
determined by the competition between their uses.

The material balance is an **equality**: every unit of a by-product must be
used, traded or consumed. There is no free disposal, and this has a consequence
worth stating plainly. A by-product that nobody wants is a *bad*: somebody must
be paid to take it, and its equilibrium price is **negative**. The model
reports it that way, since the price is the signed dual of the balance.

If free disposal is what you want instead, add a process that consumes the
by-product and produces nothing:

```julia
Process(name = :dump, inputs = [leontief(product = :chips, coeff = 1.0)],
        outputs = Pair{Symbol,Float64}[], vacost = 0.0)
```

Its price then settles at zero rather than going negative, and a positive
`vacost` on it represents a real disposal cost instead. Writing the balances as
inequalities, as a mixed-complementarity formulation would, is the other route
to the same thing; it is not what this package does.

## Trade between regions

Regions trade the *same physical product* at a cost per unit of quantity
shipped. This is the spatial price equilibrium:

* prices of a product in two regions differ by **at most** the transport cost
  between them;
* on a route that is used, the price gap **equals** the transport cost.

Transport costs are per unit of quantity because shipping a m³ of logs costs
roughly the same whether the logs are cheap or expensive. Bulky, low-value
products therefore have high transport costs relative to their value, which
limits their trade.

Since products are identical whatever their origin, two regions never ship the
same product to each other — there is no cross-hauling — and each region buys
from the cheapest source available to it. Trade patterns are therefore sharper
than observed ones: a small cost advantage can capture an entire market. This
is the default; the next section is the way out of it.

## Imperfect substitution between origins (Armington)

Observed trade does not look like that. Countries import and export the same
product at the same time, market shares move gradually when prices change, and
a region's own product keeps a share of its home market even when it is not
the cheapest. The standard answer (Armington, 1969) is to treat the varieties
of a product coming from different origins as **imperfect substitutes**:
"paper from Europe" and "paper from Asia" are related but distinct goods.

Listing a product in the [`Armington`](@ref) specifications of the economy
turns this on for it. Each destination `r` then buys quantities ``X_o`` of the
variety of every origin `o` it can buy from (its own, plus every origin with a
transport route into it) and uses them through a CES composite

```math
A_r \le \left( \sum_o \delta_{o,r}^{1/\sigma}\, X_{o,r}^{\rho} \right)^{1/\rho},
\qquad \rho = \frac{\sigma - 1}{\sigma}
```

which is the same aggregator as a CES input nest, applied to origins instead of
input products. Final demand and the processes of region `r` then draw on the
composite ``A_r``, not on the individual varieties.

Three things change.

**There are now two prices per product and region.** Local producers are paid
``\pi^s_r``, the value of their own variety, while local users pay ``\pi^c_r``,
the price of the composite, which is the CES price index of the delivered
prices of all origins:

```math
\pi^c_r = \left( \sum_o \delta_{o,r} \left(\pi^s_o + \tau_{o,r}\right)^{1-\sigma} \right)^{1/(1-\sigma)}
```

The two are reported side by side in the `prices` table as `producer_price` and
`price`. They differ because a region blends its own variety with imported
ones: an expensive producer can keep selling at a high price while its
customers pay much less on average.

**Regions cross-haul, and their prices are linked everywhere.** Because each
destination wants some of every variety, a region usually imports and exports
the same product at once, and market shares respond smoothly to relative
prices with elasticity ``\sigma``. Prices are no longer tied together by the
transport cost: a shift in supply or demand anywhere passes into every other
region's composite price, with a strength set by ``\sigma`` and by the shares,
and producer prices can differ across regions by far more than the cost of
shipping between them.

**The shares matter, and they are not the observed shares.**
``\delta_{o,r}`` is the share origin `o` would have in destination `r` *if all
delivered prices were equal*, so it carries the home bias and the trade pattern
that prices alone do not explain. Base-year delivered prices are never equal —
transport alone sees to that — so feeding an observed trade matrix in as
``\delta`` does **not** reproduce the base year. Calibration means inverting
the CES demand condition ``X_o = A \delta_o (P^c/q_o)^\sigma``:

```math
\delta_{o,r} \;\propto\; X_{o,r}\, q_{o,r}^{\,\sigma}
\qquad\text{equivalently}\qquad
\delta_{o,r} \;\propto\; s_{o,r}\, q_{o,r}^{\,\sigma-1}
```

from base-year quantities ``X`` or value shares ``s`` and delivered prices
``q``. [`armington_shares`](@ref) does this. The gap is not small: an origin
with an 80% observed share at a below-average delivered price calibrates to a
``\delta`` noticeably under 80%, because part of its success is explained by
its price rather than by preference. Left unspecified, every origin available to
a destination gets an equal share, which is neutral but rarely realistic.

Two consequences of calibrated shares are worth knowing:

* **An origin with no share is not in the composite.** A route that carried
  nothing in the base year carries nothing ever after, whatever prices do. This
  is the well-known small-shares problem of Armington models: they cannot
  create new trade relationships, only grow and shrink existing ones. Give a
  plausible small share to a route you want to be able to open.
* **The composite is not the physical sum.** ``A_r`` is a CES aggregate of the
  varieties, and since the shares sum to one it is at most their sum, with
  equality only when the mix is the one the shares describe. So reported
  consumption of an Armington product is a composite index, not a tonnage, and
  the m³ or tonne identity that holds for a homogeneous product does not hold
  for it. Production, trade and the balances are in physical units throughout;
  it is only the composite that is an index.

### Shares and elasticity answer different questions

The two parameters are easy to confuse, and they pull in different directions:

* ``\delta_{o,r}`` sets **how much** of origin `o` a destination takes when
  delivered prices are equal. A regulatory or commercial penalty against an
  origin — a standard its products meet only partly, a certification few of its
  producers hold — belongs here, as a low share.
* ``\sigma`` sets **how readily buyers switch** when relative prices move. It is
  about technical interchangeability, not about preference.

A low ``\sigma`` is therefore not a barrier against an origin. It means buyers
*cannot get away from it*: as ``\sigma \to 1`` the value shares stay at
``\delta`` whatever the prices, so an expensive origin keeps its share instead
of losing it. To represent a barrier, lower its share, or add its compliance
cost to the transport cost of that route; to represent products that are hard
to swap, lower the elasticity.

### Groups of origins

One ``\sigma`` per destination still says that every origin substitutes equally
well for every other. Real markets are not like that: two regions whose
products share grades, standards or certification are near-interchangeable,
while a third whose products need re-grading, or meet the importer's
non-tariff measures only partly, is a distant substitute for both.

A single CES cannot express this — equal pairwise substitution is what CES
*means*. The structure that can is a **nest of nests**: group the origins that
are interchangeable, give the group its own higher ``\sigma``, and let the group
as a whole substitute against the rest at the lower ``\sigma`` of the level
above. In this package that is an [`OriginNest`](@ref) inside an
[`Armington`](@ref):

```julia
Armington(product = :sawn_sw, destination = :EU, sigma = 2.5,
          nests = [OriginNest(sigma = 12, origins = [:EU, :NA])])
```

Buyers in the EU swap EU for NA sawnwood readily (``\sigma = 12``) and replace
either with Asian sawnwood only slowly (``\sigma = 2.5``).

This generalises: the groups may hold groups, and the structure is a tree whose
leaves are origins. The substitution between any two origins is then the
elasticity of the smallest group containing both. That is not an arbitrary
choice of representation, it is the **only** consistent one. An arbitrary matrix
``\sigma_{o,o'}`` is not integrable to any cost function: if two origins are
perfect substitutes they are the same good, so they must substitute identically
against any third, and more generally, among any three origins the two that
join higher in the tree must share the same elasticity. A tree enforces that by
construction; a matrix does not.

Two consequences worth knowing:

* A group whose ``\sigma`` equals its parent's is a no-op — it says the
  members are no closer to each other than to the outside. The package requires
  a group to be at least as substitutable inside as outside, and rejects the
  reverse as a specification error.
* As a group's internal ``\sigma`` grows, its members' prices are pulled
  together and the group behaves as one pooled variety: the Samuelson case,
  recovered inside a nest. At ``\sigma = \infty`` it *is* that case — the
  members are added up, the group enters the level above at the cheapest
  delivered price of its members, and between themselves they neither
  cross-haul nor differ in price by more than the freight. This is how a set of
  regions that are one market internally, such as the regions of a country, sit
  inside a world model as a single distinct variety.

### Different markets, different structures

Technical standards and non-tariff measures are set by the **importer**, so
they are asymmetric: the EU may treat Asian sawnwood as a distant substitute
while Asian buyers treat European sawnwood as a close one. A specification
therefore applies to the destinations given in its `destination` field, with an
`:all` entry as the fallback, and the elasticity, the groups and the shares may
all differ between importing markets.

A product is nonetheless either homogeneous everywhere or imperfectly
substitutable everywhere: it has one variety per origin, or none. And a region
that cannot produce a product at all — no supply curve, and no process making
it there — has no variety of it: it is nobody's origin, not even its own, and
its `producer_price` is reported as `missing` while it buys at the composite
price like anyone else.

### The limit ``\sigma \to \infty``

As ``\sigma`` grows the composite becomes a plain sum of the varieties and its
price index tends to the cheapest delivered price: the model returns to the
homogeneous, Samuelson spatial price equilibrium of the previous section. This
is not an approximation — at ``\sigma = \infty``, which is the default for
every product not listed, the engine builds exactly the homogeneous
formulation, with one balance and one price per region and product.

So the two regimes are the ends of one scale, and ``\sigma`` is where a product
sits on it: low ``\sigma`` for goods whose origin matters to buyers (branded
or quality-differentiated products, appearance-grade timber), high ``\sigma``
for commodities that are graded and interchangeable (pulp, chips, roundwood of
a given species).

### What it costs

The Armington formulation is bigger: a variety variable per (product, origin,
destination), a composite variable and two balances per (product, region),
instead of one trade variable per route and one balance, plus one further
variable and constraint per group of origins. And ``\sigma``, the groups and
the shares are extra parameters to estimate, which is why they are opt-in per
product rather than global.

## What the model does not constrain

Two limitations are inherent to the formulation rather than to the data, and
both bear on how results should be read.

**Processing location is sharp.** Processes have constant returns to scale and
a constant value-added cost, with no upper bound on activity. Nothing therefore
stops the whole of a process from concentrating in whichever region is
cheapest: location is bang-bang, the same sharpness the homogeneous trade
assumption produces, and Armington on the *outputs* only softens it indirectly
by making each region's variety wanted somewhere. Until capacity bounds exist
(they are in [Planned extensions](@ref)), read the geography of *processing* far
more cautiously than that of prices or of final demand. Restricting a process
to the regions that actually have mills, through the `regions` argument of
[`Process`](@ref), is the blunt instrument available today.

**Prices come from a continuous relaxation.** All quantities are continuous, so
the model knows nothing of indivisible investments, minimum efficient scale or
shutdown decisions.

## Dynamics

The model is static: one call to [`solve_market`](@ref) computes one equilibrium.
The intended approach to multi-period runs is **recursive dynamics**: solve a
full static equilibrium each period and link periods only through physical
state, updated from the previous period's solution:

* resource stocks and their growth shift the primary supply curves;
* installed capacity bounds process activity;
* exogenous drivers (income, population, ...) shift the demand curves.

Demand and supply remain *level* functions of price in every period,
e.g. ``D_t(P) = \alpha_t P^{-\eta}``.

This deliberately avoids the form used in some recursive models, where
demand in a period is a function of demand and price in the previous period,
``D_t / D_{t-1} = (P_t / P_{t-1})^{\eta}``. That form is a finite-difference
proxy for an elasticity, so its results depend on the path taken and its errors
compound across periods. (The ratio form is at least dimensionless, so the
objection is about path dependence, not about units.) Anchoring each period on a static
equilibrium with level functions keeps the dynamics interpretable.

## Previous formulation (v0.0.1)

Version 0.0.1 used a different formulation: constant-elasticity supply and
demand functions with cross-price elasticities, CES aggregation of supply and
demand over regions (Armington-type trade), and price–quantity equations solved
as a system of equations under a least-squares objective. That model never
solved ("too few degrees of freedom"). Version 0.0.2 replaces it with the
surplus-maximisation formulation above, in which market clearing comes from the
optimality conditions instead of being imposed equation by equation. The old
code remains in the git history.

## Planned extensions

* **Recursive dynamics**: a time loop around [`solve_market`](@ref), as
  described above.
* **Armington refinements**: a separate composite per user (final demand and
  each process) instead of one per region, as GTAP does for the
  domestic-versus-imported margin; and a CET nest on the export side, so that
  redirecting output between destinations is costly for producers too, which is
  the mirror image of the import-side nest built here.
* **Exogenous drivers** of demand and supply (income, population, resource
  availability) as explicit shifters.
* **Inputs with an exogenous price** inside a nest (e.g. a resin with a world
  price). A `SupplySpec` with a very high elasticity is an approximation
  available today.
* **Policies**: taxes, subsidies, tariffs and quotas.
* **Estimation** of the parameters from data.

## References

* Armington, P. S. (1969). A theory of demand for products distinguished by
  place of production. *IMF Staff Papers* 16(1).
* Buongiorno, J., Zhu, S., Zhang, D., Turner, J. & Tomberlin, D. (2003).
  *The Global Forest Products Model.* Academic Press.
* Samuelson, P. A. (1952). Spatial price equilibrium and linear programming.
  *American Economic Review* 42(3).
* Takayama, T. & Judge, G. G. (1971). *Spatial and Temporal Price and
  Allocation Models.* North-Holland.
