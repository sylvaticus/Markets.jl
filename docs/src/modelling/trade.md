# Taxes and trade

## Taxes, subsidies and other wedges

A tax drives a wedge between what one side of a transaction pays and what the
other receives. [`Tax`](@ref) places one on production, on consumption, on a
transformation or on a trade route, as an amount per unit or as a fraction, and
a subsidy is a negative rate.

In a surplus-maximising program a wedge enters exactly as a cost does: a tax on
a route is indistinguishable, to the optimiser, from a more expensive route.
The first-order conditions of the taxed program are then the equilibrium
conditions of the taxed market, so the quantities and prices that come out are
the right ones. A specific tariff of 10 on a route with a freight cost of 5
opens a price gap of exactly 15; a value-added tax of 25% puts the consumer
price 25% above the market price and moves the quantity accordingly.

!!! warning "The objective is no longer welfare"
    The program treats the tax as a cost consumed, when it is a transfer
    received, so the objective of a taxed model is not welfare and comparing
    it with an untaxed one overstates the loss. With per-unit rates it falls
    short by exactly the revenue; a relative rate on supply or demand scales a
    whole curve, and the gap is then not the revenue. `Results.welfare` is
    welfare in every case — the surplus at the solution, without the taxes —
    and `Results.taxes` reports what each instrument raises; who ultimately
    gets it is outside the model.

### Why ad valorem rates on trade are special

Three of the four bases take an ad valorem rate exactly, because the thing the
rate multiplies is known before the model is solved: a tax on production scales
the supply curve, one on consumption scales the demand curve, one on a process
scales its value-added cost. A tariff is different — it is a fraction of the
price **at the origin**, which is a dual of the problem and therefore not
available to the objective that produces it.

The honest routes out are to restate the model as a complementarity problem,
where prices are variables — the approach of [Mathiesen (1985)](@ref mathiesen1985)
and [Rutherford (1995)](@ref rutherford1995), available in Julia as
[[C4]](@ref C4) on the solver [[C5]](@ref C5) — or to iterate. This package iterates:
[`solve_market`](@ref) solves repeatedly, each time holding the tariff at the
money equivalent implied by the previous solution's producer prices, and stops
when the two agree. The fixed point is the ad valorem equilibrium — at it, the
price gap on a used route is freight plus the rate times the origin price,
which is what a tariff means — and everything else in the model is still solved
exactly, in one pass.

### Where along the chain to tax

An externality carried by the goods — the carbon in wood, the material that
ends up as waste — can be taxed where it arises, at the **source**, or where
it ends up, embodied in what is consumed, at the **end of the chain**.
[Gerlagh, Lorang, Pommeret & Pottier (2026)](@ref glpp2026) show that the two
are equivalent only if the prices of goods can depend on their embodied
content, so that a buyer facing the end-of-chain tax pays less for a cleaner
unit of the same product and the producer is rewarded for making it. Where
buyers see only the average content of a product group, the end-of-chain tax
falls on that average, producers have no reason to change their recipes, and
the source tax is the more efficient instrument — by a factor of ten for iron,
in the calibrated world model of [Lorang & Gerlagh (2026)](@ref lorang2026).

In this package a [`Tax`](@ref) is levied per unit, or as a fraction of the
value, of a product: a tax on primary supply (`on = "supply"`) is a source
tax, and behaves as the efficient one; a tax on consumption (`on = "demand"`)
taxes the product and not its content, which is the group-average case. Where
recipes are fixed the two coincide, as the content per unit of every good is
fixed too; where processes choose their inputs from a [`ces`](@ref) nest, or
users their origins from an [`Armington`](@ref) composite, they do not. A
consumption tax is then a second-best proxy for pricing the content, and the
model will show it as such.

## Trade between regions

Regions trade the *same physical product* at a cost per unit of quantity
shipped. This is the spatial price equilibrium:

* prices of a product in two regions differ by **at most** the transport cost
  between them;
* on a route that is used, the price gap **equals** the transport cost.

Transport costs are per unit of quantity because shipping a tonne of ore costs
roughly the same whether the ore is rich or poor. Bulky, low-value products —
gravel, cement, ore — therefore have high transport costs relative to their value, which
limits their trade.

Since products are identical whatever their origin, two regions never ship the
same product to each other — there is no cross-hauling, though a region may
carry it through on the cheapest path between two others — and each region buys
from the cheapest source available to it. Trade patterns are therefore sharper
than observed ones: a small cost advantage can capture an entire market. This
is the default; the next section is the way out of it.

## Imperfect substitution between origins (Armington)

Observed trade does not look like that. Countries import and export the same
product at the same time, market shares move gradually when prices change, and
a region's own product keeps a share of its home market even when it is not
the cheapest. The standard answer, [Armington (1969)](@ref armington1969), is
to treat the varieties
of a product coming from different origins as **imperfect substitutes**:
"wine from France" and "wine from Chile" are related but distinct goods.

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
  varieties, and with shares summing to one it is at most their sum, with
  equality only when the mix is the one the shares describe. Production, trade
  and the balances are in physical units throughout; it is only the composite
  that is an index. The `scale` of an [`Armington`](@ref) specification
  multiplies the composite, and [`armington_calibration`](@ref) sets it so
  that the composite equals the physical sum at the base year: base-year data
  in m³ then balance exactly, and the composite is priced at the average
  delivered unit value. Away from the base year the composite drifts from the
  sum again, to either side, as any CES does; the results report both — see
  [Physical accounts](@ref).

### Shares and elasticity answer different questions

The two parameters are easy to confuse, and they pull in different directions:

* ``\delta_{o,r}`` sets **how much** of origin `o` a destination takes when
  delivered prices are equal. A regulatory or commercial penalty against an
  origin — a standard its products meet only partly, a certification few of its
  producers hold — belongs here, as a low share.
* ``\sigma`` sets **how readily buyers switch** when relative prices move. It is
  about technical interchangeability, not about preference.
  [McDaniel & Balistreri (2003)](@ref mcdaniel2003) review what it has actually
  been estimated to be, and how widely those estimates differ — worth reading
  before settling on a number.

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
*means*. The structure that can is a **nest of nests**, after
[Sato (1967)](@ref sato1967): group the origins that
are interchangeable, give the group its own higher ``\sigma``, and let the group
as a whole substitute against the rest at the lower ``\sigma`` of the level
above. In this package that is an [`OriginNest`](@ref) inside an
[`Armington`](@ref):

```julia
Armington(product = "cars", destination = "EU", sigma = 2.5,
          nests = [OriginNest(sigma = 12, origins = ["EU", "US"])])
```

Buyers in the EU swap EU for US cars readily (``\sigma = 12``) and replace
either with Chinese cars only slowly (``\sigma = 2.5``).

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
they are asymmetric: the EU may treat Chinese cars as a distant substitute
while Chinese buyers treat European cars as a close one. A specification
therefore applies to the destinations given in its `destination` field, with
one for every destination (`:`) as the fallback, and the elasticity, the groups
and the shares may
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
or quality-differentiated products: cars, wine, machinery), high ``\sigma``
for commodities that are graded and interchangeable (crude oil of a given
grade, wheat, copper).

### What it costs

The Armington formulation is bigger: a variety variable per (product, origin,
destination), a composite variable and two balances per (product, region),
instead of one trade variable per route and one balance, plus one further
variable and constraint per group of origins. And ``\sigma``, the groups and
the shares are extra parameters to estimate, which is why they are opt-in per
product rather than global.
