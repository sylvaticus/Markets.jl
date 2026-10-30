# Processes and substitution

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

A nest with several products is a constant-elasticity-of-substitution bundle,
the function of [Arrow, Chenery, Minhas & Solow (1961)](@ref acms1961). The aim is **smooth** substitution: if one input becomes more
expensive, the process uses proportionally less of it, rather than switching
entirely from one input to another as a linear model with alternative
processes would do. For instance, a steel mill that uses 30% scrap at equal
prices may use 20% scrap when scrap becomes more expensive, and the share moves
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
  away from an input whose relative price rises. At 1 the mix is Cobb–Douglas
  (value shares stay constant), below 1 the inputs are poor substitutes, at 0
  the mix does not move at all, and large values make the inputs nearly
  perfect substitutes.

The shares are written in calibrated form, so that the price of the CES
composite equals the input price when all input prices are equal. As a result
the equilibrium price of an output equals its true marginal cost of
production.

### Range of the elasticity of substitution

``\sigma`` may be any value ``\ge 0``, here and in every other CES of the
package — the [`Armington`](@ref) and [`OriginNest`](@ref) groups of origins
and the [`DemandNest`](@ref) groups of final demand:

* ``\sigma > 1``: good substitutes. When an input's relative price rises, its
  quantity falls by more than its price rose, so its value share falls.
* ``\sigma = 1``: Cobb–Douglas. Each input keeps the value share
  ``\delta_i`` whatever the prices.
* ``0 < \sigma < 1``: poor substitutes, or complements. The quantity falls by
  less than the price rose, so the value share of the dearer input rises.
* ``\sigma = 0``: fixed proportions (Leontief). Input ``i`` is ``\delta_i``
  of the composite whatever the prices.

The aggregator is written with ``\rho = (\sigma - 1)/\sigma``, which cannot
be evaluated at ``\sigma = 0`` (``\rho = -\infty``), at ``\sigma = 1``
(``\rho = 0``) or at ``\sigma = \infty``. At those three values the engine
writes its limit instead, and the price of the composite is the limit of the
CES index:

| ``\sigma`` | composite of the ``q_i`` | price of the composite |
|:-----------|:-------------------------|:-----------------------|
| ``0`` | ``\min_i q_i/\delta_i`` | ``\sum_i \delta_i p_i`` |
| ``1`` | ``\prod_i (q_i/\delta_i)^{\delta_i}`` | ``\prod_i p_i^{\delta_i}`` |
| ``\infty`` | ``\sum_i q_i`` | ``\min_i p_i`` |

These are the exact limits of the calibrated CES, so an elasticity just beside
0 or 1 gives nearly the same equilibrium as 0 or 1 itself. Every one of them is
concave, so the program stays convex over the whole range.

Two practical points. Below 1 the exponent ``\rho`` is negative, so
``q_i^{\rho}`` grows without bound as an input goes to zero; the quantity floor
keeps the solver away from zero, as it already does above 1. And very close to
0 — below about ``\sigma = 0.01``, depending on the data — ``q_i^{\rho}``
leaves the range of floating-point numbers and the solve fails. Use
``\sigma = 0`` there: at ``\sigma = 0.01`` the equilibrium is already within
a fraction of a percent of it.

A nest of a single product has nothing to substitute and is a fixed
coefficient whatever its ``\sigma``; [`leontief`](@ref) builds one and stores
``\sigma = 0``. A fixed-proportion requirement of several products is
`ces(...; sigma = 0)`.

### Value-added cost

`vacost` is the cost per unit of activity of everything not modelled as an
explicit input: labour, energy, capital, chemicals. It is a constant, so
these inputs have an exogenous price and cannot substitute for the modelled
inputs. It may differ by region, as wages and energy do, and a `vacost`
calibrated on base-year margins always will.

### Why processes, and not supply and demand per transformation

There is another way to write the same economy. Instead of a transformation
with fixed coefficients whose scale the material balances determine, each
industry could carry its own estimated **output supply** and **input demand**
functions — ideally derived from a profit function, so that Hotelling's lemma
makes them consistent — and prices would become unknowns of a square system
rather than duals of a balance. That is the behavioural or equation-system
paradigm, and the choice between the two is old enough to have a standard
reference in [McCarl & Spreen (1980)](@ref mccarl1980).

Which one a model uses is largely a matter of which literature it grew up in.
Forest-sector, land-use and bioeconomy models are mostly process-based, as the
review in [Latta et al. (2013)](@ref latta2013) shows; agricultural outlook and
trade-policy models are mostly behavioural, as the classification in
[van Tongeren et al. (2001)](@ref vantongeren2001) shows; and the large models
in both camps end up hybrids, because a behavioural model that does not tie its
input demands to its output supplies by a fixed yield can imply impossible
physical conversions.

This package takes the process route, and the reasons are worth stating plainly
because they are also its limitations.

**What the choice buys.** Every balance is an equality in physical units, and
a recipe is a fixed yield: what goes into a plant comes out of it, as product
or as residue, and no combination of estimated elasticities can violate that.
The one qualification is aggregation: a CES bundle of inputs, or an Armington
composite of origins, is an aggregate rather than a sum, and the physical
account then needs one more line — see [Physical accounts](@ref). The parameters are engineering
coefficients, which exist for most transformations, where econometric estimates
on intermediate markets usually do not. A technology with no history at all can
be added by writing down its recipe. And because the problem stays a single
convex optimisation, the equilibrium exists, the welfare measure means
something, and one solve delivers both.

**What it costs.** Constant returns to scale mean a process has no supply
response of its own. At the zero-profit price it would run at any scale, so its
level is pinned entirely by the balances around it — never by its own
profitability. A transformation in this model does not decide how much to
produce; the rest of the economy decides for it.

#### What that looks like in practice

Sharper than "the cheapest region takes everything". Take two otherwise
identical regions, each able to run the same process (ten units of input per
unit of output, a value-added cost of 1300, freight of 5 on both products), and
lower one region's value-added cost step by step:

| cost in A | output in A | output in B |
|----------:|------------:|------------:|
| 1300      | 79.8        | 79.8        |
| 1247      | 80.4        | 79.6        |
| 1245      | 104.6       | 55.4        |
| 1243      | 160.0       | **0.0**     |

Two failures, not one. Over the first 53 units of cost advantage — 4% —
**nothing happens at all**; then the whole industry moves over a change of 4
units, or 0.3%. Processing location is a step function: numb, then total.

The width of the numb zone is not arbitrary. A local processor is protected by
the freight on everything it would otherwise import, which is the freight on
the *inputs multiplied by the recipe coefficient*, plus the freight on the
output: ``10 \times 5 + 5 = 55`` here, and the step duly falls at
``1300 - 55``. Raise the input freight to 60 and the protection becomes 605: a
46% cost advantage then buys nothing whatever. In a sector whose raw material
is expensive to move and whose recipes take several units of input per unit of
output — ore smelting, sugar refining, forestry — a model is normally deep
inside the numb zone, so the practical symptom is inertia rather than
volatility.

#### Why this is usually not the problem it appears to be

Because [Armington](@ref "Imperfect substitution between origins (Armington)")
restores the missing slope. If a product's regional varieties are imperfect
substitutes, every region that makes it faces a downward-sloping demand for its
own variety, and the step disappears. The same experiment, with
``\sigma = 4`` on the output:

| cost in A | output in A | output in B |
|----------:|------------:|------------:|
| 1255      | 80.3        | 79.6        |
| 1243      | 80.4        | 79.6        |
| 1200      | 81.8        | 78.4        |
| 1100      | 85.7        | 75.2        |

Smooth, and gentle: a 15% cost advantage moves the split from 50/50 to 53/47.
The [forest example](@ref "The forest sector: France in the world") declares an
`Armington` specification for every product, so between its blocs it does
not exhibit the behaviour above at all. **Read the caveat as applying to
products you have deliberately left homogeneous**, which in practice means
intermediates for which no trade shares were available, and groups of regions
pooled as perfect substitutes — where, once costs are calibrated, it takes a
sharper form still.

#### Perfect substitutes and constant returns

Calibrating `vacost` on base-year margins puts every process that ran exactly
at break-even, in every region. Between regions whose products are perfect
substitutes — a homogeneous product, or an [`OriginNest`](@ref) with
`sigma = Inf` — that turns the step above into a **plateau**: the location of
processing inside the group is not decided at all.

The mechanism is simple. Take two regions of the group that both sell the
input *and* the output into the same market — two regions of one country
both exporting a raw material and the product made from it to the same
neighbour, say. Move one unit of processing from the first to the second: the
first exports the raw material it no longer processes, the second exports
that much less, and the product changes hands the other way. Both swaps
happen inside one market, at one delivered price, so the move costs exactly
nothing, and so does the next one. Nothing in the model prefers one split to
another; the solver stops wherever its path ends, and a cost difference of a
hair moves the whole plateau. In a test on the forest economy of the
[examples](@ref "The forest sector: France in the world") with its French
regions perfectly substitutable, a 0.1% change in the costs of the French
mills emptied one region of sawmilling and nearly doubled it in another.

Whether a base year sits on such a plateau depends on its trade pattern. When
each region of the group sells its output to a *different* market, moving
processing between them shifts quantities between markets, whose demands
slope, and the split is pinned; when they share a market for both input and
output, as regions of one country selling to their neighbours commonly do, it
is not.

There are three ways out, in increasing order of what they ask for:

* **A finite elasticity inside the group**, however high. Each region's
  product then has a demand of its own, and the plateau tilts into a slope.
  The cost is that the group is no longer one market: products cross-haul
  between its regions, and producer prices may differ by more than the
  freight between them. The forest examples do this, with `sigma = 30`, and
  show both.
* **A value-added cost rising with activity**, the principled remedy described
  in the next section. Each region's mill then has a supply curve of its own,
  so the split is decided by costs rather than left open, and the group can
  stay perfectly substitutable.
* **Capacities that bind at the base year.** They pin the split, but only if
  set equal to observed output, which says every plant ran flat out — rarely
  true, and it forbids every plant from growing in any scenario. A capacity
  with room to spare leaves the plateau in place inside its bounds.

#### The fix that is not implemented, and why

The principled remedy, for the step and for the plateau alike, is not to
abandon processes but to give them decreasing returns: replace the constant `vacost` with a cost rising in the activity
level, calibrated so that the observed activity is what a processor would
choose at observed margins. This is [positive mathematical
programming](@ref howitt1995), the standard answer in agricultural sector
models, and it would sit comfortably in this package — the objective term would
take the same shape as the supply curves already do, the problem would stay
convex, and an infinite elasticity would recover today's behaviour exactly.
Between that limit and a [`Limit`](@ref) on the activity, which is the
opposite limit of a perfectly inelastic process, the parameter would simply
interpolate.

It is not implemented because of what it would need as data, not because of
what it would need as code. Per process and per region it asks for a reference
activity, a base-year margin — the difference between the output price and the
input costs, which requires consistent price data on intermediate products in
every region, the thinnest data in the sector — and an elasticity of processing
supply, which is not an estimated quantity anywhere. Choosing it would mean
replacing capacity bounds that were invented with an elasticity that was
invented.

That trade may still be worth making, and for a model used to study where
industry locates it probably is: a slope is a weaker assumption than a bound,
it has an unambiguous direction, it can be swept in sensitivity analysis, and
it makes the base year an interior optimum rather than a corner held in place
by a constraint. But it should be made deliberately, and documented as a
calibration choice rather than presented as measurement. Two details to get
right if you do: the rising cost needs a positive intercept — an isoelastic
cost alone makes the first unit free, so every region would produce a little of
everything — and the inframarginal profit it creates has to enter the welfare
measure, or welfare stops being comparable across model versions.

A full switch to the behavioural paradigm is a different and much larger
decision. Independently estimated response functions are generally not
integrable, so no objective exists and the equilibrium has to be posed as a
complementarity problem and solved as one — [Mathiesen (1985)](@ref
mathiesen1985), [Rutherford (1995)](@ref rutherford1995),
[Ferris & Pang (1997)](@ref ferris1997), and [[C4]](@ref C4) and
[[C5]](@ref C5) for the tooling. That route would handle ad valorem wedges
exactly instead of by iteration, and would accept asymmetric cross-price
effects; it would also give up the guarantee of existence, the single welfare
number, and the mass-balance consistency that is the main reason this package
is built the way it is.

## Joint products and residues

A process can have several outputs. Residues (e.g. the slag of a blast furnace,
or the fuel oil of a refinery) are ordinary products with their own material balance in each region, so they can be
consumed, traded or used as inputs to other processes, and their price is
determined by the competition between their uses.

The material balance is an **equality**: every unit of a by-product must be
used, traded or consumed. There is no free disposal, and this has a consequence
worth stating plainly. A by-product that nobody wants is a *bad*: somebody must
be paid to take it, and its equilibrium price is **negative**. The model
reports it that way, since the price is the signed dual of the balance.

That is the right answer to the question as posed, but it is rarely the
question meant, and under a large supply shock it can take a local market
somewhere surprising. If free disposal is what you want instead, give the
by-product somewhere to go — [`disposal`](@ref) writes the process for you:

```julia
processes = [refinery, steel_mill,
             disposal(product = "fuel_oil"),          # free: price ≥ 0
             disposal(product = "slag", cost = 12)]   # paid: price ≥ −12
```

Its price then settles at zero rather than going negative, and a positive
`cost` represents a real disposal charge — a landfill fee, an incinerator, a
tailings dam. Deciding which by-products have a
way out is a modelling choice the model cannot make for you, which is why there
is no default. Writing the balances as
inequalities, as a mixed-complementarity formulation would, is the other route
to the same thing; it is not what this package does.

## Substitution in final demand

Final demands are independent of each other by default: each product's curve
responds to its own price alone. That is often too strong — butter and
margarine go on the same bread, petrol and diesel cars do the same job on the
road — and the question is how to let one product's price move another's
demand.

The tempting answer, a matrix of cross-price elasticities, is **not available**,
and for the same reason an arbitrary matrix of Armington elasticities is not.
Demand here is the gradient of the consumer-benefit term of the objective;
that is what makes the surplus-maximisation equivalence hold in the first
place. A demand system that is not the gradient of anything — an arbitrary
matrix is not, unless it satisfies the symmetry conditions — has no benefit
function behind it, and the program it would be put into would no longer be the
market it claims to represent.

What *is* available is the same answer as for origins: a **tree**, through
[`DemandNest`](@ref). Products that substitute closely sit in a group with a
high ``\sigma``; the group as a whole substitutes against the rest at the lower
``\sigma`` of the level above; and the whole bundle faces one
[`DemandCurve`](@ref) curve, which describes the composite. So

```julia
DemandCurve(region = "EU", p0 = 30000, q0 = 12, elasticity = 0.5,
           product = DemandNest(sigma = 1.6, products = ["electric_car"], shares = [0.25],
                       nests = [DemandNest(sigma = 5.0, share = 0.75,
                                           products = ["petrol_car", "diesel_car"],
                                           shares   = [0.65, 0.35])]))
```

says that petrol and diesel cars are close substitutes, that an electric car is
a more distant substitute for either, and that cars as a whole have an
own-price elasticity of 0.5.

The cross-price effects this implies are the CES ones: within a group, the
compensated cross-price elasticity of product ``i`` with respect to the price of
``j`` is ``\sigma s_j``, with ``s_j`` the value share of ``j``. A large
supplier therefore transmits its price more strongly than a small one at the
same ``\sigma`` — the shares matter as much as the elasticity, exactly as they
do between origins.

### Substitutes outside the model

Cars compete with public transport, whose fares this model does not
determine. That needs none of the machinery above: an exogenous price is a
constant, so its effect on demand is a **shift of the whole curve**. Put the
cross-price elasticity in [`exogenous_shift`](@ref) as the driver elasticity:

```julia
DemandCurve(product = cars, region = "EU", p0 = ..., q0 = ...,
           elasticity = 0.4,
           shift = exogenous_shift(fares = (level = 3.2, reference = 2.0,
                                            elasticity = 0.35)))
```

Fares 60% dearer than at calibration, with a cross-price elasticity of 0.35,
raise demand for the car bundle by ``1.6^{0.35} \approx 1.18``. The
endogenous substitution inside the bundle and the exogenous substitution
against public transport are different mechanisms, and it is worth keeping them visibly
apart.
