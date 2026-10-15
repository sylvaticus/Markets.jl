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
of a single optimisation that maximises **net social surplus**. The problem was
posed by [Enke (1951)](@ref enke1951) — who solved it by building an electrical
circuit — turned into an optimisation by [Samuelson (1952)](@ref samuelson1952),
and given the programming formulation everything since has used by
[Takayama & Judge (1964)](@ref takayama1964). The same formulation, applied to
the global forest sector, is what [[C3]](@ref C3) solves:

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
matters because that is how demand for most forest products is estimated — see
[Buongiorno (2015)](@ref buongiorno2015) for the orders of magnitude:
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

A demand curve on its own depends on its own price alone. Products that compete
for the same use are instead grouped under one curve, which is what gives
cross-price effects between them — see
[Substitution in final demand](@ref). Either way the demand system is the
gradient of the objective's benefit term, and that integrability is what makes
the surplus-maximisation equivalence legitimate in the first place; it is also
why the cross-price structure has to be a tree rather than a free matrix.

## Everything other than the price

A curve says how much is offered or wanted *at a given price*. Its position
depends on everything else: how much timber is standing, how dense the forest
road network is, how many crews and machines there are, how many people live in
the region and what they earn. Those are exogenous to a partial-equilibrium
model and have to enter it somehow. There are three standard devices, and they
are not interchangeable.

### 1. A multiplicative shifter — the workhorse

Scale the whole curve by a factor built from the drivers and the elasticity of
the curve with respect to each:

```math
S(P) = \underbrace{\prod_j \left(\frac{Z_j}{Z_j^0}\right)^{\gamma_j}}_{\text{shift}}
       \; q_0 \left(\frac{P}{p_0}\right)^{\varepsilon}
```

This is what [`SupplySpec`](@ref)'s and [`DemandSpec`](@ref)'s `shift` field
holds, and what [`exogenous_shift`](@ref) builds. It is the device used by
essentially every forest-sector partial-equilibrium model for the link between
the forest and the harvest, and there are three good reasons for that.

It **comes free with the estimation**. The regression that gives the price
elasticity — ``\log S`` on ``\log P`` — gives the ``\gamma_j`` in the same
step if the drivers are in it. No separate exercise is needed.

It is **safe**: a positive factor times a positive curve stays positive, at any
level of any driver. An additive shifter, ``S = q_0 (P/p_0)^\varepsilon + A``,
can drive supply negative and destroys the constant-elasticity reading of
``\varepsilon``. Use addition only where something genuinely arrives in a
fixed quantity regardless of price — salvage after a storm, recovered paper
collected under a mandate, a quota of imports from outside the model. In this
package that is a separate inelastic [`SupplySpec`](@ref) with a `capacity`,
not a shifter.

It **separates calibration from scenario**. `p0` and `q0` keep holding the
base-year point the curve was fitted to; the scenario lives in `shift`. Rewrite
`q0` instead and the calibration is gone, which matters most in a
recursive-dynamic run, where the time loop updates the shifter from the state of the
forest each period and must not touch the fitted curve.

!!! note "Quantity shifters and cost shifters are the same thing here"
    A road that opens a stand does not create timber, it lowers the cost of
    reaching it — a *cost* shifter, not a quantity one. With an iso-elastic
    curve the distinction collapses: multiplying the quantity by ``m`` is
    exactly multiplying the marginal cost by ``m^{-1/\varepsilon}``, since
    ``b = p_0 (m q_0)^{-1/\varepsilon}``. With ``\varepsilon = 0.5``, offering
    10% more wood at every price is the same as harvesting 17% more cheaply at
    every volume. What matters is only that the elasticity you apply is the one
    your ``\gamma`` was estimated against.

### 2. A ceiling — for a limit that binds

A shifter scales a curve; it never stops it. When the constraint is a genuine
limit — an allowable annual cut, a licensed quota, the capacity of the crews
and machines that exist this year — the right device is a bound:

```math
S \le \bar S(Z)
```

which is [`SupplySpec`](@ref)'s `capacity`. The difference is not cosmetic. At
a binding ceiling the price parts company with the marginal cost on the curve,
and the gap is the **scarcity rent** of the resource — the economically
interesting quantity, and the one a pure shifter cannot produce. Doubling a
shifter doubles supply at every price; a ceiling says that beyond a point no
price will buy more.

The rent is the price producers get minus what the curve says the last unit
cost, taxes on production included:

```julia
rent = producer_price - ((1 + t_rel) * s.p0 * (S / (s.shift * s.q0))^(1 / s.elasticity) + t_abs)
```

which is the multiplier of the ceiling, and is how `Results.supply` reports it:
zero whenever the ceiling does not bind.

A ceiling set above what would have been produced anyway changes nothing, so it
is safe to carry one on every curve that has a known physical limit.

### 3. A structural supply — when the factors are the point

The two devices above are reduced-form: the curve is estimated, and the drivers
move it. The alternative is to derive supply from the technology instead, by
making the resource a product and harvesting a process:

* the standing stock becomes a product with its own curve or endowment, bounded
  by `capacity`;
* harvesting becomes a [`Process`](@ref) that consumes standing timber and
  whatever else the model prices, and yields roundwood;
* labour, machinery and access then enter as priced inputs or as bounds, and
  the *shape* of the resulting supply curve is derived rather than assumed.

This is how a computable general equilibrium model treats factors of
production, and it is the honest route when the question is about the factors
themselves — the effect of a labour shortage, of a machinery subsidy, of a road
programme. The price is the data: factor shares and factor prices by region,
which a reduced-form shifter does not need.

### Choosing

| Driver | Device |
|:---|:---|
| Growing stock, forest area | shifter, ``\gamma`` typically well below 1 |
| Allowable cut, quota | ceiling |
| Road density, accessibility | shifter (a cost shifter in disguise) |
| Harvesting crews and machines | ceiling, or structural if that is the question |
| Income, population | shifter on demand |
| Storm salvage, recovered fibre | a separate inelastic curve with a ceiling |

The sensible default for a sector model is the first column of devices, with
ceilings wherever a physical limit is known, and the structural route reserved
for the questions that are really about factors.

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

**What the choice buys.** Mass and yield consistency hold by construction: what
goes into a mill comes out of it, as product or as residue, and no combination
of estimated elasticities can violate that. The parameters are engineering
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
output — forestry, exactly — a model is normally deep inside the numb zone, so
the practical symptom is inertia rather than volatility.

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
`Armington` specification for every product, so it does not exhibit the
behaviour above at all. **Read the caveat as applying to products you have
deliberately left homogeneous**, which in practice means intermediates for
which no trade shares were available.

#### The fix that is not implemented, and why

The principled remedy is not to abandon processes but to give them decreasing
returns: replace the constant `vacost` with a cost rising in the activity
level, calibrated so that the observed activity is what a processor would
choose at observed margins. This is [positive mathematical
programming](@ref howitt1995), the standard answer in agricultural sector
models, and it would sit comfortably in this package — the objective term would
take the same shape as the supply curves already do, the problem would stay
convex, and an infinite elasticity would recover today's behaviour exactly.
Between that limit and `Process`'s `capacity`, which is the opposite limit of a
perfectly inelastic process, the parameter would simply interpolate.

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

A process can have several outputs. Residues (e.g. sawmill chips) are ordinary
products with their own material balance in each region, so they can be
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
processes = [sawmill, panelmill,
             disposal(product = "chips"),                   # free: price ≥ 0
             disposal(product = "blackliquor", cost = 12)]   # paid: price ≥ −12
```

Its price then settles at zero rather than going negative, and a positive
`cost` represents a real disposal charge — a landfill fee, an incinerator, the
cost of leaving residues on the forest floor. Deciding which by-products have a
way out is a modelling choice the model cannot make for you, which is why there
is no default. Writing the balances as
inequalities, as a mixed-complementarity formulation would, is the other route
to the same thing; it is not what this package does.

## Substitution in final demand

Final demands are independent of each other by default: each product's curve
responds to its own price alone. That is often too strong — pellets and
firewood are the same heat, sawn softwood and panels do the same job in a
building — and the question is how to let one product's price move another's
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
[`DemandSpec`](@ref) curve, which describes the composite. So

```julia
DemandSpec(region = "EU", p0 = 260, q0 = 110, elasticity = 0.5,
           product = DemandNest(sigma = 1.6, products = ["sawn_hw"], shares = [0.25],
                       nests = [DemandNest(sigma = 5.0, share = 0.75,
                                           products = ["sawn_sw", "panel"],
                                           shares   = [0.65, 0.35])]))
```

says that softwood lumber and panels are close substitutes, that hardwood
lumber is a more distant substitute for either, and that the construction
bundle as a whole has an own-price elasticity of 0.5.

The cross-price effects this implies are the CES ones: within a group, the
compensated cross-price elasticity of product ``i`` with respect to the price of
``j`` is ``\sigma s_j``, with ``s_j`` the value share of ``j``. A large
supplier therefore transmits its price more strongly than a small one at the
same ``\sigma`` — the shares matter as much as the elasticity, exactly as they
do between origins.

### Substitutes outside the model

Wood fuel competes with gas and oil, whose prices this model does not
determine. That needs none of the machinery above: an exogenous price is a
constant, so its effect on demand is a **shift of the whole curve**. Put the
cross-price elasticity in [`exogenous_shift`](@ref) as the driver elasticity:

```julia
DemandSpec(product = pellets_and_firewood, region = "EU", p0 = ..., q0 = ...,
           elasticity = 0.4,
           shift = exogenous_shift(gas = (level = 48.0, reference = 30.0,
                                          elasticity = 0.35)))
```

Gas 60% dearer than at calibration, with a cross-price elasticity of 0.35,
raises demand for the wood-fuel bundle by ``1.6^{0.35} \approx 1.18``. The
endogenous substitution inside the bundle and the exogenous substitution
against gas are different mechanisms, and it is worth keeping them visibly
apart.

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
    The program treats the revenue as a cost consumed, when it is a transfer
    received. The objective of a taxed model is therefore **welfare minus the
    revenue**, and comparing it with an untaxed one overstates the loss by the
    whole amount raised. `Results.taxes` reports what each instrument raises
    precisely so the transfer can be added back; who ultimately gets it is
    outside the model.

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
the cheapest. The standard answer, [Armington (1969)](@ref armington1969), is
to treat the varieties
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
Armington(product = "sawn_sw", destination = "EU", sigma = 2.5,
          nests = [OriginNest(sigma = 12, origins = ["EU", "NA"])])
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
`"all"` entry as the fallback, and the elasticity, the groups and the shares may
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

**Processing location does not respond smoothly.** Processes have constant
returns to scale and a constant value-added cost, so a process has no supply
response of its own: below a threshold set by freight it does not move at all,
and above it the whole industry moves at once.
[Why processes, and not supply and demand per transformation](@ref) measures
this and explains where it bites — chiefly on products left homogeneous, since
an [`Armington`](@ref) specification on the output restores a smooth response.

Otherwise the remedy is to say what capacity exists. [`Process`](@ref)'s
`capacity` bounds the activity level, region by region if need be, and a
binding bound earns a quasi-rent reported in `Results.activity` — what another
unit of mill would be worth there. Without such a bound, and without Armington
on its output, read the geography of *processing* far more cautiously than that
of prices or of final demand.

**Prices come from a continuous relaxation.** All quantities are continuous, so
the model knows nothing of indivisible investments, minimum efficient scale or
shutdown decisions.

## Calibration

A reasonable first reaction to all of this is: if I already know what was
produced last year, and at what price, what is the model for?

The answer is that an observation is a **point**, and a point with an
elasticity attached is a **curve**. The model never tells you the base year —
you measured it. It tells you where the curve's intersection moves when
something changes, and it does so for every product, region and process at
once, with material balance, zero profit and no-arbitrage all holding
simultaneously. That joint consistency is the thing no observation can give
you, because it is a property of the system rather than of any one market.

### A consistent benchmark is a fixed point

There is a sharp way to see what calibration does. If the data you anchor on
*is* an equilibrium, the model reproduces it exactly and adds nothing:

```
A  lumber   price 163.8222 → 163.8222    relative difference 2.6e-15
B  lumber   price 171.8222 → 171.8222    relative difference 1.5e-12
A  logs     price  56.9111 →  56.9111    relative difference 1.1e-12
```

That is *benchmark replication*, and reproducing it is the standard test that
the data and the structure agree with each other.

Real data never is an equilibrium — it carries measurement error, aggregation
over products and periods, margins the model does not represent, and markets
that were not clearing. Anchor one price 5% away from its consistent value and
the model refuses to return it:

```
A  lumber   anchored at 172.013   model returns 166.62   gap 3.1%
A  logs     anchored at  56.911   model returns  58.31   gap 2.5%
```

It splits the difference, because your numbers cannot all be true at once. That
gap is not a nuisance: it is the model telling you so, and its size is a
measure of how far your data is from being a coherent description of a market.

### What is actually calibrated

Entering `p0` and `q0` is not calibration, it is data entry. Calibration is
choosing the parameters you do **not** observe so that the benchmark becomes an
equilibrium:

* `vacost` of each process, from the observed output price less the observed
  cost of its inputs — which makes zero profit hold exactly at base-year
  prices;
* `transport_costs` on each route, from the observed price gap between its ends on
  routes that were actually used;
* [`armington_shares`](@ref) for the trade shares, from the base-year bilateral
  matrix at base-year delivered prices — not from the raw observed shares, for
  the reason given under
  [Imperfect substitution between origins (Armington)](@ref);
* `shift` on the curves, to absorb what remains.

The example in [The forest sector: France in the world](@ref) does none of
this. Its reference points are plausible figures, not a jointly consistent
benchmark, which is exactly why its equilibrium prices sit some way from them.
It is an illustration of the machinery, not a calibrated model.

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
moves between periods moves through `shift` and `capacity`, computed from state
and from exogenous drivers. The discipline in one sentence:

> Use model output to update **physical state**, never **behavioural
> parameters**.

Updating a forest stock by `stock + growth − harvest` is an accounting
identity, and the harvest may perfectly well come from the model. Updating `q0`
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
the same engine sit under a forest model with annual steps, a capacity-planning
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
intended seam: the loop updates `shift` and `capacity` between periods and
never touches the fitted curves.

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

These are things the engine itself might grow. A time loop is not among them,
for the reason given under [Dynamics are somebody else's job](@ref).

* **Armington refinements**: a separate composite per user (final demand and
  each process) instead of one per region, as GTAP does for the
  domestic-versus-imported margin; and a CET nest on the export side, so that
  redirecting output between destinations is costly for producers too, which is
  the mirror image of the import-side nest built here.
* **Inputs with an exogenous price** inside a nest (e.g. a resin with a world
  price). A `SupplySpec` with a very high elasticity is an approximation
  available today.
* **Ad valorem taxes without iteration**, which would need the model restated
  as a complementarity problem rather than an optimisation.
* **Estimation** of the parameters from data.

## Where these ideas come from

The works behind the choices on this page, and the code they can be read in,
are collected on the [References](@ref) page: a theory list cited as
*Author (year)*, and an implementation list cited as *[C1]*, *[C2]* and so on,
pointing at repositories where the same concept is written out and runs.
