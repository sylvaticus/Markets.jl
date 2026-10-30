# Equilibrium, curves and welfare

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
and given the programming formulation everything since has used, by
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
  price and the supply curve (the rent of the resource owners, e.g. mine or
  land owners);
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
matters because that is how demand for most fuels, foods and basic materials is
estimated — for forest products, for example, see
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
term is still increasing and concave and the program is still convex. What is
lost is the *level* of welfare, and with it some of the comparisons between
scenarios: see [Welfare when demand is inelastic](@ref).

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

### Straight lines

A curve can instead be a straight line, ``D = \alpha + \beta p`` with
``\beta < 0`` or ``S = \alpha + \beta p`` with ``\beta > 0``: the form of
the classic spatial-equilibrium models from
[Samuelson (1952)](https://www.jstor.org/stable/1810137) and Takayama and
Judge onwards, and of most mathematical-programming sector models. Its integral
is quadratic, so the objective stays concave and the program convex.

The two forms agree at the reference point and in the slope there, so small
changes are the same in both. They part on large ones. A line has a finite
*choke price*, ``-\alpha/\beta``, above which nothing is demanded; a
constant-elasticity demand never reaches zero. And the elasticity of a line
grows as the price rises: a demand that is inelastic at the reference point
becomes elastic further up. A large supply shock therefore raises the price
less along a line than along a constant-elasticity curve calibrated to the
same point. Neither is more correct; the line is what a linear regression
estimates, and what a model built on one should be ported with.

A line through a reference point is the tangent there to the
constant-elasticity curve, `form = :linear` in [`DemandCurve`](@ref) and
[`SupplyCurve`](@ref). The flat curve, `elasticity = Inf`, is the same in both
forms: a price given from outside, as world prices are for a small country.

## Everything other than the price

A curve says how much is offered or wanted *at a given price*. Its position
depends on everything else: how much ore is in the ground or how many fish
are in the sea, how good the roads and ports are, how many workers and machines
there are, how many people live in the region and what they earn. Those are exogenous to a partial-equilibrium
model and have to enter it somehow. There are three standard devices, and they
are not interchangeable.

### A multiplicative shifter — the workhorse

Scale the whole curve by a factor built from the drivers and the elasticity of
the curve with respect to each:

```math
S(P) = \underbrace{\prod_j \left(\frac{Z_j}{Z_j^0}\right)^{\gamma_j}}_{\text{shift}}
       \; q_0 \left(\frac{P}{p_0}\right)^{\varepsilon}
```

This is what [`SupplyCurve`](@ref)'s and [`DemandCurve`](@ref)'s `shift` field
holds, and what [`exogenous_shift`](@ref) builds. It is the device most
sector partial-equilibrium models use to link a resource to what is taken from
it — the forest to the harvest, the fish stock to the catch — and there are
three good reasons for that.

It **comes free with the estimation**. The regression that gives the price
elasticity — ``\log S`` on ``\log P`` — gives the ``\gamma_j`` in the same
step if the drivers are in it. No separate exercise is needed.

It is **safe**: a positive factor times a positive curve stays positive, at any
level of any driver. An additive shifter, ``S = q_0 (P/p_0)^\varepsilon + A``,
can drive supply negative and destroys the constant-elasticity reading of
``\varepsilon``. Use addition only where something genuinely arrives in a
fixed quantity regardless of price — a stockpile released by a government,
scrap collected under a recycling mandate, a quota of imports from outside the
model. In this
package that is a separate [`SupplyCurve`](@ref) held to that quantity by a
[`Limit`](@ref), not a shifter.

It **separates calibration from scenario**. `p0` and `q0` keep holding the
base-year point the curve was fitted to; the scenario lives in `shift`. Rewrite
`q0` instead and the calibration is gone, which matters most in a
recursive-dynamic run, where the time loop updates the shifter from the state of the
resource each period and must not touch the fitted curve.

!!! note "Quantity shifters and cost shifters are the same thing here"
    A new road to a mine does not create ore, it lowers the cost of getting
    it out — a *cost* shifter, not a quantity one. With an isoelastic
    curve the distinction collapses: multiplying the quantity by ``m`` is
    exactly multiplying the marginal cost by ``m^{-1/\varepsilon}``, since
    ``b = p_0 (m q_0)^{-1/\varepsilon}``. With ``\varepsilon = 0.5``, offering
    10% more ore at every price is the same as mining it 17% more cheaply at
    every volume. What matters is only that the elasticity you apply is the one
    your ``\gamma`` was estimated against.

    On a straight line the two part. Multiplying the quantity rotates the line
    about the price at which supply starts, ``-\alpha/\beta``, which stays
    where it was; making every unit cheaper would lower that price too. A
    driver that works through cost belongs in ``\alpha`` on a line, which is
    also where a regression in levels puts it.

### A ceiling — for a limit that binds

A shifter scales a curve; it never stops it. When the constraint is a genuine
limit — a fishing quota, a mining licence, the capacity of the rigs, boats
and crews that exist this year — the right device is a bound:

```math
S \le \bar S(Z)
```

which is a [`Limit`](@ref) on the supply:

```julia
@limit "cod quota" supply(region = "EU", product = "cod") <= 350
```

The difference is not cosmetic. At
a binding ceiling the price parts company with the marginal cost on the curve,
and the gap is the **scarcity rent** of the resource — the economically
interesting quantity, and the one a pure shifter cannot produce. Doubling a
shifter doubles supply at every price; a ceiling says that beyond a point no
price will buy more.

The rent is the price producers get minus what the curve says the last unit
cost, taxes on production included:

```julia
rent = producer_price - ((1 + t_rel) * Markets.curve_price(s, S) + t_abs)
```

which is the multiplier of the ceiling, and is how `Results.supply` and
`Results.limits` report it: zero whenever the ceiling does not bind.

A ceiling set above what would have been produced anyway changes nothing, so it
is safe to carry one on every curve that has a known physical limit.

A ceiling can also be **soft**: crossed at a given cost per unit, which caps
its rent at that cost. That is a device for holding a quantity near a
reference — trade near its recent levels, in the tradition of the recursive
forest sector models — rather than a physical limit, and its cost is not part
of welfare. See [Soft limits](@ref).

### A structural supply — when the factors are the point

The two devices above are reduced-form: the curve is estimated, and the drivers
move it. The alternative is to derive supply from the technology instead, by
making the resource a product and its extraction a process. In a fishery
model, for example:

* the fish stock becomes a product with its own curve or endowment, bounded by
  a [`Limit`](@ref);
* fishing becomes a [`Process`](@ref) that consumes the stock, fuel and
  whatever else the model prices, and yields landed fish;
* labour, boats and access then enter as priced inputs or as bounds, and
  the *shape* of the resulting supply curve is derived rather than assumed.

This is how a computable general equilibrium model treats factors of
production, and it is the honest route when the question is about the factors
themselves — the effect of a labour shortage, of a fuel subsidy, of a fleet
buy-back. The price to pay is the data: factor shares and factor prices by
region, which a reduced-form shifter does not need.

The [forest example](@ref "The forest sector: France in the world") takes the
reduced-form route instead: its roundwood comes from supply curves, and its
storm scenario moves one of them through its `shift` rather than modelling the
forest itself.

### Choosing

| Driver | Device |
|:---|:---|
| Resource stock: reserves, fish stock, farmland | shifter, ``\gamma`` typically well below 1 |
| Quota, licence, legal cap | ceiling |
| Roads, ports, accessibility | shifter (a cost shifter in disguise) |
| Workers, rigs and machines | ceiling, or structural if that is the question |
| Income, population | shifter on demand |
| Stockpile release, recycled material | a separate inelastic curve with a ceiling |

The sensible default for a sector model is the second column of devices, with
ceilings wherever a physical limit is known, and the structural route reserved
for the questions that are really about factors.

## Welfare when demand is inelastic

`Results.welfare` is net social surplus at the solution: the terms of the
objective evaluated with every tax left out, a tax being a transfer, and
without the penalties of [soft limits](@ref "Soft limits"). The `objective`
itself is not welfare as soon as there is a tax or a soft limit. Whether the
level of welfare means anything, and which differences between scenarios do,
depends on the demand curves.

**When every demand curve has ``\eta > 1``** the antiderivative in the
objective is zero at zero, so it is the area under the curve from zero.
Welfare is then a true level, and two scenarios can be compared whatever
differs between them. A straight line or a flat curve counts as such a curve:
the area under it from zero is finite.

**When at least one demand curve has ``\eta \le 1``** the area under it from
zero is infinite: as the quantity falls towards zero the price rises fast
enough that the area never stops growing. The model keeps the antiderivative
of [Demand and supply curves](@ref) and leaves out the infinite remainder, so

```math
W_{\text{true}} = K + W_{\text{reported}}
```

where ``K`` is infinite and depends on the curve: its elasticity, its reference
point and its shifter. Three things follow.

### The level means nothing

Neither the size nor the sign of `welfare` or `objective` says anything on its
own. The [forest example](@ref "The forest sector: France in the world"),
whose demand elasticities lie between 0.3 and 0.45, reports a negative
objective, and that is not a market destroying value.

### Differences are exact when the inelastic curves are the same

Between two runs whose inelastic demand curves are identical, ``K`` is the
same in both and cancels:

```math
(K + W_2) - (K + W_1) = W_2 - W_1
```

The difference is then the exact change in consumer surplus, producer surplus
and tax revenue together, in money: its size as well as its sign. When an
inelastic curve itself differs, the two runs leave out different infinite
amounts, the difference between their `welfare` is meaningless, and a
different but equally valid antiderivative could even reverse its sign:

| What differs between the two runs | Difference in `welfare` | Prices, quantities, trade, activity, rents |
|:---|:---|:---|
| supply curves, their shifters included; value-added and transport costs; limits; taxes; demand curves with ``\eta > 1`` | exact, in size and sign | exact |
| on a demand curve with ``\eta \le 1``: its elasticity, its shifter, its reference point (`p0`, `q0`) or its form | meaningless, in size and in sign | exact |

The last column holds in every case. Prices and quantities come from the
first-order conditions of the program, which ``K`` does not enter, so a
projection that moves the demand curves with their shifters is exact in
everything it reports except welfare.

### Only differences, not ratios

A difference removes ``K``; a ratio keeps it. The relative change
``(W_2 - W_1) / (K + W_1)`` has a numerator that is right and a denominator
that depends on where the zero was put. A fall in reported welfare from −1000
to −1020 is a loss of 20 whatever ``K`` is, but as a percentage:

| Where the zero is put | Before | After | Change |
|:---|---:|---:|---:|
| as reported, ``K = 0`` | −1000 | −1020 | +2%, a "rise", the base being negative |
| ``K = 5\,000`` | 4000 | 3980 | −0.5% |
| ``K = 1\,000\,000`` | 999,000 | 998,980 | −0.002% |

Every line describes the same loss, and none is more correct than another:
with the true, infinite ``K`` every change is 0%. It is the difference between
Celsius and Kelvin, the same degree with the zero elsewhere: warming from 10 °C
to 20 °C is +10 degrees on both scales, but +100% in Celsius and +3.5% in
Kelvin. To express a welfare change in relative terms, divide it by something
that has a true zero: the value of the market, ``\pi D``, producer surplus, or
the revenue a tax raises.

### When the inelastic curves differ

When the scenarios move an inelastic demand curve, the gap in consumer surplus
between the two curves is itself infinite — a demand shifted outwards lies
above the old one at every quantity — so no figure the model could report
would compare them. What remains comparable is:

* **producer surplus**, ``\pi S - \int_0^S P^S``, since the area under a supply
  curve is finite at every elasticity;
* **consumer surplus over a common range of prices**: the area under the
  demand curve between the price ``p`` consumers pay and a ceiling ``\bar p``
  chosen once for all the scenarios. For a curve
  ``D = \text{shift} \cdot q_0 (p/p_0)^{-\eta}`` on one product it is

  ```math
  \int_p^{\bar p} D = \text{shift} \cdot q_0\, p_0^{\eta}\,
      \frac{\bar p^{\,1-\eta} - p^{\,1-\eta}}{1-\eta} \quad (\eta \ne 1),
  \qquad \text{shift} \cdot q_0\, p_0 \log\frac{\bar p}{p} \quad (\eta = 1),
  ```

  finite for every ``\eta``. The ceiling is a convention, such as the price of
  a backstop substitute; the measure depends on it, so it has to be reported
  with the results.
