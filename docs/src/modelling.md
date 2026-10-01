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

Demand elasticities must be **greater than 1**. With ``\eta \le 1`` the area
under the demand curve from zero, ``\int_0^D P^D``, is infinite, and the
gross consumer benefit in the objective would be undefined. Supply elasticities
only need to be positive.

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

The material balance is an equality: every unit of a by-product must be used,
traded or consumed. There is no free disposal.

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
same product to each other (no cross-hauling), and trade patterns
can be more extreme than observed ones. Making products from different
origins imperfect substitutes (Armington, 1969) is a planned extension; see
[Planned extensions](@ref).

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
proxy for an elasticity, so its results depend on the path and the units, and
its errors compound across periods. Anchoring each period on a static
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

* **Armington trade**: imperfect substitution between origins, through a second
  CES nest over origins in each region's consumption.
* **Recursive dynamics**: a time loop around [`solve_market`](@ref), as
  described above.
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
