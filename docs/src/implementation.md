# Code implementation

This page describes the mathematical program built by [`solve_market`](@ref),
how results are extracted, and the layout of the code. The reasons behind the
formulation are on the [Modelling choices](@ref) page.

## Source layout

| File | Content |
|:-----|:--------|
| `src/Markets.jl`  | module definition, exports |
| `src/types.jl`    | data schema: [`DemandSpec`](@ref), [`SupplySpec`](@ref), [`Nest`](@ref), [`leontief`](@ref), [`ces`](@ref), [`Process`](@ref), [`Armington`](@ref), [`MarketData`](@ref) |
| `src/model.jl`    | [`solve_market`](@ref): builds and solves the JuMP model |
| `src/results.jl`  | [`Results`](@ref) and [`for_region`](@ref): extraction of the result tables |
| `examples/forest/` | the forest-products example economy (`example_data.jl`) and a driver printing all tables (`run_example.jl`) |
| `test/runtests.jl` | test suite |
| `docs/`           | this documentation (Documenter.jl) |

Every type of the data schema is declared with `Base.@kwdef`, so it is built
with keyword arguments and its optional fields have defaults; the field
documentation on the [API reference](@ref) page is generated from the struct
definitions with
[DocStringExtensions](https://github.com/JuliaDocs/DocStringExtensions.jl).
The two nest builders, [`leontief`](@ref) and [`ces`](@ref), take keyword
arguments too.

The model is written with [JuMP](https://jump.dev) and solved by default with
[Ipopt](https://github.com/coin-or/Ipopt), an interior-point solver for
nonlinear programs. Variables and constraints are created only for the
combinations present in the data, and stored in dictionaries keyed by symbols
(e.g. `D[(region, product)]`, `T[(product, from, to)]`), not in dense arrays.

## The mathematical program

### Sets and data

| Symbol | Meaning |
|:-------|:--------|
| ``r, r' \in R`` | regions |
| ``p \in P`` | products |
| ``k \in K`` | processes; ``R_k \subseteq R`` the regions where ``k`` operates |
| ``n \in N_k`` | input nests of process ``k``; ``I_n`` the products of nest ``n`` |
| ``\mathcal{D}, \mathcal{S} \subseteq R \times P`` | (region, product) pairs with a demand / supply curve |
| ``\mathcal{T}`` | trade routes ``(p, r, r')`` present in `transport` |
| ``P^A \subseteq P`` | products with a finite Armington elasticity; ``P^H = P \setminus P^A`` the homogeneous ones |
| ``O_{p,r} \subseteq R`` | origins destination ``r`` buys ``p \in P^A`` from: itself, plus every origin with a route into ``r`` and a positive share |
| ``p_0, q_0, \eta, \varepsilon`` | reference price, reference quantity, demand and supply elasticities |
| ``y_{k,p}`` | yield of output ``p`` per unit of activity of ``k`` |
| ``\bar a_n, \delta_{n,i}, \sigma_n, \phi_n`` | composite requirement, shares, elasticity of substitution and scale of nest ``n`` |
| ``\sigma_p, \delta_{p,o,r}`` | Armington elasticity of ``p \in P^A`` and the value share of origin ``o`` in destination ``r`` (normalised over ``O_{p,r}``) |
| ``c_k`` | value-added cost per unit of activity of ``k`` (`vacost`) |
| ``\tau_{p,r,r'}`` | unit transport cost on a route |

### Variables

| Variable | Domain | Meaning |
|:---------|:-------|:--------|
| ``D_{r,p} \ge q_{min}`` | ``(r,p) \in \mathcal{D}`` | final demand |
| ``S_{r,p} \ge q_{min}`` | ``(r,p) \in \mathcal{S}`` | primary supply |
| ``z_{r,k} \ge 0`` | ``k \in K, r \in R_k`` | process activity |
| ``x_{r,k,n,i} \ge q_{min}`` | CES nests ``n``, ``i \in I_n`` | quantity of input ``i`` used in nest ``n`` |
| ``T_{p,r,r'} \ge 0`` | ``(p,r,r') \in \mathcal{T}``, ``p \in P^H`` | trade flow from ``r`` to ``r'`` |
| ``X_{p,o,r} \ge q_{min}`` | ``p \in P^A``, ``o \in O_{p,r}`` | quantity of the variety of origin ``o`` used in ``r`` |
| ``A_{p,r} \ge q_{min}`` | ``p \in P^A``, ``r \in R`` | composite of ``p`` available to the users of ``r`` |

``q_{min} = 10^{-4}`` is the constant `QFLOOR` (see
[Numerical details](@ref)).

### Objective

```math
\begin{aligned}
\max\; & \sum_{(r,p) \in \mathcal{D}} \frac{a_{r,p}}{1 - 1/\eta_{r,p}}\, D_{r,p}^{\,1 - 1/\eta_{r,p}}
       \;-\; \sum_{(r,p) \in \mathcal{S}} \frac{b_{r,p}}{1 + 1/\varepsilon_{r,p}}\, S_{r,p}^{\,1 + 1/\varepsilon_{r,p}} \\
       & \;-\; \sum_{k} \sum_{r \in R_k} c_k\, z_{r,k}
       \;-\; \sum_{(p,r,r') \in \mathcal{T}} \tau_{p,r,r'}\, T_{p,r,r'}
       \;-\; \sum_{p \in P^A} \sum_{r} \sum_{o \in O_{p,r},\, o \ne r} \tau_{p,o,r}\, X_{p,o,r}
\end{aligned}
```

with ``a = p_0\, q_0^{1/\eta}`` and ``b = p_0\, q_0^{-1/\varepsilon}``. The
first two sums are the closed-form integrals of the inverse demand curve
``P^D(D) = a D^{-1/\eta}`` and of the inverse supply curve
``P^S(S) = b S^{1/\varepsilon}``.

### CES nest constraints

For each process ``k``, region ``r \in R_k`` and nest ``n`` with more than one
product, with ``\rho_n = (\sigma_n - 1)/\sigma_n``:

```math
\phi_n \left( \sum_{i \in I_n} \delta_{n,i}^{\,1/\sigma_n}\, x_{r,k,n,i}^{\,\rho_n} \right)^{1/\rho_n} \;\ge\; \bar a_n\, z_{r,k}
```

The constraint is an inequality, which keeps the feasible set convex; at the
optimum it is binding because inputs are costly. A nest with a single product
(Leontief) creates no variable and no constraint: the input use
``\bar a_n z_{r,k}`` enters the material balance directly.

### Material balance of a homogeneous product

For every region ``r \in R`` and product ``p \in P^H``:

```math
\begin{aligned}
& S_{r,p} + \sum_{k} y_{k,p}\, z_{r,k} + \sum_{r'} T_{p,r',r} \\
=\; & D_{r,p} + \sum_{k} \sum_{\substack{n \in N_k \\ I_n = \{p\}}} \bar a_n\, z_{r,k}
      + \sum_{k} \sum_{\substack{n \in N_k \\ p \in I_n,\ |I_n| > 1}} x_{r,k,n,p}
      + \sum_{r'} T_{p,r,r'}
\end{aligned}
```

Terms are present only where defined (e.g. ``S_{r,p}`` only if
``(r,p) \in \mathcal{S}``). In the code the constraint is written as
`sources + imports − uses − exports == 0`.

### Balances of an Armington product

For ``p \in P^A`` the single balance splits in two, because what a region
produces and what its users consume are no longer the same good. Writing
``\text{src}_{r,p}`` and ``\text{use}_{r,p}`` for the source and use sides
above (without the trade terms), for every ``r \in R``:

```math
\text{src}_{r,p} = \sum_{r' :\, r \in O_{p,r'}} X_{p,r,r'}
\qquad\text{(origin balance: all local output is shipped somewhere)}
```

```math
A_{p,r} = \text{use}_{r,p}
\qquad\text{(absorption balance: the composite covers all local use)}
```

```math
\left( \sum_{o \in O_{p,r}} \delta_{p,o,r}^{\,1/\sigma_p}\, X_{p,o,r}^{\,\rho_p} \right)^{1/\rho_p} \;\ge\; A_{p,r},
\qquad \rho_p = \frac{\sigma_p - 1}{\sigma_p}
\qquad\text{(CES over origins)}
```

Note that ``o = r`` is in ``O_{p,r}``: a region's own variety goes through the
composite like any other, and ``X_{p,r,r}`` is its domestic use. A destination
with a single origin gets the linear constraint ``X_{p,o,r} \ge A_{p,r}``
instead, avoiding the powers.

The aggregator is the same calibrated CES as the input nests, so the implied
composite price is the standard price index

```math
\pi^c_{p,r} = \left( \sum_{o \in O_{p,r}} \delta_{p,o,r}\,
              \bigl(\pi^s_{p,o} + \tau_{p,o,r}\bigr)^{1-\sigma_p} \right)^{1/(1-\sigma_p)}
```

which tends to ``\min_o (\pi^s_{p,o} + \tau_{p,o,r})`` as ``\sigma_p \to
\infty``. At ``\sigma_p = \infty`` the engine does not build these constraints
at all: the product is treated as homogeneous, which is the same model.

### Convexity

The first term of the objective is concave (exponent ``1 - 1/\eta`` in
``(0,1)`` because ``\eta > 1``), the second is convex and subtracted, and the
others are linear, so the objective is concave. The CES aggregator is concave
for ``\rho < 1``, so the set where it is above a linear function is convex. The
material balances are linear. The problem is therefore a convex program, and
the local optimum found by Ipopt is the global one.

## From the solution to prices

The balance faced by the users of a region is stored in `balance[(r, p)]`, and
the price those users pay is the dual of that constraint:

```julia
abs(dual(balance[(r, p)]))        # the `price` column
```

The dual is the change in the objective per extra unit of product available to
them, i.e. the competitive price. In the forest example solved with Ipopt the
duals are already positive; `abs` guards against differences in sign
conventions between solvers.

For an Armington product the origin balance is stored as well, in
`origin_balance[(r, p)]`, and its dual is what local producers are paid:

```julia
abs(dual(origin_balance[(r, p)])) # the `producer_price` column
```

The two prices differ by the composition of the CES bundle: local users pay the
price index of all the varieties they buy, local producers are paid the value
of their own. For a homogeneous product there is no origin balance and
`producer_price` repeats `price`.

!!! note
    Because there is no free disposal, the economic price of a product in
    excess supply (e.g. a by-product with no profitable use) can be negative.
    `abs` would report such a price as positive. This does not happen in the
    forest example.

## Results extraction

`build_results` (internal) reads the variable values after the solve and fills
a [`Results`](@ref) with the tables:

* `production`: ``S_{r,p}`` plus ``\sum_k y_{k,p} z_{r,k}``, for all region and
  product pairs with a value above ``10^{-6}``;
* `consumption`: ``D_{r,p}`` for every pair in ``\mathcal{D}``;
* `trade`: ``T_{p,r,r'}`` and the cross-border ``X_{p,o,r}`` (``o \ne r``) above ``10^{-6}``;
* `net_trade`: exports and imports summed over routes, for pairs with any trade;
* `prices`: one row for every region and product, if the solver returned duals,
  with `price` and `producer_price` as above;
* `activity`: ``z_{r,k}`` above ``10^{-6}``.

Intermediate input use (the ``x`` variables and the Leontief uses) and the
domestic Armington varieties ``X_{p,r,r}`` are not reported as tables; they can
be recovered from `res.model`.

## Numerical details

* **Lower bound `QFLOOR`.** Demand, supply and CES input quantities are bounded
  below by ``10^{-4}`` instead of 0. The derivatives of ``D^{1-1/\eta}`` and of
  ``x^{\rho}`` are infinite at zero, which the interior-point solver cannot
  handle. As a consequence every demand and every input of an active CES nest
  takes at least ``10^{-4}`` units.
* **Starting point.** ``D`` and ``S`` start at their reference quantities
  ``q_0``, ``z`` at 1, CES inputs at ``\bar a_n``, and trade flows at 0.
* **Solver status.** If the termination status is neither `OPTIMAL` nor
  `LOCALLY_SOLVED`, [`solve_market`](@ref) emits a warning and still returns the
  results, which may then be unreliable.
* **Other solvers.** Any JuMP solver that supports nonlinear constraints and
  duals can be passed with the `optimizer` keyword.

## Tests

`test/runtests.jl` checks the solution against conditions that hold at the
equilibrium:

* a one-product, two-region market: the exporting region is the cheap one, the
  price gap equals the transport cost, prices equal the inverse demand and
  supply curves, and regional balances hold;
* a two-stage chain (a CES process feeding a Leontief process): the Leontief
  output price equals coefficient × input price + value-added cost, the CES
  input ratio satisfies the cost-minimisation condition, and the CES output price
  equals the composite unit cost + value-added cost;
* Armington trade: ``\sigma = \infty`` reproduces the homogeneous solution
  exactly, the error against it falls monotonically as ``\sigma`` grows, a
  finite ``\sigma`` produces cross-hauling and producer-price gaps wider than
  the transport cost, the composite price equals the CES price index of the
  delivered prices, the variety mix follows the CES demand condition, a supply
  shock in one region moves prices and flows in the other, and the
  specification errors are caught;
* the forest example: the solve succeeds, all prices are positive, no route has
  a price gap above its transport cost, used routes have a gap equal to it, and
  world exports equal world imports for every product.

Run them with:

```julia
using Pkg; Pkg.test("Markets")
```

## Building this documentation

From the package folder:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The pages are written to `docs/build/`. On GitHub, the CI workflow builds them
and deploys them to the `gh-pages` branch.
