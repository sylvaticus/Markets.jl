# Code implementation

This page describes the mathematical program built by [`solve_market`](@ref),
how results are extracted, and the layout of the code. The reasons behind the
formulation are on the [Modelling choices](@ref) page.

## Source layout

| File | Content |
|:-----|:--------|
| `src/Markets.jl`  | module definition, exports |
| `src/types.jl`    | data schema: [`DemandSpec`](@ref), [`SupplySpec`](@ref), [`Nest`](@ref), [`leontief`](@ref), [`ces`](@ref), [`Process`](@ref), [`MarketData`](@ref) |
| `src/model.jl`    | [`solve_market`](@ref): builds and solves the JuMP model |
| `src/results.jl`  | [`Results`](@ref) and [`for_region`](@ref): extraction of the result tables |
| `examples/forest/` | the forest-products example economy (`example_data.jl`) and a driver printing all tables (`run_example.jl`) |
| `test/runtests.jl` | test suite |
| `docs/`           | this documentation (Documenter.jl) |

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
| ``p_0, q_0, \eta, \varepsilon`` | reference price, reference quantity, demand and supply elasticities |
| ``y_{k,p}`` | yield of output ``p`` per unit of activity of ``k`` |
| ``\bar a_n, \delta_{n,i}, \sigma_n, \phi_n`` | composite requirement, shares, elasticity of substitution and scale of nest ``n`` |
| ``c_k`` | value-added cost per unit of activity of ``k`` (`vacost`) |
| ``\tau_{p,r,r'}`` | unit transport cost on a route |

### Variables

| Variable | Domain | Meaning |
|:---------|:-------|:--------|
| ``D_{r,p} \ge q_{min}`` | ``(r,p) \in \mathcal{D}`` | final demand |
| ``S_{r,p} \ge q_{min}`` | ``(r,p) \in \mathcal{S}`` | primary supply |
| ``z_{r,k} \ge 0`` | ``k \in K, r \in R_k`` | process activity |
| ``x_{r,k,n,i} \ge q_{min}`` | CES nests ``n``, ``i \in I_n`` | quantity of input ``i`` used in nest ``n`` |
| ``T_{p,r,r'} \ge 0`` | ``(p,r,r') \in \mathcal{T}`` | trade flow from ``r`` to ``r'`` |

``q_{min} = 10^{-4}`` is the constant `QFLOOR` (see
[Numerical details](@ref)).

### Objective

```math
\begin{aligned}
\max\; & \sum_{(r,p) \in \mathcal{D}} \frac{a_{r,p}}{1 - 1/\eta_{r,p}}\, D_{r,p}^{\,1 - 1/\eta_{r,p}}
       \;-\; \sum_{(r,p) \in \mathcal{S}} \frac{b_{r,p}}{1 + 1/\varepsilon_{r,p}}\, S_{r,p}^{\,1 + 1/\varepsilon_{r,p}} \\
       & \;-\; \sum_{k} \sum_{r \in R_k} c_k\, z_{r,k}
       \;-\; \sum_{(p,r,r') \in \mathcal{T}} \tau_{p,r,r'}\, T_{p,r,r'}
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

### Material balance

For every region ``r \in R`` and product ``p \in P``:

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

### Convexity

The first term of the objective is concave (exponent ``1 - 1/\eta`` in
``(0,1)`` because ``\eta > 1``), the second is convex and subtracted, and the
others are linear, so the objective is concave. The CES aggregator is concave
for ``\rho < 1``, so the set where it is above a linear function is convex. The
material balances are linear. The problem is therefore a convex program, and
the local optimum found by Ipopt is the global one.

## From the solution to prices

Each material-balance constraint is stored in `balance[(r, p)]`, and the
price of product ``p`` in region ``r`` is read as the dual of that constraint:

```julia
abs(dual(balance[(r, p)]))
```

The dual is the change in the objective per extra unit of product available in
that region, i.e. the competitive price. In the forest example solved with
Ipopt the duals are already positive; `abs` guards against differences in
sign conventions between solvers.

!!! note
    Because there is no free disposal, the economic price of a product in
    excess supply (e.g. a by-product with no profitable use) can be negative.
    `abs` would report such a price as positive. This does not happen in the
    forest example.

## Results extraction

The [`Results`](@ref) constructor reads the variable values after the solve and
builds the tables:

* `production`: ``S_{r,p}`` plus ``\sum_k y_{k,p} z_{r,k}``, for all region and
  product pairs with a value above ``10^{-6}``;
* `consumption`: ``D_{r,p}`` for every pair in ``\mathcal{D}``;
* `trade`: ``T_{p,r,r'}`` above ``10^{-6}``;
* `net_trade`: exports and imports summed over routes, for pairs with any trade;
* `prices`: one row for every region and product, if the solver returned duals;
* `activity`: ``z_{r,k}`` above ``10^{-6}``.

Intermediate input use (the ``x`` variables and the Leontief uses) is not
reported as a table; it can be recovered from `res.model`.

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
