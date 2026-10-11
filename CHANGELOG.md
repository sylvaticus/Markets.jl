# Changelog

All notable changes to Markets.jl. Versions follow [semantic versioning](https://semver.org),
with the caveat that while the major version is 0 the **minor** version carries
breaking changes.

## [Unreleased]

### Added
- `validate(::MarketData)`, run by `solve_market` before anything is built. It
  rejects names the economy does not declare, duplicate regions, products,
  process names and curves, transport routes for products left out of
  `tradable`, tax and trade-limit selectors that match nothing or that the tax
  base ignores, and processes that would produce without limit.
- `Results` carries `status`, `objective` and `welfare` (the objective with tax
  revenue added back), and a `supply` table reporting primary supply with its
  capacity and the scarcity rent a binding one earns.
- Every variable in the JuMP model is named, so `res.model` can be read.

### Changed
- `solve_market` **throws** by default when the solver returns anything but an
  optimal or locally optimal solution; `on_failure = :warn` or `:ignore`
  restores the previous behaviour. The two `ALMOST_` statuses warn and return.
- Ipopt's feasibility tolerance is set relative to `qfloor`. Its absolute
  default could be coarser than the quantities of an economy written in small
  units, which made some infeasible models report as solved.
- `ces` and `Nest` throw `ArgumentError` rather than `AssertionError`, and
  `Nest` now validates its arguments instead of trusting them.
- `consumption` is filtered by `tol` like the other quantity tables.
- Rents are reported with the sign the solver gives them rather than through
  `abs`, and `for_region` includes the `supply` table.

### Fixed
- A per-unit tax on the demand for a group of substitutes was levied on the
  group's composite — an index number — instead of on each product's own
  quantity. The two agree only when the members cost the same, so the wedge
  charged and the revenue reported could differ. Ad valorem rates were, and
  remain, exact.
- `DemandNest(sigma = Inf)` built a model containing `NaN`. Infinitely
  substitutable groups are now pooled by addition, as origin groups already
  were, in demand nests and process input nests alike.

## [0.1.0] — unreleased

A complete rewrite. The equation-system formulation of 0.0.x is replaced by a
spatial partial equilibrium solved as net-social-surplus maximisation, with
Armington trade, substitution in final demand, exogenous drivers, bounds and a
tax framework. Nothing of the 0.0.x API survives.

## [0.0.2] — 2022-06-09

The last state of the original formulation, tagged after the fact so that its
documentation has a version of its own. Never registered.

## [0.0.1] — 2022-05-30

The only version in the General registry. A different model from the one this
package now contains.
