# Changelog

All notable changes to Markets.jl. Versions follow [semantic versioning](https://semver.org),
with the caveat that while the major version is 0 the **minor** version carries
breaking changes.

## [Unreleased]

### Changed
- `MarketData.transport` is now `MarketData.transport_costs`, which is what it
  always held: the cost per unit of quantity carried on a route.
- `MarketData.tradable` is **removed**. It never reached the solver — a product
  was tradable because it had a route, and the list only existed to be checked
  against those routes. Writing a route in `transport_costs` is now the whole
  declaration, so there are no longer two lists to keep in step. A product that
  may move at no cost takes a route priced `0.0`.

### Added
- A minimal worked example — two regions, two products, one dairy, and a milk
  that is both drunk and turned into cheese — as the first of an **Examples**
  section, in the README and as the opening documentation page. The forest
  model follows it as the advanced one.
- A **Calibration** section in the modelling documentation: what calibration
  buys you, that a consistent benchmark is a fixed point while inconsistent
  data is not, which parameters are actually calibrated, and why a multi-period
  run must never re-anchor `p0`/`q0` on its own output.
- A **References** documentation page carrying two lists: the theory behind
  each concept, cited as *Author (year)*, and open-source implementations of
  the same concept, cited as *[C1]*, *[C2]*… Citations are placed at the
  concepts throughout the documentation, and every DOI and repository was
  checked against CrossRef or the GitHub API.
- `validate(::MarketData)`, run by `solve_market` before anything is built. It
  rejects names the economy does not declare, duplicate regions, products,
  process names and curves, routes that go nowhere, tax and trade-limit
  selectors that match nothing or that the tax base ignores, and processes that
  would produce without limit.
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

### Changed (robustness)
- The quantity floor is now set **per variable**, at a millionth of that
  variable's own reference magnitude, instead of one absolute floor for the
  whole economy. A single obscure product with a tiny reference quantity used
  to drag every other variable's floor down with it.
- The fixed-point loop for ad valorem trade taxes halves its step whenever a
  pass increases the error rather than reducing it. Successive substitution
  contracts quickly here — several-hundred-per-cent rates settle in three to
  seven passes — but it has no guarantee, and the damping costs nothing when it
  is not needed.
- `disposal(product = ..., cost = ...)` writes the process that gives a
  by-product somewhere to go, so that free (or paid) disposal no longer has to
  be hand-rolled.

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
