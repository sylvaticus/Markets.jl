# Changelog

All notable changes to Markets.jl. Versions follow [semantic versioning](https://semver.org),
with the caveat that while the major version is 0 the **minor** version carries
breaking changes.

## [Unreleased]

### Fixed
- A process restricted to a single region given as a bare name,
  `Process(regions = "B", ...)`, was read as the letters of that name and
  rejected by `validate`. It is now that one region.
- The scarcity `rent` in `Results.supply` was the price minus the *untaxed*
  marginal cost, read off the price *users* pay. Under a tax on production it
  reported the tax as rent, even with the ceiling slack; for an Armington
  product it reported the gap between the user and the producer price. It is
  now the multiplier of the ceiling: the producer price minus the taxed
  marginal cost, zero unless the ceiling binds.

### Changed
- The solver is chosen with the `solver` keyword of `solve_market`, which
  replaces `optimizer`.
- **Names are strings.** Every region, product, process and tax is named with
  a `String` — `"EU"` rather than `:EU` — in the data types, in the result
  tables and throughout the documentation. A `Symbol` is still accepted
  wherever a name is wanted and converted on the way in, so the two spellings
  may be mixed and existing scripts keep working; what is stored, compared and
  reported back is always the string. The `on` field of a `Tax`, the `"all"`
  selector sentinel and the `on_failure` option of `solve_market` follow the
  same rule.

  Symbols would make the model's dictionary lookups about 7 ns cheaper each.
  A model of the size this package is for spends two seconds in the solver, so
  even a hundred thousand lookups would save under a millisecond of it — while
  the price would be that `spec.region == "EU"` quietly answers `false` for
  anyone who reads a field back off their own data.
- `MarketData.transport` is now `MarketData.transport_costs`, which is what it
  always held: the cost per unit of quantity carried on a route.
- `MarketData.tradable` is **removed**. It never reached the solver — a product
  was tradable because it had a route, and the list only existed to be checked
  against those routes. Writing a route in `transport_costs` is now the whole
  declaration, so there are no longer two lists to keep in step. A product that
  may move at no cost takes a route priced `0.0`.

### Added
- `ces_calibration(; quantities, prices, sigma)` calibrates any CES group —
  substitutes in a `DemandNest`, inputs of a `ces` nest, origins — from the
  base-year quantities and prices of its members. It returns the shares at
  equal prices, the CES price index and the composite quantity, which are the
  `shares`, `p0` and `q0` a demand curve over a bundle needs, and it takes the
  limiting forms at `sigma = 0`, `1` and `Inf`. `armington_shares` is now
  `ces_calibration` applied to each destination, and the forest example uses
  it in place of its own `calibrate`.
- **Elasticities of substitution of 1 and below.** Every CES — `ces` input
  nests, `Armington`, `OriginNest` and `DemandNest` — accepts any `sigma ≥ 0`
  instead of only `sigma > 1`, and so does `armington_shares` (finite values).
  Below 1 the members are poor substitutes, or complements. `sigma = 1` is
  written as the Cobb–Douglas `Πᵢ (qᵢ/δᵢ)^δᵢ` and `sigma = 0` as the Leontief
  `minᵢ qᵢ/δᵢ` (one linear constraint per member, with no quantity floor), the
  exact limits of the calibrated CES there, so `ces(...; sigma = 0)` is a
  fixed-proportion requirement of several products. Below about `sigma = 0.01`
  the CES formula leaves floating-point range: use 0.
- `jump_model(d)` returns the JuMP model `solve_market` would solve, built but
  not solved, so `print(jump_model(d))` shows the program in full. It takes the
  same solver keywords as `solve_market`. Every constraint is now named
  (`balance[r,p]`, `trade_limit[name]`, …), as the variables already were.
- A `MarketData` displays as a summary of the economy: what is supplied,
  consumed, made, used and traded of each product, the processes as a chart
  ordered by stage, the routes and their costs, the Armington specifications,
  the trade limits and the taxes and subsidies.
- `solve_market(d; solver_options = ...)` passes options to the solver, as a
  `NamedTuple`, a `Dict` or a vector of pairs, keyed by the solver's own option
  names or by MOI attributes. They are set last, so they override the
  tolerances Markets sets for Ipopt. An Ipopt wrapped in
  `optimizer_with_attributes` is still recognised as Ipopt.
- `Results` reports how the solve went: `raw_status` (the solver's own
  message), `solve_time` in seconds, `iterations` and the number of `solves`,
  the last three summed over the repeated solves an ad valorem tariff needs.
- A minimal worked example — two regions, two products, one dairy, and a milk
  that is both drunk and turned into cheese — as the first of an **Examples**
  section, in the README and as the opening documentation page. The forest
  model follows it as the advanced one.
- A modelling-choices section, **Why processes, and not supply and demand per
  transformation**, on the paradigm this package takes: what fixed-coefficient
  transformations buy in mass-balance consistency and data, what constant
  returns cost in behavioural realism, a measurement of how sharply processing
  location switches and why an `Armington` output removes it, and what a
  decreasing-returns process would cost in modelling and in data if it were
  added. Five references were added with it.
- A **Calibration** section in the modelling documentation: what calibration
  buys you, that a consistent benchmark is a fixed point while inconsistent
  data are not, which parameters are actually calibrated, and why a multi-period
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
