# # The forest sector: France in the world
#
# This example builds a complete spatial partial-equilibrium model of the
# forest-products sector with [Markets.jl](https://github.com/sylvaticus/Markets.jl)
# and solves it. It is meant to be read top to bottom: every number the model
# needs is introduced where it is used, and the result tables are shown as they
# come out of the solver.
#
# The economy has **eight regions**. France is split into four, so that the
# model can say something about French regional markets, and the rest of the
# world into four blocs:
#
# | | |
# |:---|:---|
# | `SEF` | South-east France — the Alps and the Massif Central, spruce and fir |
# | `SWF` | South-west France — the Landes, Europe's largest maritime-pine plantation |
# | `GEF` | Grand Est — the Vosges, and the oak and beech of the north-east |
# | `OF`  | the rest of France |
# | `EU`  | the rest of the European Union |
# | `NA`  | North America |
# | `AS`  | Asia-Pacific |
# | `RW`  | the rest of the world |
#
# What makes the regional detail work is the way the varieties of a product
# substitute for each other, which is the subject of
# [the Armington section](@ref "Who substitutes for whom") below: French
# products are **perfect substitutes for each other**, close substitutes for
# the rest of the EU, and more distant substitutes for the other blocs.
#
# All the figures are illustrative — plausible orders of magnitude, not
# calibrated statistics. Quantities are in million m³ for wood products and
# million tonnes for pulp and paper; prices are in USD per m³ or per tonne.

using Markets, DataFrames

# ## Regions and products

france  = [:SEF, :SWF, :GEF, :OF]
regions = [france..., :EU, :NA, :AS, :RW]
nothing # hide

# The product list follows the FAOSTAT / UNECE chain, splitting roundwood by
# **coniferous** (softwood, `swr`) and **non-coniferous** (hardwood, `hwr`),
# because the two lead to different products.
#
# `smallwood` is the third assortment the forest yields: small-diameter
# roundwood from thinnings and tops, too poor for sawing. In French statistics
# it is the *bois d'industrie et bois énergie* category, and the name says what
# matters here — the same log can be burnt in a stove or sent to a mill, and
# which it is depends on what the two uses are willing to pay.

products = [:swr, :hwr, :smallwood, :chips, :pulp,
            :sawn_sw, :sawn_hw, :panel, :paper, :pellets]
nothing # hide

# They are connected by the transformations set up further down:
#
# ```
#   FOREST ── swr ───────┬─ sawmill_sw ─→ sawn_sw + chips
#          ── hwr ───────┼─ sawmill_hw ─→ sawn_hw + chips
#          ── smallwood ─┤
#                        ├─ panelmill  : CES(chips, smallwood, swr, hwr) ─→ panel
#                        ├─ pulpmill   : CES(swr, hwr, chips, smallwood) ─→ pulp ─→ papermill ─→ paper
#                        └─ pelletmill : CES(chips, smallwood, swr, hwr) ─→ pellets
#
#                           smallwood ─────────────────────────────────→ burnt as firewood
# ```
#
# Two things in that picture are worth more than a passing glance.
#
# `chips` are a sawmill by-product *and* an input to three other mills, so the
# competition for residues is part of the equilibrium rather than an assumption
# about it.
#
# `smallwood` is **both a final product and an intermediate one**. It comes out
# of the forest, and from there it either goes into a household stove as
# firewood — final demand, like sawnwood or paper — or into a mill as a raw
# material. The model does not need to be told which: both uses face the same
# local price, and the split between them is whatever that price makes optimal.
# A product in this engine may be primary-supplied, manufactured, consumed and
# used as an input all at once; its balance simply sums the sources and the
# uses.

# ## Primary supply
#
# Roundwood coming out of the forest, as a constant-elasticity curve through a
# reference point. The French regions harvest some 27 Mm³ of industrial
# roundwood between them: the Landes (`SWF`) dominate softwood, while `GEF`
# and `OF` hold most of the hardwood. Their supply elasticities are lower than
# the blocs' — French forests are fragmented between many small private owners,
# which makes the harvest slow to respond to price.

supply = SupplySpec[
    ## softwood roundwood
    SupplySpec(product = :swr, region = :SEF, p0 = 72, q0 =   3.2, elasticity = 0.40),
    SupplySpec(product = :swr, region = :SWF, p0 = 65, q0 =   7.5, elasticity = 0.45),
    SupplySpec(product = :swr, region = :GEF, p0 = 68, q0 =   2.6, elasticity = 0.40),
    SupplySpec(product = :swr, region = :OF,  p0 = 70, q0 =   3.4, elasticity = 0.40),
    SupplySpec(product = :swr, region = :EU,  p0 = 68, q0 = 270.0, elasticity = 0.60),
    SupplySpec(product = :swr, region = :NA,  p0 = 62, q0 = 380.0, elasticity = 0.60),
    SupplySpec(product = :swr, region = :AS,  p0 = 80, q0 = 150.0, elasticity = 0.55),
    SupplySpec(product = :swr, region = :RW,  p0 = 58, q0 = 190.0, elasticity = 0.55),
    ## hardwood roundwood
    SupplySpec(product = :hwr, region = :SEF, p0 = 78, q0 =   1.4, elasticity = 0.35),
    SupplySpec(product = :hwr, region = :SWF, p0 = 75, q0 =   1.2, elasticity = 0.35),
    SupplySpec(product = :hwr, region = :GEF, p0 = 80, q0 =   3.4, elasticity = 0.35),
    SupplySpec(product = :hwr, region = :OF,  p0 = 78, q0 =   4.6, elasticity = 0.35),
    SupplySpec(product = :hwr, region = :EU,  p0 = 76, q0 =  85.0, elasticity = 0.50),
    SupplySpec(product = :hwr, region = :NA,  p0 = 72, q0 =  70.0, elasticity = 0.50),
    SupplySpec(product = :hwr, region = :AS,  p0 = 74, q0 = 180.0, elasticity = 0.50),
    SupplySpec(product = :hwr, region = :RW,  p0 = 66, q0 = 130.0, elasticity = 0.50),
    ## smallwood: thinnings and tops, worth a fraction of a sawlog
    SupplySpec(product = :smallwood, region = :SEF, p0 = 42, q0 =   6.0, elasticity = 0.50),
    SupplySpec(product = :smallwood, region = :SWF, p0 = 38, q0 =   3.5, elasticity = 0.55),
    SupplySpec(product = :smallwood, region = :GEF, p0 = 40, q0 =   6.5, elasticity = 0.50),
    SupplySpec(product = :smallwood, region = :OF,  p0 = 41, q0 =   9.0, elasticity = 0.50),
    SupplySpec(product = :smallwood, region = :EU,  p0 = 39, q0 = 150.0, elasticity = 0.60),
    SupplySpec(product = :smallwood, region = :NA,  p0 = 34, q0 =  70.0, elasticity = 0.60),
    SupplySpec(product = :smallwood, region = :AS,  p0 = 45, q0 = 200.0, elasticity = 0.55),
    SupplySpec(product = :smallwood, region = :RW,  p0 = 30, q0 = 250.0, elasticity = 0.55),
]
nothing # hide

# ## Final demand
#
# Only finished products are consumed: construction lumber, appearance-grade
# hardwood lumber, panels and paper. The reference quantities are roughly
# proportional to population and income, so `OF` — which holds two thirds of
# the French population — is much the largest French market.
#
# The own-price elasticities are **below 1**: demand for forest products is
# usually found to be price-inelastic, construction lumber and paper most of
# all, since their cost is a small part of the building or the printed product
# they end up in. An inelastic demand makes prices, rather than quantities,
# absorb a supply shock — which matters for the storm scenario at the end.

reference_demand = Dict(
    ## product  => (price, quantities by region, elasticity)
    :sawn_sw => (250.0,  Dict(:SEF => 1.90, :SWF => 0.90, :GEF => 0.80, :OF => 5.60,
                              :EU => 75.0, :NA => 120.0, :AS => 80.0,  :RW => 45.0), 0.45),
    :sawn_hw => (300.0,  Dict(:SEF => 0.25, :SWF => 0.12, :GEF => 0.11, :OF => 0.72,
                              :EU => 12.0, :NA =>  20.0, :AS => 60.0,  :RW => 15.0), 0.55),
    :panel   => (350.0,  Dict(:SEF => 0.90, :SWF => 0.40, :GEF => 0.40, :OF => 2.80,
                              :EU => 45.0, :NA =>  40.0, :AS => 90.0,  :RW => 25.0), 0.50),
    :paper   => (1000.0, Dict(:SEF => 1.80, :SWF => 0.80, :GEF => 0.75, :OF => 5.60,
                              :EU => 70.0, :NA =>  80.0, :AS => 130.0, :RW => 40.0), 0.30),
    ## firewood: smallwood bought and burnt as it comes, the oldest use of all
    :smallwood => (55.0, Dict(:SEF => 3.20, :SWF => 1.50, :GEF => 2.60, :OF => 7.00,
                              :EU => 90.0, :NA =>  35.0, :AS => 150.0, :RW => 190.0), 0.35),
    ## pellets: the same energy, dried, ground and compressed
    :pellets   => (250.0, Dict(:SEF => 0.50, :SWF => 0.25, :GEF => 0.30, :OF => 1.40,
                              :EU => 22.0, :NA =>   9.0, :AS =>  12.0, :RW =>   4.0), 0.60),
)

demand = [DemandSpec(product = p, region = r, p0 = price, q0 = quantities[r], elasticity = η)
          for (p, (price, quantities, η)) in reference_demand for r in regions]
nothing # hide

# ## Processes
#
# Five transformations. Sawmilling is fixed-proportion joint production: a log
# yields lumber *and* chips. Panel and pulp mills instead draw on a CES bundle
# of wood inputs, so when chips become scarce — because the other mill is
# bidding for them — each mill shifts smoothly toward roundwood.
#
# Pulp mills are capital-heavy and few: in France only the Landes (`SWF`) and a
# handful of sites elsewhere (`OF`) have the capacity, which the `regions`
# argument expresses. `SEF` and `GEF` therefore have no pulp of their own and
# must buy it from somewhere else.

pulp_regions = [:SWF, :OF, :EU, :NA, :AS, :RW]

processes = [
    Process(name    = :sawmill_sw,
            inputs  = [leontief(product = :swr, coeff = 1.0)],
            outputs = [:sawn_sw => 0.50, :chips => 0.35],
            vacost  = 40),
    Process(name    = :sawmill_hw,
            inputs  = [leontief(product = :hwr, coeff = 1.0)],
            outputs = [:sawn_hw => 0.45, :chips => 0.35],
            vacost  = 45),
    Process(name    = :panelmill,
            inputs  = [ces(composite = 1.30,
                           products  = [:chips, :smallwood, :swr, :hwr],
                           shares    = [0.40, 0.25, 0.23, 0.12],
                           sigma     = 2.5)],
            outputs = [:panel => 1.0],
            vacost  = 120),
    Process(name    = :pulpmill,
            regions = pulp_regions,
            inputs  = [ces(composite = 4.0,
                           products  = [:swr, :hwr, :chips, :smallwood],
                           shares    = [0.35, 0.22, 0.23, 0.20],
                           sigma     = 1.8)],
            outputs = [:pulp => 1.0],
            vacost  = 300),
    Process(name    = :papermill,
            inputs  = [leontief(product = :pulp, coeff = 1.10)],
            outputs = [:paper => 1.0],
            vacost  = 200),
    ## the pellet mill bids for the same wood as the panel and pulp mills, and
    ## for the same smallwood that households would otherwise burn as it comes
    Process(name    = :pelletmill,
            inputs  = [ces(composite = 2.2,            # m³ of wood per tonne
                           products  = [:chips, :smallwood, :swr, :hwr],
                           shares    = [0.40, 0.40, 0.12, 0.08],
                           sigma     = 3.0)],
            outputs = [:pellets => 1.0],
            vacost  = 90),
]
nothing # hide

# ## Distance and transport costs
#
# Distances are relative, with the Europe–North America crossing as the unit.
# The French regions are an order of magnitude closer to each other than to
# anywhere else, which is what will keep most French wood inside France.

dist = Dict(
    (:SEF, :SWF) => 0.12, (:SEF, :GEF) => 0.10, (:SEF, :OF) => 0.08,
    (:SWF, :GEF) => 0.14, (:SWF, :OF)  => 0.08, (:GEF, :OF) => 0.07,
    (:SEF, :EU)  => 0.22, (:SWF, :EU)  => 0.26, (:GEF, :EU) => 0.18, (:OF, :EU) => 0.22,
    (:SEF, :NA)  => 1.10, (:SWF, :NA)  => 1.00, (:GEF, :NA) => 1.12, (:OF, :NA) => 1.05,
    (:SEF, :AS)  => 1.28, (:SWF, :AS)  => 1.35, (:GEF, :AS) => 1.30, (:OF, :AS) => 1.32,
    (:SEF, :RW)  => 0.82, (:SWF, :RW)  => 0.85, (:GEF, :RW) => 0.88, (:OF, :RW) => 0.85,
    (:EU, :NA)   => 1.00, (:EU, :AS)   => 1.20, (:EU, :RW)  => 0.80,
    (:NA, :AS)   => 1.10, (:NA, :RW)   => 1.00, (:AS, :RW)  => 0.95,
)
distance(a, b) = a == b ? 0.0 : get(dist, (a, b), get(dist, (b, a), 0.0))
nothing # hide

# Freight is charged per unit of quantity, not of value, because shipping a m³
# of logs costs much the same whether the logs are cheap or dear. Bulky,
# low-value products therefore carry a high transport cost relative to their
# price, which is what keeps roundwood trade regional and lets paper travel.

freight = Dict(:swr => 12.0, :hwr => 12.0, :smallwood => 13.0, :chips => 14.0,
               :pulp => 35.0, :sawn_sw => 30.0, :sawn_hw => 32.0, :panel => 35.0,
               :paper => 55.0, :pellets => 28.0)

tradable  = products
transport = Dict((p, a, b) => freight[p] * distance(a, b)
                 for p in tradable, a in regions, b in regions if a != b)
length(transport)

# ## Who substitutes for whom
#
# This is where the regional detail earns its place. If the varieties of a
# product were perfect substitutes everywhere, each market would simply buy
# from the cheapest source and a small cost advantage would take everything.
# If instead they were all imperfect substitutes with one elasticity, French
# oak would be no closer to French pine than to Asian hardwood, which is just
# as wrong.
#
# What the data says is a **hierarchy**:
#
# * inside France the product of one region is the *same good* as that of
#   another — a m³ of Landes pine and a m³ of Vosges pine of the same grade are
#   interchangeable, and only transport separates them;
# * between France and the rest of the EU the products are close but not
#   identical: shared standards and certification, short distances, but
#   different species mixes and supply contracts;
# * between the big blocs they are more distant substitutes: different species,
#   grading rules, phytosanitary requirements and commercial habits.
#
# That hierarchy is a tree, which is exactly what nested
# [`OriginNest`](@ref)s express: France pooled by `sigma = Inf`, then the
# French pool and the rest of the EU together at a high elasticity, then Europe
# as a whole against the other blocs at a lower one.

σ_blocs = Dict(:swr => 6.0, :hwr => 6.0, :smallwood => 7.0, :chips => 7.0,
               :pulp => 5.0, :sawn_sw => 3.5, :sawn_hw => 3.0, :panel => 3.0,
               :paper => 4.0, :pellets => 6.0)

# Commodities graded to a standard (chips, roundwood, pulp) substitute readily
# between blocs; manufactured products, where origin carries a reputation, do
# so less. Within Europe every product substitutes three times more easily than
# between blocs:

σ_europe = Dict(p => 3σ for (p, σ) in σ_blocs)
nothing # hide

europe(p) = OriginNest(sigma   = σ_europe[p],
                       origins = [:EU],
                       nests   = [OriginNest(sigma = Inf, origins = france)])
nothing # hide

# Note the direction of the inequality the package enforces: a group must hold
# together *at least as tightly* as it holds to the outside, `Inf ≥ 3σ ≥ σ`.
# The reverse would say that members of a group are worse substitutes for each
# other than for outsiders, which no cost function can represent.
#
# The shares are the other half of the specification: they say how much of each
# origin a market would take **if all delivered prices were equal**, so they
# carry the home bias and the trade pattern that prices alone do not explain.
# Calibrating them properly means fitting a base-year trade matrix; here they
# are generated by a gravity rule — big suppliers and near ones weigh more,
# and a region leans toward itself.

size_of = Dict(r => sum(s.q0 for s in supply if s.region == r) for r in regions)
decay, home_bias = 1.4, 3.0
shares = Dict((o, r) => size_of[o] * exp(-decay * distance(o, r)) * (o == r ? home_bias : 1.0)
              for o in regions, r in regions)
nothing # hide

# Shares inside the French group do not matter, as perfect substitutes are
# simply added up; what the group contributes to the level above is their sum.
#
# One [`Armington`](@ref) specification per product then assembles the three
# levels. The same structure applies to every destination here, but it need not
# — a `destination` argument would let the EU market treat a bloc differently
# from how that bloc treats the EU, which is what non-tariff measures do.

armington = [Armington(product = p,
                       sigma   = σ_blocs[p],
                       shares  = shares,
                       nests   = [europe(p)])
             for p in products]
nothing # hide

# ## Solving

example_market = MarketData(; regions, products, demand, supply,
                              processes, tradable, transport, armington)

res = solve_market(example_market)
nothing # hide

# ## Prices
#
# Because the products are differentiated by origin, two prices coexist in
# every region: what local users pay for the composite they consume, and what
# local producers are paid for their own variety.

unstack(res.prices, :region, :product, :price)

#-

unstack(res.prices, :region, :product, :producer_price)

# The `missing` entries are not a failure: `SEF` and `GEF` have no pulp mill,
# so they have no pulp of their own to price. They buy it, and the first table
# gives what they pay.
#
# Notice how close the four French rows are to each other in both tables, and
# how far they are from `NA` or `AS`. That is the nest at work.

# ## Quantities

unstack(res.production, :region, :product, :quantity)

#-

unstack(res.consumption, :region, :product, :quantity)

# ## Firewood or raw material?
#
# `smallwood` is the product that has to choose. Every region's harvest of it
# is split between the stove and the mills by nothing more than the local
# price, and the split comes out very differently from place to place:

quantity_of(table, r, p) = (rows = table[(table.region .== r) .& (table.product .== p),
                                            :quantity];
                            isempty(rows) ? 0.0 : only(rows))
shipped(direction, r) = sum(res.trade[(res.trade.product .== :smallwood) .&
                                      (res.trade[!, direction] .== r), :quantity]; init = 0.0)

smallwood_use = DataFrame(
    (region      = r,
     harvest     = quantity_of(res.production, r, :smallwood),
     firewood    = quantity_of(res.consumption, r, :smallwood),
     net_exports = shipped(:from, r) - shipped(:to, r))
    for r in regions)
smallwood_use.to_the_mills =
    smallwood_use.harvest .- smallwood_use.firewood .- smallwood_use.net_exports
smallwood_use

# The shares are worth reading off. North America burns under a quarter of what
# it cuts and sends the rest to industry; Asia-Pacific and the rest of the world
# burn over two fifths. The four French regions burn between a quarter (`GEF`,
# `SWF`) and a half (`OF`, where most of the population is) — and **export most
# of what is left rather than milling it**, because French industry is small
# next to its neighbours' and the wood is worth more across the border.
#
# None of that was imposed. The same wood is offered to the stove and to the
# mill at one price, and the quantities follow from which use outbids the other
# where.
#
# The pellet mills are bidding for that same wood, which is what makes the
# competition interesting — a pellet is wood that has been dried and pressed
# rather than burnt as it came, so the two uses are close substitutes in energy
# terms but reach the consumer as different products:

res.activity[res.activity.process .== :pelletmill, :]

# France produces little of its own pellets here even though it has the wood,
# because French wood is dearer than EU or RW wood and the mills have no
# capacity bound to hold them in place. That is the sharpness that constant
# returns to scale always produces in processing location — see
# [What the model does not constrain](@ref) — and it is the part of this
# example to treat with the most caution.

# ## Trade
#
# The ten largest flows:

first(sort(res.trade, :quantity, rev = true), 10)

# France as a whole, by product — what the four regions together buy from and
# sell to the rest of the world:

french_trade = combine(groupby(
    res.trade[in.(res.trade.from, Ref(france)) .⊻ in.(res.trade.to, Ref(france)), :],
    [:product]),
    :quantity => sum => :crossing_the_border)

# ## Is France really one market?
#
# The French regions were declared perfect substitutes for each other, so they
# should behave between themselves exactly as the homogeneous model of
# Samuelson does. That means two things, and both are worth checking rather
# than assuming.
#
# First, **no cross-hauling**: two regions never ship the same product to each
# other at once, because for an identical good that would be pure waste of
# freight. Indeed, the only flows inside France are of pulp, and they run one
# way:

res.trade[in.(res.trade.from, Ref(france)) .& in.(res.trade.to, Ref(france)), :]

# `SEF` and `GEF` have no pulp mill, so they buy pulp from `OF`, which has one.
# Nothing else moves between French regions: each covers its own needs first
# and sends what is left abroad, since carrying an identical good from one
# French region to another costs freight that the sale cannot repay.
#
# Second, **prices differ by no more than transport**, the no-arbitrage
# condition of a homogeneous market. Across every product and every pair of
# French regions, the largest gap still sits below the cost of shipping:

pp(r, p) = only(res.prices[(res.prices.region .== r) .& (res.prices.product .== p),
                           :producer_price])
gaps = DataFrame((product = p, from = a, to = b,
                  price_gap      = abs(pp(a, p) - pp(b, p)),
                  transport_cost = freight[p] * distance(a, b))
                 for p in products for a in france for b in france
                 if a < b && !ismissing(pp(a, p)) && !ismissing(pp(b, p)))
gaps.slack = gaps.transport_cost .- gaps.price_gap
(pairs = nrow(gaps), all_arbitrage_free = all(gaps.slack .>= -1e-6))

# The pairs where it binds most tightly — the ones where shipping would almost
# have been worth it:

first(sort(gaps, :slack), 6)

# ## A storm in the Landes
#
# In January 2009 the storm Klaus felled some 40 Mm³ of maritime pine in
# south-west France, several years of harvest at once. The salvage has to be
# sold, so the region's supply curve shifts far out.
#
# That is what the `shift` field of a [`SupplySpec`](@ref) is for. Everything
# about the curve other than its price response lives there — the standing
# volume, the road network, the labour and machines available — while `p0` and
# `q0` go on holding the calibration. Tripling `q0` would have worked here too,
# but it would have destroyed the reference point the curve was fitted to;
# shifting leaves it intact, which is what lets a scenario be run, reported and
# reversed. For a driver with a known elasticity,
# [`exogenous_shift`](@ref) turns levels into the multiplier:
#
# ```julia
# shift = exogenous_shift(growing_stock = (level = 118.0, reference = 100.0,
#                                          elasticity = 0.6))
# ```
#
# Here the salvage simply triples what the Landes offer at any price, and
# everything else is left alone:

reconfigure(d; kwargs...) =
    MarketData(; regions = d.regions, products = d.products, demand = d.demand,
                 supply = d.supply, processes = d.processes, tradable = d.tradable,
                 transport = d.transport, armington = d.armington, kwargs...)
nothing # hide

salvage = [s.region == :SWF && s.product == :swr ?
               SupplySpec(product = :swr, region = :SWF, p0 = s.p0, q0 = s.q0,
                          elasticity = s.elasticity, shift = 3.0) : s
           for s in supply]
nothing # hide

storm = solve_market(reconfigure(example_market, supply = salvage))
nothing # hide

# The question a spatial model exists to answer is how far the shock travels.
# Comparing the softwood price before and after, region by region:

price_of(r, p, result) = only(result.prices[(result.prices.region .== r) .&
                                            (result.prices.product .== p), :price])
nothing # hide

DataFrame((region = r,
           before = price_of(r, :swr, res),
           after  = price_of(r, :swr, storm),
           change_pct = 100 * (price_of(r, :swr, storm) / price_of(r, :swr, res) - 1))
          for r in regions)

# The gradient is the whole point. The salvage lands hardest on France, where
# the varieties are the same good and the extra wood simply adds to one pool;
# it passes with some loss into the rest of the EU, a close substitute; and it
# barely reaches North America or Asia, which buy a different good. Had every
# product been homogeneous worldwide, a single world price would have moved
# everywhere at once; had every variety been equally distant, the shock would
# have stayed bottled up in `SWF`. Neither is what the forest sector looks
# like.

# ## Where to go from here
#
# * Change the elasticities in `σ_blocs` and `σ_europe` and watch the gradient
#   above widen or flatten.
# * Replace the gravity shares with a base-year trade matrix. The `shares` are
#   the shares *at equal delivered prices*, not the observed ones, so pass the
#   matrix through [`armington_shares`](@ref) rather than using it directly.
# * Give the EU market its own `destination` specification to represent a
#   non-tariff measure applying to one bloc only.
# * Make transport concave in distance, or compute distances between the
#   centres of gravity of the forest resource and of the population rather than
#   between region centroids.
# * Add capacity bounds to the mills, and link periods through forest growth,
#   for a recursive-dynamic run.

#src ---------------------------------------------------------------------------
#src When run as a script rather than rendered as a page, print the main tables.
if abspath(PROGRAM_FILE) == @__FILE__                                        #src
    section(t) = (println("\n", "="^70); println("  ", t); println("="^70))   #src
    section("PRICES paid by users");   show(unstack(res.prices, :region, :product, :price); allcols = true)  #src
    section("PRICES got by producers"); show(unstack(res.prices, :region, :product, :producer_price); allcols = true) #src
    section("PRODUCTION");  show(unstack(res.production, :region, :product, :quantity); allcols = true)      #src
    section("CONSUMPTION"); show(unstack(res.consumption, :region, :product, :quantity); allcols = true)     #src
    section("LARGEST TRADE FLOWS"); show(first(sort(res.trade, :quantity, rev = true), 15))                  #src
    println()                                                                 #src
end                                                                           #src
