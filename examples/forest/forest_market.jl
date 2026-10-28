# # The forest sector: France in the world
#
# This example builds a complete spatial partial-equilibrium model of the
# forest-products sector with [Markets.jl](https://github.com/sylvaticus/Markets.jl),
# solves it, and runs three scenarios on it. It is the first of three pages:
#
# 1. **this one**, on how to write the model and use it, with the kind of
#    numbers that are at hand before any calibration effort;
# 2. [Calibrating with partial data](@ref "The forest sector: calibrating with partial data"),
#    on turning incomplete, inconsistent statistics into the same model;
# 3. [Complete data, and what calibration changes](@ref "The forest sector: complete data, and what calibration changes"),
#    on what perfect information would take, and how much the scenarios
#    depend on which of the three you have.
#
# !!! warning "For illustrative purposes only"
#     The data in this example are plausible, but they are not statistics,
#     and are not calibrated: what that means, and why it matters, is the
#     subject of the last section. For simplicity, the structure of the model
#     and its data are also written together here; the next pages read the
#     data from files, as an application should.
#
# The economy has **eight regions**. France is split into four regions and the
# rest of the world into four blocs:
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
# products are **near-perfect substitutes for each other**, close substitutes
# for the rest of the EU, and more distant substitutes for the other blocs.
#
# The sector modelled here is the one the Global Forest Products Model
# established as the reference, [Buongiorno et al. (2003)](@ref gfpm2003); for
# the same economics in working open-source code, with different choices about
# dynamics and data, see TiMBA, [[C3]](@ref C3).
#
# Quantities are in million m³ for wood products and million tonnes for pulp,
# paper and pellets; prices are in USD per m³ or per tonne.

using Markets, DataFrames

# ## Regions and products

france  = ["SEF", "SWF", "GEF", "OF"]
regions = [france..., "EU", "NA", "AS", "RW"]
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

products = ["swr", "hwr", "smallwood", "chips", "pulp",
            "sawn_sw", "sawn_hw", "panel", "plywood", "paper", "pellets"]
nothing # hide

# They are connected by the transformations set up further down:
#
# ```
#   FOREST ── swr ───────┬─ sawmill_sw ─→ sawn_sw + chips
#          ── hwr ───────┼─ sawmill_hw ─→ sawn_hw + chips
#          ── smallwood ─┤
#                        ├─ plymill    : CES(swr, hwr) ───────────────────→ plywood + chips
#                        ├─ panelmill  : CES(chips, smallwood, swr, hwr) ─→ panel
#                        ├─ pulpmill   : CES(swr, hwr, chips, smallwood) ─→ pulp ─→ papermill ─→ paper
#                        └─ pelletmill : CES(chips, smallwood, swr, hwr) ─→ pellets
#
#                           smallwood ─────────────────────────────────→ burnt as firewood
# ```
#
# Two things in that picture are worth more than a passing glance.
#
# `chips` are a sawmill and plymill by-product *and* an input to three other mills, so the
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
# reference point: each region's harvest at its roadside price, which forest
# statistics give. The French regions harvest some 31 Mm³ of industrial
# roundwood between them: the Landes (`SWF`) dominate softwood, while `GEF`
# and `OF` hold most of the hardwood. Their supply elasticities are lower than
# the blocs' — French forests are fragmented between many small private owners,
# which makes the harvest slow to respond to price.

supply = SupplySpec[
    ## softwood roundwood
    SupplySpec(product = "swr", region = "SEF", p0 =  95, q0 =   3.58, elasticity = 0.40),
    SupplySpec(product = "swr", region = "SWF", p0 =  93, q0 =   8.79, elasticity = 0.45),
    SupplySpec(product = "swr", region = "GEF", p0 =  96, q0 =   2.99, elasticity = 0.40),
    SupplySpec(product = "swr", region = "OF",  p0 =  96, q0 =   3.86, elasticity = 0.40),
    SupplySpec(product = "swr", region = "EU",  p0 =  98, q0 = 335.0,  elasticity = 0.60),
    SupplySpec(product = "swr", region = "NA",  p0 =  94, q0 = 488.0,  elasticity = 0.60),
    SupplySpec(product = "swr", region = "AS",  p0 = 135, q0 = 200.0,  elasticity = 0.55),
    SupplySpec(product = "swr", region = "RW",  p0 = 103, q0 = 260.0,  elasticity = 0.55),
    ## hardwood roundwood
    SupplySpec(product = "hwr", region = "SEF", p0 = 111, q0 =   1.58, elasticity = 0.35),
    SupplySpec(product = "hwr", region = "SWF", p0 = 112, q0 =   1.38, elasticity = 0.35),
    SupplySpec(product = "hwr", region = "GEF", p0 = 108, q0 =   3.78, elasticity = 0.35),
    SupplySpec(product = "hwr", region = "OF",  p0 = 108, q0 =   5.16, elasticity = 0.35),
    SupplySpec(product = "hwr", region = "EU",  p0 = 117, q0 = 105.0,  elasticity = 0.50),
    SupplySpec(product = "hwr", region = "NA",  p0 = 132, q0 =  94.7,  elasticity = 0.50),
    SupplySpec(product = "hwr", region = "AS",  p0 = 131, q0 = 240.0,  elasticity = 0.50),
    SupplySpec(product = "hwr", region = "RW",  p0 = 107, q0 = 166.0,  elasticity = 0.50),
    ## smallwood: thinnings and tops, worth less than a sawlog
    SupplySpec(product = "smallwood", region = "SEF", p0 = 61, q0 =   8.31, elasticity = 0.50),
    SupplySpec(product = "smallwood", region = "SWF", p0 = 63, q0 =   5.35, elasticity = 0.55),
    SupplySpec(product = "smallwood", region = "GEF", p0 = 62, q0 =   9.17, elasticity = 0.50),
    SupplySpec(product = "smallwood", region = "OF",  p0 = 62, q0 =  12.5,  elasticity = 0.50),
    SupplySpec(product = "smallwood", region = "EU",  p0 = 68, q0 = 244.0,  elasticity = 0.60),
    SupplySpec(product = "smallwood", region = "NA",  p0 = 80, q0 = 137.0,  elasticity = 0.60),
    SupplySpec(product = "smallwood", region = "AS",  p0 = 79, q0 = 317.0,  elasticity = 0.55),
    SupplySpec(product = "smallwood", region = "RW",  p0 = 64, q0 = 439.0,  elasticity = 0.55),
]
nothing # hide

# ## Final demand
#
# Only finished products are consumed. Each region's consumption is known, and
# `OF` — which holds two thirds of the French population — is much the largest
# French market. For prices, what is usually at hand is **one price per
# product**: a world or national average, the same for every region.

reference = Dict(
    ## product    => (price, quantities by region)
    "sawn_sw"   => (205.0, Dict("SEF" => 1.95,  "SWF" => 0.931,  "GEF" => 0.815,  "OF" => 5.66,
                               "EU"  => 80.4,  "NA"  => 134.0,  "AS"  => 45.5,   "RW" => 40.3)),
    "sawn_hw"   => (280.0, Dict("SEF" => 0.228, "SWF" => 0.108,  "GEF" => 0.101,  "OF" => 0.658,
                               "EU"  => 10.7,  "NA"  => 15.4,   "AS"  => 53.3,   "RW" => 15.4)),
    "panel"     => (240.0, Dict("SEF" => 1.03,  "SWF" => 0.461,  "GEF" => 0.465,  "OF" => 3.29,
                               "EU"  => 48.9,  "NA"  => 43.0,   "AS"  => 110.0,  "RW" => 29.1)),
    ## plywood: veneers glued crosswise, dearer than particleboard and used
    ## where a panel has to carry a load
    "plywood"   => (285.0, Dict("SEF" => 0.196, "SWF" => 0.0952, "GEF" => 0.0811, "OF" => 0.562,
                               "EU"  => 10.2,  "NA"  => 19.8,   "AS"  => 93.2,   "RW" => 19.4)),
    "paper"     => (935.0, Dict("SEF" => 1.83,  "SWF" => 0.814,  "GEF" => 0.763,  "OF" => 5.71,
                               "EU"  => 71.1,  "NA"  => 81.4,   "AS"  => 132.0,  "RW" => 41.8)),
    ## firewood: smallwood bought and burnt as it comes, the oldest use of all
    "smallwood" => (70.0,  Dict("SEF" => 2.36,  "SWF" => 1.07,   "GEF" => 1.99,   "OF" => 5.05,
                               "EU"  => 59.3,  "NA"  => 17.1,   "AS"  => 97.6,   "RW" => 166.0)),
    ## pellets: the same energy, dried, ground and compressed
    "pellets"   => (285.0, Dict("SEF" => 0.637, "SWF" => 0.322,  "GEF" => 0.395,  "OF" => 1.73,
                               "EU"  => 27.5,  "NA"  => 12.3,   "AS"  => 20.3,   "RW" => 7.26)),
)
nothing # hide

# The own-price elasticities are **below 1**: demand for forest products is
# usually found to be price-inelastic, construction lumber and paper most of
# all, since their cost is a small part of the building or the printed product
# they end up in. An inelastic demand makes prices, rather than quantities,
# absorb a supply shock — which matters for the storm scenario at the end.
#
# ### Products that compete for the same job
#
# Those seven products are not seven separate markets. Softwood lumber and panels
# do the same job in a building and a builder swaps one for the other readily;
# hardwood lumber is a more distant substitute for either. Pellets and firewood
# are the same heat in different packaging. Paper competes with nothing else
# here.
#
# That is a [`DemandNest`](@ref) each: a tree of products with one demand curve
# on the bundle. It has to be a tree rather than a matrix of cross-price
# elasticities, for the reason set out in
# [Substitution in final demand](@ref) — a free matrix is not the gradient of
# any benefit function, and the model would stop being the market it claims to
# represent.
#
# ### Calibrating a bundle
#
# A demand curve over a bundle needs three numbers that no statistic reports:
# the **shares** of its members, and the **price** and **quantity** of the
# bundle itself, which are the `p0` and `q0` of its curve. What the `reference`
# table gives instead is a price and a quantity for each product.
# [`ces_calibration`](@ref) turns the second into the first, choosing the three
# so that at the base-year prices the model buys exactly the base-year
# quantities. For the wood fuel of `OF`, 5.05 Mm³ of firewood at 70 and 1.73
# Mt of pellets at 285:

ces_calibration(quantities = [5.05, 1.73], prices = [70.0, 285.0], sigma = 3)

# * The **shares** are what each product's share would be if both cost the
#   same. Pellets get 96% of the bundle although households buy three times
#   more firewood: they are bought little *because* they are four times
#   dearer, and would be bought far more at an equal price.
# * The **price** is the CES price index of the bundle — the cost of one unit of
#   it at the base-year prices.
# * The **quantity** is the base-year spending, 847, divided by that price. It
#   is an index number, not a volume: m³ of firewood and tonnes of pellets
#   would mean nothing added up, while their values at a common price do. That
#   is what a CES composite is.
#
# Feed those back into the CES demand condition at the base-year prices and
# out come exactly 5.05 and 1.73.
#
# ### Wood in buildings
#
# Construction is a nest inside a nest — softwood lumber, particleboard and
# plywood together at `σ = 5`, hardwood lumber against that group at `σ = 1.6`.
# The three in the inner group are what a builder actually chooses between for
# a floor or a wall; appearance-grade hardwood is a different kind of purchase.
#
# A nested bundle is calibrated from the inside out. The inner group is
# calibrated on its three products; its price and quantity then make it **one
# member** of the outer bundle, next to hardwood lumber; and the outer price
# and quantity are the reference point of the curve:

σ_lumber, σ_construction, σ_fuel = 5.0, 1.6, 3.0
lumber = ["sawn_sw", "panel", "plywood"]

function construction_demand(r)
    q(p) = reference[p][2][r]                    # base-year quantity in region r
    price(p) = reference[p][1]                   # base-year price
    inner = ces_calibration(quantities = q.(lumber), prices = price.(lumber),
                            sigma = σ_lumber)
    outer = ces_calibration(quantities = [q("sawn_hw"), inner.quantity],
                            prices     = [price("sawn_hw"), inner.price],
                            sigma = σ_construction)
    return DemandSpec(region = r, p0 = outer.price, q0 = outer.quantity, elasticity = 0.45,
        product = DemandNest(sigma = σ_construction,
            products = ["sawn_hw"], shares = [outer.shares[1]],
            nests = [DemandNest(sigma = σ_lumber, share = outer.shares[2],
                                products = lumber, shares = inner.shares)]))
end

# The `elasticity` of 0.45 is the own-price elasticity of the whole bundle:
# how much less wood goes into buildings when all of it gets dearer. How the
# bundle splits when only one product gets dearer is the job of the two `σ`.
#
# ### Wood for heat
#
# Wood fuel is a single nest, and the one place where a price from outside the
# model matters: firewood and pellets both compete with gas and oil, whose
# prices this sector does not determine. An exogenous price needs no nest at
# all — being a constant, it shifts the whole curve — so the cross-price
# elasticity goes into [`exogenous_shift`](@ref). Here gas sits at its
# reference level, which leaves the shifter at exactly 1; the scenario at the
# end moves it.

gas_price, gas_reference, gas_cross_elasticity = 30.0, 30.0, 0.35
fuels = ["smallwood", "pellets"]

function fuel_demand(r, gas = gas_price)
    fuel = ces_calibration(quantities = [reference[p][2][r] for p in fuels],
                           prices     = [reference[p][1] for p in fuels], sigma = σ_fuel)
    return DemandSpec(region = r, p0 = fuel.price, q0 = fuel.quantity, elasticity = 0.35,
        shift = exogenous_shift(gas = (level = gas, reference = gas_reference,
                                       elasticity = gas_cross_elasticity)),
        product = DemandNest(sigma = σ_fuel, products = fuels, shares = fuel.shares))
end

# ### All of final demand
#
# Paper competes with nothing here, so it keeps a plain curve through its own
# reference point. The whole demand side is then one curve per bundle and
# region. It is built by a function of the gas price rather than written out
# once, so that the gas scenario at the end can rebuild it at another price:

build_demand(gas = gas_price) =
    [[construction_demand(r) for r in regions];
     [fuel_demand(r, gas) for r in regions];
     [DemandSpec(product = "paper", region = r, p0 = reference["paper"][1],
                 q0 = reference["paper"][2][r], elasticity = 0.30) for r in regions]]

demand = build_demand()
nothing # hide

# ## Processes
#
# Six transformations, of three kinds.
#
# **Fixed proportions with a by-product.** Sawmilling turns a log into lumber
# *and* chips in set amounts. The plywood mill is the same shape — veneer is
# peeled from the log, and what is left, nearly half of it, comes out as chips
# — except that its input is a CES bundle of the two roundwood grades, because
# a veneer mill can work softwood or hardwood and will lean toward whichever is
# cheaper.
#
# **A CES bundle of wood.** Panel, pulp and pellet mills draw on a mix of
# chips, smallwood and roundwood, so when chips become scarce — because the
# other mills are bidding for them — each shifts smoothly toward roundwood.
#
# Plywood sits on both sides of that: it sheds chips, which helps the pulp
# mill, and bids for the same logs, which does not. Which effect wins is not
# something to guess at; it is the kind of link a sector model exists to
# settle.
#
# The recipes below — the shares of each wood in a bundle, the value-added
# costs — are reasonable guesses, the same in every region. The next page
# works them out from data instead.
#
# **Fixed proportions, single input.** The paper mill takes pulp and nothing
# else.
#
# Pulp mills are capital-heavy and few: in France only the Landes (`SWF`) and a
# handful of sites elsewhere (`OF`) have the capacity, which the `regions`
# argument expresses. `SEF` and `GEF` therefore have no pulp of their own and
# must buy it from somewhere else.

pulp_regions = ["SWF", "OF", "EU", "NA", "AS", "RW"]

# Pellet mills, by contrast, exist everywhere — but only so many of them. A
# [`Process`](@ref) with constant returns would pile into whichever region makes
# pellets most cheaply; the installed plant is what stops it, region by region.
# The plants' capacity, in Mt of pellets a year, becomes a limit on the mills'
# activity further down:

pellet_capacity = Dict("SEF" => 0.427, "SWF" => 0.50, "GEF" => 0.349, "OF" => 0.645,
                       "EU"  => 28.3,  "NA"  => 15.3, "AS"  => 19.9,  "RW" => 14.1)

processes = [
    Process(name    = "sawmill_sw",
            inputs  = [leontief(product = "swr", coeff = 1.0)],
            outputs = ["sawn_sw" => 0.50, "chips" => 0.35],
            vacost  = 40),
    Process(name    = "sawmill_hw",
            inputs  = [leontief(product = "hwr", coeff = 1.0)],
            outputs = ["sawn_hw" => 0.45, "chips" => 0.35],
            vacost  = 45),
    ## Plywood is peeled from veneer logs, so it takes roundwood and nothing
    ## else — and leaves a great deal of it behind as chips, which the panel,
    ## pulp and pellet mills then bid for.
    Process(name    = "plymill",
            inputs  = [ces(composite = 2.0,
                           products  = ["swr", "hwr"],
                           shares    = [0.45, 0.55],
                           sigma     = 1.8)],
            outputs = ["plywood" => 1.0, "chips" => 0.8],
            vacost  = 180),
    Process(name    = "panelmill",
            inputs  = [ces(composite = 1.30,
                           products  = ["chips", "smallwood", "swr", "hwr"],
                           shares    = [0.40, 0.25, 0.23, 0.12],
                           sigma     = 2.5)],
            outputs = ["panel" => 1.0],
            vacost  = 120),
    Process(name    = "pulpmill",
            regions = pulp_regions,
            inputs  = [ces(composite = 4.0,
                           products  = ["swr", "hwr", "chips", "smallwood"],
                           shares    = [0.35, 0.22, 0.23, 0.20],
                           sigma     = 1.8)],
            outputs = ["pulp" => 1.0],
            vacost  = 300),
    Process(name    = "papermill",
            inputs  = [leontief(product = "pulp", coeff = 1.10)],
            outputs = ["paper" => 1.0],
            vacost  = 200),
    ## the pellet mill bids for the same wood as the panel and pulp mills, and
    ## for the same smallwood that households would otherwise burn as it comes
    Process(name    = "pelletmill",
            inputs  = [ces(composite = 2.2,            # m³ of wood per tonne
                           products  = ["chips", "smallwood", "swr", "hwr"],
                           shares    = [0.40, 0.40, 0.12, 0.08],
                           sigma     = 3.0)],
            outputs = ["pellets" => 1.0],
            vacost  = 90),
]
nothing # hide

# ## Distance and transport costs
#
# Distances are relative, with the Europe–North America crossing as the unit.
# The French regions are an order of magnitude closer to each other than to
# anywhere else, which is what will keep most French wood inside France.

dist = Dict(
    ("SEF", "SWF") => 0.12, ("SEF", "GEF") => 0.10, ("SEF", "OF") => 0.08,
    ("SWF", "GEF") => 0.14, ("SWF", "OF")  => 0.08, ("GEF", "OF") => 0.07,
    ("SEF", "EU")  => 0.22, ("SWF", "EU")  => 0.26, ("GEF", "EU") => 0.18, ("OF", "EU") => 0.22,
    ("SEF", "NA")  => 1.10, ("SWF", "NA")  => 1.00, ("GEF", "NA") => 1.12, ("OF", "NA") => 1.05,
    ("SEF", "AS")  => 1.28, ("SWF", "AS")  => 1.35, ("GEF", "AS") => 1.30, ("OF", "AS") => 1.32,
    ("SEF", "RW")  => 0.82, ("SWF", "RW")  => 0.85, ("GEF", "RW") => 0.88, ("OF", "RW") => 0.85,
    ("EU", "NA")   => 1.00, ("EU", "AS")   => 1.20, ("EU", "RW")  => 0.80,
    ("NA", "AS")   => 1.10, ("NA", "RW")   => 1.00, ("AS", "RW")  => 0.95,
)
distance(a, b) = a == b ? 0.0 : get(dist, (a, b), get(dist, (b, a), 0.0))
nothing # hide

# Freight is charged per unit of quantity, not of value, because shipping a m³
# of logs costs much the same whether the logs are cheap or dear. Bulky,
# low-value products therefore carry a high transport cost relative to their
# price, which is what keeps roundwood trade regional and lets paper travel.

freight = Dict("swr" => 12.0, "hwr" => 12.0, "smallwood" => 13.0, "chips" => 14.0,
               "pulp" => 35.0, "sawn_sw" => 30.0, "sawn_hw" => 32.0, "panel" => 35.0,
               "plywood" => 38.0, "paper" => 55.0, "pellets" => 28.0)

transport_costs = Dict((p, a, b) => freight[p] * distance(a, b)
                       for p in products, a in regions, b in regions if a != b)
length(transport_costs)

# ## Who substitutes for whom
#
# This is where the regional detail earns its place. If the varieties of a
# product were perfect substitutes everywhere, each market would simply buy
# from the cheapest source and a small cost advantage would take everything.
# If instead they were all imperfect substitutes with one elasticity, French
# oak would be no closer to French pine than to Asian hardwood, which is just
# as wrong.
#
# What the data say is a **hierarchy**:
#
# * inside France the product of one region is nearly the *same good* as that
#   of another — a m³ of Landes pine and a m³ of Vosges pine of the same grade
#   are all but interchangeable, and mostly transport separates them;
# * between France and the rest of the EU the products are close but not
#   identical: shared standards and certification, short distances, but
#   different species mixes and supply contracts;
# * between the big blocs they are more distant substitutes: different species,
#   grading rules, phytosanitary requirements and commercial habits.
#
# That hierarchy is a tree, which is exactly what nested
# [`OriginNest`](@ref)s express: the French regions together at a very high
# elasticity, then France and the rest of the EU at a high one, then Europe as
# a whole against the other blocs at a lower one.

σ_blocs = Dict("swr" => 6.0, "hwr" => 6.0, "smallwood" => 7.0, "chips" => 7.0,
               "pulp" => 5.0, "sawn_sw" => 3.5, "sawn_hw" => 3.0, "panel" => 3.0,
               "plywood" => 2.5, "paper" => 4.0, "pellets" => 6.0)

# Commodities graded to a standard (chips, roundwood, pulp) substitute readily
# between blocs; manufactured products, where origin carries a reputation, do
# so less. Within Europe every product substitutes three times more easily than
# between blocs, and within France more easily still:

σ_europe = Dict(p => 3σ for (p, σ) in σ_blocs)
σ_france = 30.0

europe(p) = OriginNest(sigma   = σ_europe[p],
                       origins = ["EU"],
                       nests   = [OriginNest(sigma = σ_france, origins = france)])
nothing # hide

# Note the direction of the inequality the package enforces: a group must hold
# together *at least as tightly* as it holds to the outside, `30 ≥ 3σ ≥ σ`.
# The reverse would say that members of a group are worse substitutes for each
# other than for outsiders, which no cost function can represent.
#
# Why not make the French regions *perfect* substitutes, `sigma = Inf`, and
# France a single market up to transport? Because with mills that have
# constant returns, that can leave the model unable to say *where in France*
# the milling happens — see [Perfect substitutes and constant returns](@ref)
# for why, and for the remedies. A finite elasticity, however high, gives each
# region's product a demand of its own, which settles it.
#
# The shares are the other half of the specification: they say how much of each
# origin a market would take **if all delivered prices were equal**, so they
# carry the home bias and the trade pattern that prices alone do not explain.
# Calibrating them properly means fitting a base-year trade matrix, which is
# what the next page does; here they come from a gravity rule — big suppliers
# and near ones weigh more, and a region leans toward itself.

size_of = Dict(r => sum(s.q0 for s in supply if s.region == r) for r in regions)
decay, home_bias = 1.4, 3.0
shares = Dict((o, r) => size_of[o] * exp(-decay * distance(o, r)) * (o == r ? home_bias : 1.0)
              for o in regions, r in regions)
nothing # hide

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

# ## What cannot move, and what is taxed
#
# Three final things stand between the economy as described and the
# equilibrium.
#
# ### What the forest holds
#
# The supply curves say how much wood owners are willing to sell at each
# price. They say nothing of how much there is: a year's harvest cannot take
# more than the stands hold. That limit is not one per product, because the
# products share their stands. Softwood logs come from conifers, hardwood logs
# from broadleaves, and smallwood — thinnings and tops — from both, so a m³ of
# conifer can be sold as a log or as smallwood, but not as both.
#
# Here is what each region can cut in a year, in Mm³, by kind of stand —
# about a fifth more than it cuts in a normal year:

stock = Dict("SEF" => (conifers =  11.0, broadleaves =   6.0),
             "SWF" => (conifers =  16.0, broadleaves =   4.0),
             "GEF" => (conifers =   7.0, broadleaves =  13.0),
             "OF"  => (conifers =   9.0, broadleaves =  18.0),
             "EU"  => (conifers = 520.0, broadleaves = 340.0),
             "NA"  => (conifers = 700.0, broadleaves = 190.0),
             "AS"  => (conifers = 360.0, broadleaves = 570.0),
             "RW"  => (conifers = 480.0, broadleaves = 610.0))
stands_of = Dict("swr" => [:conifers], "hwr" => [:broadleaves],
                 "smallwood" => [:conifers, :broadleaves])
primary = ["swr", "hwr", "smallwood"]
nothing # hide

# The rule that keeps every m³ counted once is the one of the Forest Sector
# Model FFSM: for **every combination** of the primary products, their total
# harvest cannot exceed the stock of all the stands at least one of them can
# come from. Softwood logs alone are bounded by the conifers; softwood logs
# and smallwood together by conifers and broadleaves, since smallwood can come
# from either. Checking each product on its own would not do — two products
# may each fit while together drawing on the same stands beyond what they hold
# — and checking every combination is exactly enough: by Hall's theorem, the
# harvests then fit into the stands without any m³ being counted twice.
#
# With three products that is seven combinations per region. Each is a
# [`@limit`](@ref): a linear constraint on the quantities of the model, here on
# the sum of the supplies of a combination. Inside `@limit`, `supply(…)` names
# the supply of the model, not the curves of the same name above:

combinations(xs) = [xs[[isodd(i >> (j - 1)) for j in eachindex(xs)]]
                    for i in 1:(2^length(xs) - 1)]
available(r, C, stock) = sum(getfield(stock[r], s) for s in union((stands_of[p] for p in C)...))

forest_limits(stock) =
    [@limit("$r: " * join(C, " + "),
            sum(supply(region = r, product = p) for p in C) <= available(r, C, stock))
     for r in regions for C in combinations(primary)]

forest_limits(stock)[1:7]

# ### Plants, logistics and policy
#
# The pellet plants are the simplest kind of limit: one per region, on the
# activity of its pellet mills. A family of limits, one per region, is what
# [`@limits`](@ref) writes — `pellet_plant[SEF]`, `pellet_plant[SWF]`, and so on:

pellet_limits = @limits "pellet_plant" [r in regions] activity(region = r, process = "pelletmill") <=
                                                      pellet_capacity[r]
pellet_limits[1:2]

# Wood leaving the Landes has to be carried on the roads and through the ports
# that are there. One limit caps the **sum** of what `SWF` exports of the two
# products that travel by road, wherever it goes. At 18 Mm³ it is slack in a
# normal year — the region ships about 13 — and the storm scenario at the end
# is what makes it bite.

evacuation = @limit("landes_evacuation",
                    exports(region = "SWF", product = ["swr", "smallwood"]) <= 18.0)
nothing # hide

# All the limits of the economy, for a given stock of the stands — the
# scenarios below change the stock and nothing else:

all_limits(stock) = [forest_limits(stock); pellet_limits; evacuation]
nothing # hide

# And the rest of the world taxes its log exports, as several timber exporters
# do, to keep the wood at home for its own mills. The rate is **ad valorem**, a
# quarter of the price at the border, which is the one kind of wedge this model
# cannot write into its objective — the price is a dual of the problem. It is
# found by iteration instead, which [`solve_market`](@ref) does on its own; see
# [Why ad valorem rates on trade are special](@ref).

taxes = [Tax(name = "rw_log_duty", on = "trade", product = ["swr", "hwr"], from = "RW",
             rate = 0.25, relative = true)]

# ## Solving

example_market = MarketData(; regions, products, demand, supply, processes,
                              transport_costs, armington, taxes,
                              limits = all_limits(stock))

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

# ## What binds, and what it is worth
#
# Several of the things set up above are limits rather than curves, and the
# solution says what each is worth at the margin.
#
# The pellet mills first. Where the installed plant is enough, the capacity
# shows no rent; where it is not, the bound binds and the rent is what another
# tonne of annual capacity would earn its owner — the signal that would, in a
# model with investment, build the next mill. In the reference run every plant
# has room to spare:

res.activity[res.activity.process .== "pelletmill", :]

# The limits come next, the evacuation ceiling of the Landes and the stands of
# every region together. Their `rent` is what one more unit of the bound would
# be worth; all are slack in a normal year, so all are zero:

res.limits[res.limits.rent .!= 0, :]

# The evacuation ceiling, for one: the region ships about 13 Mm³ of roundwood
# against a ceiling of 18, so it costs nothing and earns nothing.

res.limits[res.limits.limit .== "landes_evacuation", :]

# And the export duty on logs leaving the rest of the world raises this, at a
# quarter of the border price on every m³ that still crosses:

res.taxes

# That revenue is a **transfer**, not a cost: the program subtracts it like any
# other cost, so the objective of this model is welfare minus the revenue
# rather than welfare. Adding it back is what makes a taxed run comparable with
# an untaxed one.

# ## Firewood or raw material?
#
# `smallwood` is the product that has to choose. Every region's harvest of it
# is split between the stove and the mills by nothing more than the local
# price, and the split comes out very differently from place to place:

quantity_of(table, r, p) = (rows = table[(table.region .== r) .& (table.product .== p),
                                            :quantity];
                            isempty(rows) ? 0.0 : only(rows))
shipped(direction, r) = sum(res.trade[(res.trade.product .== "smallwood") .&
                                      (res.trade[!, direction] .== r), :quantity]; init = 0.0)

smallwood_use = DataFrame(
    (region      = r,
     harvest     = quantity_of(res.production, r, "smallwood"),
     firewood    = quantity_of(res.consumption, r, "smallwood"),
     net_exports = shipped(:from, r) - shipped(:to, r))
    for r in regions)
smallwood_use.to_the_mills =
    smallwood_use.harvest .- smallwood_use.firewood .- smallwood_use.net_exports
smallwood_use

# The shares are worth reading off. North America burns a tenth of what it
# cuts, and imports wood for its mills on top; Asia-Pacific and the rest of the
# world burn between a quarter and two fifths. The French regions burn between
# a fifth (`SWF`, `GEF`) and two fifths (`OF`, where most of the population
# is) — and **export most of what is left rather than milling it**, because
# French industry is small next to its neighbours' and the wood is worth more
# across the border.
#
# None of that was imposed. The same wood is offered to the stove and to the
# mill at one price, and the quantities follow from which use outbids the other
# where.
#
# The pellet mills are bidding for that same wood, which is what makes the
# competition interesting — a pellet is wood that has been dried and pressed
# rather than burnt as it came, so the two uses are close substitutes in energy
# terms but reach the consumer as different products:

res.activity[res.activity.process .== "pelletmill", :]

# ## Trade
#
# The ten largest flows:

first(sort(res.trade, :quantity, rev = true), 10)

# France as a whole, by product — what the four regions together buy from and
# sell to the rest of the world:

french_trade = combine(groupby(
    res.trade[in.(res.trade.from, Ref(france)) .⊻ in.(res.trade.to, Ref(france)), :],
    [:product]),
    :quantity => sum => "crossing_the_border")

# ## How close is France to one market?
#
# The French regions were declared near-perfect substitutes, not perfect ones,
# and the difference is visible in two places.
#
# First, wood moves between French regions **in both directions** — the
# cross-hauling a homogeneous market would never do, because for an identical
# good it would be a pure waste of freight. Here it is a modest share of what
# France makes:

intra = res.trade[in.(res.trade.from, Ref(france)) .& in.(res.trade.to, Ref(france)), :]
made_in_france = combine(groupby(res.production[in.(res.production.region, Ref(france)), :],
                                 :product), :quantity => sum => :made_in_france)
leftjoin(combine(groupby(intra, :product), :quantity => sum => :between_regions),
         made_in_france, on = :product)

# Second, **producer prices can differ by more than the freight** between two
# French regions, which no-arbitrage would forbid for an identical good:

pp(r, p) = only(res.prices[(res.prices.region .== r) .& (res.prices.product .== p),
                           :producer_price])
price_gaps = DataFrame((product = p, from = a, to = b,
                        price_gap      = abs(pp(a, p) - pp(b, p)),
                        transport_cost = transport_costs[(p, a, b)])
                       for p in products for a in france for b in france
                       if a < b && !ismissing(pp(a, p)) && !ismissing(pp(b, p)))
price_gaps.ratio = price_gaps.price_gap ./ price_gaps.transport_cost
first(sort(price_gaps, :ratio, rev = true), 6)

# Ten of the sixty-one pairs have a gap above the freight, by a few units of
# price at most — 4 on a softwood log worth about 108, between `OF` and `SWF`.
#
# That is the price of `σ_france` being finite, and it is worth paying. With
# `sigma = Inf` both would vanish — no cross-hauling, no gap above freight —
# but where the French mills operate could be left undetermined, as
# [Perfect substitutes and constant returns](@ref) explains. The higher
# `σ_france`, the closer France comes to one market and the closer the model
# comes to that indeterminacy; 30 keeps it well away while leaving the four
# regions nearly interchangeable.

# ## Implementing scenarios
#
# A scenario is the same economy with something changed: a curve shifted, a
# stock smaller, a price outside the model higher. The economy is plain data,
# so a scenario is a copy of it with some fields replaced, solved again, and
# read against the reference run. One small helper does the copying:

reconfigure(d; kwargs...) =
    MarketData(; (f => get(kwargs, f, getfield(d, f)) for f in fieldnames(MarketData))...)
nothing # hide

# Three scenarios follow: a storm that floods the Landes with salvaged pine, a
# dearer gas that sends households to wood fuel, and a third of the French
# forests set aside.
#
# ### A storm in the Landes
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
# Here the salvage simply triples what the Landes offer at any price. The
# felled trees are also wood that can now be cut, so the storm adds them to
# the conifer stock of the Landes, as dead timber that can be salvaged — some
# 30 Mm³ of the 40 Klaus brought down. Everything else is left alone:

salvage = [s.region == "SWF" && s.product == "swr" ?
               SupplySpec(s; shift = 3.0) : s
           for s in supply]
nothing # hide

felled = merge(stock, Dict("SWF" => (conifers = stock["SWF"].conifers + 30.0,
                                     broadleaves = stock["SWF"].broadleaves)))
storm = solve_market(reconfigure(example_market, supply = salvage,
                                 limits = all_limits(felled)))
nothing # hide

# The question a spatial model exists to answer is how far the shock travels.
# Comparing the softwood price before and after, region by region:

price_of(r, p, result) = only(result.prices[(result.prices.region .== r) .&
                                            (result.prices.product .== p), :price])
nothing # hide

DataFrame((region = r,
           before = price_of(r, "swr", res),
           after  = price_of(r, "swr", storm),
           change_pct = 100 * (price_of(r, "swr", storm) / price_of(r, "swr", res) - 1))
          for r in regions)

# Two things are going on at once, and they are worth separating.
#
# The **gradient across the blocs** is the Armington structure at work. The
# salvage lands on France, passes with some loss into the rest of the EU, a
# close substitute, and barely reaches North America or Asia, which buy a
# different good. Had every product been homogeneous worldwide, a single world
# price would have moved everywhere at once; had every variety been equally
# distant, the shock would have stayed bottled up in `SWF`. Neither is what the
# forest sector looks like.
#
# The **gap inside France** is the evacuation limit. `SWF` loses seven times as
# much as its neighbours, which are near-perfect substitutes for it and should
# share in the glut. What holds the glut back is that the wood cannot get out:

storm.limits[storm.limits.limit .== "landes_evacuation", :]

# The ceiling that cost nothing in a normal year is now full, and its
# congestion rent is the price of a m³ of capacity to move wood out of the
# Landes. That rent is what drives a wedge between `SWF` and the rest of the
# French market, and it is why a salvage price collapse is a local event: not
# because the wood is worth less in principle, but because it cannot reach
# anyone who would pay more for it.

# ### When gas gets dearer
#
# The last of the exogenous things is the price of the fuels wood competes
# with. Doubling the gas price does not touch any curve in the model — gas is
# not in it — but it shifts the wood-fuel bundle's demand by
# `(60/30)^0.35 ≈ 1.27`, which [`exogenous_shift`](@ref) works out from the
# cross-price elasticity:

gas = solve_market(reconfigure(example_market, demand = build_demand(60.0)))

DataFrame((region = r,
           firewood = 100 * (quantity_of(gas.consumption, r, "smallwood") /
                             quantity_of(res.consumption, r, "smallwood") - 1),
           pellets  = 100 * (quantity_of(gas.consumption, r, "pellets") /
                             quantity_of(res.consumption, r, "pellets") - 1),
           roundwood_price = 100 * (price_of(r, "swr", gas) / price_of(r, "swr", res) - 1))
          for r in regions)

# Households burn much more wood, and the pull reaches back up the chain to the
# price of standing timber — which is how an energy policy becomes a forestry
# question.
#
# But the two columns do not move together, and the reason is the mills. The
# shifter raises demand for the *bundle*; which half of it grows depends on
# which half can grow. Pellets are made in plants whose capacity was fixed
# above:

gas.activity[gas.activity.process .== "pelletmill", :]

# Every one of them now runs full, and its capacity earns a rent. Pellets can
# grow only as far as the room the plants had spare, so most of the extra
# demand lands on the fuel that needs no plant: a log in a stove.
#
# Three of the features set up in this example meet in that one table — a
# bundle of substitutes in demand, a capacity bound on transformation, and a
# price from outside the model — and the result is not what any of them would
# have given alone.

# ### Setting forests aside
#
# The stands bind once there are fewer of them. Suppose a third of every
# French forest is set aside — protected, out of reach of the harvest — and
# nothing else changes:

protected = Dict(r => r in france ? (conifers = 2 / 3 * s.conifers, broadleaves = 2 / 3 * s.broadleaves) : s
                 for (r, s) in stock)
aside = solve_market(reconfigure(example_market, limits = all_limits(protected)))
nothing # hide

# Now the French regions hit their stands, and the limits that bind are those
# of all three products together: every m³ either kind of stand can give is
# cut.

aside.limits[aside.limits.rent .> 0, :]

# The rent of each is what one more m³ of standing timber would be worth to
# the economy: the stumpage value the scarcity creates, on top of what the
# supply curves ask. What it does to the harvest and to the price of logs:

DataFrame((region = r,
           harvest_before = sum(quantity_of(res.production, r, p) for p in primary),
           harvest_after  = sum(quantity_of(aside.production, r, p) for p in primary),
           log_price_change = 100 * (price_of(r, "swr", aside) / price_of(r, "swr", res) - 1))
          for r in regions)

# France cuts a sixth less wood, as much as its stands now allow, yet the
# price of its logs barely moves: French logs compete with those of the rest
# of Europe, a close substitute, which makes up most of the difference. The
# scarcity shows up instead in the rent — the gap between what a m³ fetches
# and what its owner would have sold it for — which is the value the protected
# stands take out of the harvest. A region whose wood had no close substitute
# would see it in its prices too.

# ## The limits of uncalibrated data
#
# Every curve in this model passes through a reference point taken from the
# data: each supply curve through a region's harvest at its roadside price,
# each demand curve through a region's consumption at the price given for the
# product. Yet the equilibrium does not land on those points. Here are the
# consumer prices it settles at, against the one price per product it was
# given:

DataFrame((product = p, given = reference[p][1],
           lowest  = minimum(price_of(r, p, res) for r in regions),
           highest = maximum(price_of(r, p, res) for r in regions))
          for p in sort(collect(keys(reference))))

# and the harvests, against those the supply curves were drawn through:

harvest_gap = DataFrame((product = s.product, region = s.region, given = s.q0,
                         model = quantity_of(res.production, s.region, s.product))
                        for s in supply)
harvest_gap.gap_pct = 100 .* (harvest_gap.model ./ harvest_gap.given .- 1)
first(sort(harvest_gap, :gap_pct, by = abs, rev = true), 6)

# The harvests land within 2% of the data: logs trade freely, and the supply
# curves hold them close to their reference points. The prices do not. They
# miss the one price given by up to a fifth, and by different amounts in
# different regions — plywood costs 306 to 339 where 285 was given, softwood
# lumber 209 to 245 against 205.
#
# The reason is that the numbers are not consistent with each other. One price
# per product says every region pays the same, while freight and the
# preference for local products say they cannot. The trade shares come from a
# rule, not from what was shipped, so the flows the model needs to feed each
# region are not those that made the observed harvests worth cutting. The
# mills' costs are guesses, so the margins they leave do not match the prices
# the mills paid and got. Each piece of data pulls the equilibrium its own way,
# and the equilibrium is wherever they balance — which is not the base year.
#
# Does it matter? A scenario is read as a difference from the model's own
# reference run, not from the data, and for much of it the answer is no:
# constant-elasticity curves give the same percentage response wherever they
# are evaluated. But a difference also depends on what the shock runs into —
# how close a limit is, how large each market share is, where the mills run —
# and those the uncalibrated reference run gets wrong in its own way. How much
# that changes the storm and the gas scenarios is measured on the third page.
#
# Getting the reference run onto the data is **calibration**: choosing the
# numbers nobody observes — trade shares, costs, the recipe of each mill, the
# reference point of each curve — so that the observed year is an equilibrium.
# The next page, [Calibrating with partial data](@ref "The forest sector: calibrating with partial data"),
# does it with the incomplete and inconsistent statistics a real application
# has to work with.
#
# The results of this page are kept, for that comparison:

results_dir = joinpath(pkgdir(Markets), "examples", "forest", "results")

function outcome(result, scenario)
    t = outerjoin(select(result.prices, :region, :product, :price, :producer_price),
                  rename(result.production, :quantity => :production),
                  select(result.consumption, :region, :product, :quantity => :consumption),
                  on = [:region, :product])
    return insertcols!(t, 1, :scenario => scenario)
end

function save_outcomes(file, runs)
    t = vcat((outcome(r, s) for (s, r) in runs)...)
    cell(x) = ismissing(x) ? "" : x isa Real ? string(round(x; sigdigits = 6)) : string(x)
    mkpath(results_dir)
    open(joinpath(results_dir, file), "w") do io
        println(io, join(names(t), ","))
        foreach(row -> println(io, join(cell.(collect(row)), ",")), eachrow(t))
    end
end

save_outcomes("uncalibrated.csv", ["reference" => res, "storm" => storm, "gas" => gas])

# ## Where to go from here
#
# * Change the elasticities in `σ_blocs` and `σ_europe` and watch the gradient
#   above widen or flatten, or those in `σ_lumber`, `σ_construction` and
#   `σ_fuel` and watch demand move between products instead of between
#   regions.
# * Raise or remove the pellet capacities and see how much of the energy shock
#   moves from firewood to pellets.
# * Put a tariff on something and compare the objective with and without the
#   revenue added back; the difference is the transfer the program counts as a
#   cost.
# * Put real elasticities in. [Buongiorno (2015)](@ref buongiorno2015) for the
#   demand side, [McDaniel & Balistreri (2003)](@ref mcdaniel2003) for the
#   Armington ones.
# * Add limits of your own with [`@limit`](@ref): a minimum share of
#   domestic logs in the French mills, a quota on exports, a cap on the
#   smallwood that may be burnt — any linear constraint on the quantities of
#   the model, with its rent reported in `Results.limits`.
# * Calibrate it: the [next page](@ref "The forest sector: calibrating with partial data")
#   replaces the gravity shares, the single prices and the guessed costs with
#   numbers worked out from data.
# * Give the EU market its own `destination` specification to represent a
#   non-tariff measure applying to one bloc only.
# * Make transport concave in distance, or compute distances between the
#   centres of gravity of the forest resource and of the population rather than
#   between region centroids.
# * Drive a multi-period run from outside: call `solve_market` once per period
#   and update each region's supply `shift` from the standing volume left after
#   the harvest, and each mill's capacity limit from what was built. The engine
#   stays static on purpose — the time loop belongs to whatever uses it.

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
