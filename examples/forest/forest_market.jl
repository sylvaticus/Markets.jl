# # The forest sector: France in the world
#
# This example builds a complete spatial partial-equilibrium model of the
# forest-products sector with [Markets.jl](https://github.com/sylvaticus/Markets.jl)
# and solves it.
#
# !!! warning "For illustrative purposes only"
#     The data in this example are plausible, but they have not been
#     calibrated to any real-world case.
#
#     Also, for simplicity, the model structure and the data are defined
#     together here. In a real application, you will want to keep the
#     _structure_ of your model separate from the _data_, and read the data
#     from a separate file.
#
# The economy has **eight regions**. France is split into four regions, within
# which products are homogeneous, and the rest of the world into four blocs:
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
# The sector modelled here is the one the Global Forest Products Model
# established as the reference, [Buongiorno et al. (2003)](@ref gfpm2003); for
# the same economics in working open-source code, with different choices about
# dynamics and data, see TiMBA, [[C3]](@ref C3).
#
# All the figures are illustrative — plausible orders of magnitude, not
# calibrated statistics. Quantities are in million m³ for wood products and
# million tonnes for pulp and paper; prices are in USD per m³ or per tonne.

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
# reference point. The French regions harvest some 27 Mm³ of industrial
# roundwood between them: the Landes (`SWF`) dominate softwood, while `GEF`
# and `OF` hold most of the hardwood. Their supply elasticities are lower than
# the blocs' — French forests are fragmented between many small private owners,
# which makes the harvest slow to respond to price.

supply = SupplySpec[
    ## softwood roundwood
    SupplySpec(product = "swr", region = "SEF", p0 = 72, q0 =   3.2, elasticity = 0.40),
    SupplySpec(product = "swr", region = "SWF", p0 = 65, q0 =   7.5, elasticity = 0.45),
    SupplySpec(product = "swr", region = "GEF", p0 = 68, q0 =   2.6, elasticity = 0.40),
    SupplySpec(product = "swr", region = "OF",  p0 = 70, q0 =   3.4, elasticity = 0.40),
    SupplySpec(product = "swr", region = "EU",  p0 = 68, q0 = 270.0, elasticity = 0.60),
    SupplySpec(product = "swr", region = "NA",  p0 = 62, q0 = 380.0, elasticity = 0.60),
    SupplySpec(product = "swr", region = "AS",  p0 = 80, q0 = 150.0, elasticity = 0.55),
    SupplySpec(product = "swr", region = "RW",  p0 = 58, q0 = 190.0, elasticity = 0.55),
    ## hardwood roundwood
    SupplySpec(product = "hwr", region = "SEF", p0 = 78, q0 =   1.4, elasticity = 0.35),
    SupplySpec(product = "hwr", region = "SWF", p0 = 75, q0 =   1.2, elasticity = 0.35),
    SupplySpec(product = "hwr", region = "GEF", p0 = 80, q0 =   3.4, elasticity = 0.35),
    SupplySpec(product = "hwr", region = "OF",  p0 = 78, q0 =   4.6, elasticity = 0.35),
    SupplySpec(product = "hwr", region = "EU",  p0 = 76, q0 =  85.0, elasticity = 0.50),
    SupplySpec(product = "hwr", region = "NA",  p0 = 72, q0 =  70.0, elasticity = 0.50),
    SupplySpec(product = "hwr", region = "AS",  p0 = 74, q0 = 180.0, elasticity = 0.50),
    SupplySpec(product = "hwr", region = "RW",  p0 = 66, q0 = 130.0, elasticity = 0.50),
    ## smallwood: thinnings and tops, worth a fraction of a sawlog
    SupplySpec(product = "smallwood", region = "SEF", p0 = 42, q0 =   6.0, elasticity = 0.50),
    SupplySpec(product = "smallwood", region = "SWF", p0 = 38, q0 =   3.5, elasticity = 0.55),
    SupplySpec(product = "smallwood", region = "GEF", p0 = 40, q0 =   6.5, elasticity = 0.50),
    SupplySpec(product = "smallwood", region = "OF",  p0 = 41, q0 =   9.0, elasticity = 0.50),
    SupplySpec(product = "smallwood", region = "EU",  p0 = 39, q0 = 150.0, elasticity = 0.60),
    SupplySpec(product = "smallwood", region = "NA",  p0 = 34, q0 =  70.0, elasticity = 0.60),
    SupplySpec(product = "smallwood", region = "AS",  p0 = 45, q0 = 200.0, elasticity = 0.55),
    SupplySpec(product = "smallwood", region = "RW",  p0 = 30, q0 = 250.0, elasticity = 0.55),
]
nothing # hide

# ## Final demand
#
# Only finished products are consumed. The reference quantities are roughly
# proportional to population and income, so `OF` — which holds two thirds of
# the French population — is much the largest French market; the reference
# prices are those of the base year.

reference = Dict(
    ## product    => (price, quantities by region)
    "sawn_sw"   => (250.0,  Dict("SEF" => 1.90, "SWF" => 0.90, "GEF" => 0.80, "OF" => 5.60,
                                "EU" => 75.0, "NA" => 120.0, "AS" => 80.0,  "RW" => 45.0)),
    "sawn_hw"   => (300.0,  Dict("SEF" => 0.25, "SWF" => 0.12, "GEF" => 0.11, "OF" => 0.72,
                                "EU" => 12.0, "NA" =>  20.0, "AS" => 60.0,  "RW" => 15.0)),
    "panel"     => (350.0,  Dict("SEF" => 0.90, "SWF" => 0.40, "GEF" => 0.40, "OF" => 2.80,
                                "EU" => 45.0, "NA" =>  40.0, "AS" => 90.0,  "RW" => 25.0)),
    ## plywood: veneers glued crosswise, dearer than particleboard and used
    ## where a panel has to carry a load
    "plywood"   => (450.0,  Dict("SEF" => 0.12, "SWF" => 0.06, "GEF" => 0.05, "OF" => 0.35,
                                "EU" =>  6.0, "NA" =>  14.0, "AS" => 60.0,  "RW" =>  9.0)),
    "paper"     => (1000.0, Dict("SEF" => 1.80, "SWF" => 0.80, "GEF" => 0.75, "OF" => 5.60,
                                "EU" => 70.0, "NA" =>  80.0, "AS" => 130.0, "RW" => 40.0)),
    ## firewood: smallwood bought and burnt as it comes, the oldest use of all
    "smallwood" => (55.0,   Dict("SEF" => 3.20, "SWF" => 1.50, "GEF" => 2.60, "OF" => 7.00,
                                "EU" => 90.0, "NA" =>  35.0, "AS" => 150.0, "RW" => 190.0)),
    ## pellets: the same energy, dried, ground and compressed
    "pellets"   => (250.0,  Dict("SEF" => 0.50, "SWF" => 0.25, "GEF" => 0.30, "OF" => 1.40,
                                "EU" => 22.0, "NA" =>   9.0, "AS" =>  12.0, "RW" =>   4.0)),
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
# quantities. For the wood fuel of `OF`, 7 Mm³ of firewood at 55 and 1.4 of
# pellets at 250:

ces_calibration(quantities = [7.0, 1.4], prices = [55.0, 250.0], sigma = 3)

# * The **shares** are what each product's share would be if both cost the
#   same. Pellets get 95% of the bundle although they are a sixth of the
#   quantity bought: they are bought little *because* they are four and a half
#   times dearer, and would be bought far more at an equal price.
# * The **price** is the CES price index of the bundle — the cost of one unit of
#   it at the base-year prices.
# * The **quantity** is the base-year spending, 735, divided by that price. It
#   is an index number, not a volume: m³ of firewood and tonnes of pellets
#   would mean nothing added up, while their values at a common price do. That
#   is what a CES composite is.
#
# Feed those back into the CES demand condition at the base-year prices and
# out come exactly 7 and 1.4.
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
# something to guess at — doubling plywood demand in this example leaves the
# chip price almost exactly where it was, the extra residues offsetting the
# extra competition for them, while hardwood logs rise about 10% and pulp
# about 2%. The log market wins. That is the kind of link a sector model
# exists to settle.
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
# [`Process`](@ref) with constant returns and no `capacity` would pile into
# whichever region makes pellets most cheaply; the installed plant is what
# stops it, region by region:

pellet_capacity = Dict("SEF" => 0.35, "SWF" => 0.50, "GEF" => 0.35, "OF" => 1.10,
                       "EU"  => 20.0, "NA"  => 12.0, "AS"  => 15.0, "RW" => 10.0)

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
            capacity = pellet_capacity,
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

σ_blocs = Dict("swr" => 6.0, "hwr" => 6.0, "smallwood" => 7.0, "chips" => 7.0,
               "pulp" => 5.0, "sawn_sw" => 3.5, "sawn_hw" => 3.0, "panel" => 3.0,
               "plywood" => 2.5, "paper" => 4.0, "pellets" => 6.0)

# Commodities graded to a standard (chips, roundwood, pulp) substitute readily
# between blocs; manufactured products, where origin carries a reputation, do
# so less. Within Europe every product substitutes three times more easily than
# between blocs:

σ_europe = Dict(p => 3σ for (p, σ) in σ_blocs)
nothing # hide

europe(p) = OriginNest(sigma   = σ_europe[p],
                       origins = ["EU"],
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

# ## What cannot move, and what is taxed
#
# Two final things stand between the economy as described and the equilibrium.
#
# Wood leaving the Landes has to be carried on the roads and through the ports
# that are there. A [`TradeLimit`](@ref) caps the **sum** of the flows it
# covers, so one entry limits all the roundwood `SWF` can send anywhere. At
# 18 Mm³ it is slack in a normal year — the region ships about 12 — and the
# storm scenario at the end is what makes it bite.

trade_limits = [TradeLimit(name = "landes_evacuation", products = ["swr", "smallwood"],
                           from = "SWF", capacity = 18.0)]

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
                              transport_costs, armington, trade_limits, taxes)

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
# Three of the things set up above are limits rather than curves, and the
# solution says what each is worth at the margin.
#
# The pellet mills first. Where the installed plant is enough, the capacity
# shows no rent; where it is not, the bound binds and the rent is what another
# tonne of annual capacity would earn its owner — the signal that would, in a
# model with investment, build the next mill:

res.activity[res.activity.process .== "pelletmill", :]

# The Landes evacuation limit is slack in a normal year: the region ships about
# 12 Mm³ of roundwood against a ceiling of 18, so the ceiling costs nothing and
# earns nothing.

res.limits

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

# ## Is France really one market?
#
# The French regions were declared perfect substitutes for each other, so they
# should behave between themselves exactly as the homogeneous model of
# Samuelson does. That means two things, and both are worth checking rather
# than assuming.
#
# First, **no cross-hauling**: two regions never ship the same product to each
# other at once, because for an identical good that would be a pure waste of
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
    MarketData(; (f => get(kwargs, f, getfield(d, f)) for f in fieldnames(MarketData))...)
nothing # hide

salvage = [s.region == "SWF" && s.product == "swr" ?
               SupplySpec(product = "swr", region = "SWF", p0 = s.p0, q0 = s.q0,
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
# The **gap inside France** is the evacuation limit. `SWF` loses far more than
# its neighbours, which were declared perfect substitutes for it — a pool that
# should share one price up to transport. What breaks the pool is that the wood
# cannot get out:

storm.limits

# The ceiling that cost nothing in a normal year is now full, and its
# congestion rent is the price of a m³ of capacity to move wood out of the
# Landes. That rent is exactly what drives a wedge between `SWF` and the rest
# of the French market, and it is why a salvage price collapse is a local
# event: not because the wood is worth less in principle, but because it cannot
# reach anyone who would pay more for it.

# ## When gas gets dearer
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
# above, and in the regions where that capacity already binds they can barely
# expand at all, so most of the extra demand lands on the fuel that needs no
# plant: a log in a stove. In the rest of the world, where the mills are
# capped and the competition for wood is fiercest, pellet consumption actually
# *falls*.
#
# Three of the features set up in this example meet in that one table — a
# bundle of substitutes in demand, a capacity bound on transformation, and a
# price from outside the model — and the result is not what any of them would
# have given alone.

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
# * Replace the gravity shares with a base-year trade matrix. The `shares` are
#   the shares *at equal delivered prices*, not the observed ones, so pass the
#   matrix through [`armington_shares`](@ref) rather than using it directly.
# * Give the EU market its own `destination` specification to represent a
#   non-tariff measure applying to one bloc only.
# * Make transport concave in distance, or compute distances between the
#   centres of gravity of the forest resource and of the population rather than
#   between region centroids.
# * Drive a multi-period run from outside: call `solve_market` once per period
#   and update each region's supply `shift` from the standing volume left after
#   the harvest, and each mill's `capacity` from what was built. The engine
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
