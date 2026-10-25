# # The forest sector: France in the world
#
# This example builds a complete spatial partial-equilibrium model of the
# forest-products sector with [Markets.jl](https://github.com/sylvaticus/Markets.jl)
# the way a real application would: starting from **observed base-year data**,
# adding the behavioural parameters one year of data cannot tell, deriving
# everything else by calibration, and then checking how close the model lands
# to the data it was built from before asking it anything new.
#
# !!! warning "For illustrative purposes only"
#     The data are plausible, but they are not statistics. They are rounded
#     the way published figures are, small flows are missing as they are from
#     trade statistics, and they do not balance exactly — like real data, and
#     unlike an equilibrium. The model therefore lands near them rather than
#     on them, and the page shows by how much.
#
# The structure of the model lives in this script and the data in the CSV
# files of `examples/forest/data`, as they should in any application: the same
# script runs on another base year by pointing it at other files.
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
# substitute for each other, the subject of
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

using Markets, DataFrames, DelimitedFiles

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

# They are connected by the mills set up further down:
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
# `chips` are a sawmill and plymill by-product *and* an input to three other
# mills, so the competition for residues is part of the equilibrium rather than
# an assumption about it.
#
# `smallwood` is **both a final product and an intermediate one**. It comes out
# of the forest, and from there it either goes into a household stove as
# firewood — final demand, like sawnwood or paper — or into a mill as a raw
# material. The model does not need to be told which: both uses face the same
# local price, and the split between them is whatever that price makes optimal.
# A product in this engine may be primary-supplied, manufactured, consumed and
# used as an input all at once; its balance simply sums the sources and the
# uses.

# ## Three kinds of numbers
#
# Everything the model needs falls into one of three groups, and keeping them
# apart is most of the discipline of building one.
#
# 1. **Observed**: what happened in the base year — what each region produced,
#    consumed and shipped, at what prices, with what plants. It comes from
#    statistics, and is read from files below.
# 2. **Chosen**: how the economy *responds*, which a single year cannot show —
#    elasticities of supply and demand, and how easily one product or one
#    origin replaces another. They come from the literature or from your own
#    estimation, and are written in this script, where they can be argued
#    about.
# 3. **Calibrated**: everything else — the reference points of the curves,
#    the shares of every CES, the costs of the mills. These are never guessed:
#    they are whatever makes the observed year an equilibrium of the chosen
#    behaviour.
#
# The model is then a machine for answering *what if*: it reproduces the base
# year, which you knew, and moves away from it consistently, which you did not.

# ## The base year, as observed
#
# The data are long-format CSV files, one observation per row, which is the
# shape statistics come in. `DelimitedFiles`, from the standard library, is
# enough to read them; in a real application you would more likely reach for
# CSV.jl.

data = pkgdir(Markets, "examples", "forest", "data")

function read_table(name)
    cells, header = readdlm(joinpath(data, name * ".csv"), ',', header = true)
    text(c) = c isa AbstractString ? String(c) : c
    return DataFrame([String(h) => identity.(text.(cells[:, j]))
                      for (j, h) in enumerate(vec(header))])
end

## a column of a table as a dictionary, keyed by the columns named
lookup(table, value, keys...) =
    Dict(Tuple(row[k] for k in keys) => row[value] for row in eachrow(table))
nothing # hide

# ### What was produced
#
# The harvest of each assortment and the output of each mill, by region. The
# French regions harvest some 19 Mm³ of softwood between them, the Landes
# (`SWF`) nearly half of it, while `GEF` and `OF` hold most of the hardwood.
# `SEF` and `GEF` make no pulp: there is no pulp mill there.

production = read_table("production")
unstack(production, :region, :product, :quantity)

# ### What producers were paid
#
# One price per product and region: roadside for logs, mill gate for the rest.
# This is what the region's **own variety** fetched, so there is no price
# where the region makes none.

producer_prices = read_table("producer_prices")
unstack(producer_prices, :region, :product, :price)

# ### What was shipped
#
# Bilateral flows between regions, as trade statistics report them — between
# blocs from customs data, between French regions from road-freight surveys,
# which show wood moving between them in both directions. Flows under a
# thousand m³ or tonnes are not reported. What a region keeps of its own
# production is not in the table either: no statistic reports it, and it
# follows from the other two, as we shall see.

trade = read_table("trade")
first(sort(trade, :quantity, rev = true), 8)

# ### What was consumed
#
# Final consumption of the finished products. Firewood (`smallwood` burnt as it
# comes) needs its own source — household energy surveys — because most of it
# never passes through a market that statistics see.

consumption = read_table("consumption")
unstack(consumption, :region, :product, :quantity)

# ### What it costs to carry
#
# Freight is charged per unit of quantity, not of value, because shipping a m³
# of logs costs much the same whether the logs are cheap or dear. Bulky,
# low-value products therefore carry a high transport cost relative to their
# price, which is what keeps roundwood trade regional and lets paper travel.
# Distances are relative, with the Europe–North America crossing as the unit.

freight   = read_table("freight")
distances = read_table("distance")

dist = lookup(distances, :distance, :from, :to)
distance(a, b) = get(dist, (a, b), get(dist, (b, a), 0.0))
rate = Dict(zip(freight.product, freight.cost))
transport_costs = Dict((p, a, b) => rate[p] * distance(a, b)
                       for p in products, a in regions, b in regions if a != b)
nothing # hide

# ### How the mills work
#
# The conversion factors are engineering, and the same everywhere: a m³ of
# sawlog yields half a m³ of lumber and a third of a m³ of chips, a tonne of
# paper takes 1.1 tonnes of pulp. They are written into the processes below.
#
# Mills that can choose between kinds of wood need one more thing: what they
# actually used. One region's mix is enough — the industry statistics of the
# rest of the EU here — because, as we shall see, the mix in any other region
# follows from its own prices.

input_mix = read_table("input_mix")

# And the plant that exists, where it limits anything: the pellet plants, which
# are recent, few, and ran with some room to spare.

capacities = read_table("capacity")

# ### Policies and logistics
#
# The rest of the world taxes its log exports, as several timber exporters do,
# to keep the wood at home for its own mills. The rate is **ad valorem**, a
# quarter of the price at the border, which is the one kind of wedge this model
# cannot write into its objective — the price is a dual of the problem. It is
# found by iteration instead, which [`solve_market`](@ref) does on its own; see
# [Why ad valorem rates on trade are special](@ref).

log_duty = 0.25
taxes = [Tax(name = "rw_log_duty", on = "trade", product = ["swr", "hwr"], from = "RW",
             rate = log_duty, relative = true)]
duty(p, from) = from == "RW" && p in ("swr", "hwr") ? log_duty : 0.0
nothing # hide

# Wood leaving the Landes has to be carried on the roads and through the ports
# that are there. A [`TradeLimit`](@ref) caps the **sum** of the flows it
# covers, so one entry limits all the roundwood `SWF` can send anywhere. At
# 18 Mm³ it is slack in the base year, and the storm scenario at the end is
# what makes it bite.

trade_limits = [TradeLimit(name = "landes_evacuation", products = ["swr", "smallwood"],
                           from = "SWF", capacity = 18.0)]
nothing # hide

# Finally the price of gas, which wood fuel competes with but this sector does
# not determine. Only its base-year level matters here: it is the reference the
# demand for wood fuel is calibrated at.

gas_reference = 30.0
nothing # hide

# ## What one year cannot tell
#
# A single observation of each market is a point. The curves through those
# points need slopes, and the slopes come from elsewhere.
#
# **Supply.** How much more wood the forest gives when it pays more. French
# forests are fragmented between many small private owners, which makes the
# harvest slow to respond to price, so their elasticities are lower than the
# blocs'.

supply_elasticity = Dict(
    "swr"       => Dict("SEF" => 0.40, "SWF" => 0.45, "GEF" => 0.40, "OF" => 0.40,
                        "EU"  => 0.60, "NA"  => 0.60, "AS"  => 0.55, "RW" => 0.55),
    "hwr"       => Dict("SEF" => 0.35, "SWF" => 0.35, "GEF" => 0.35, "OF" => 0.35,
                        "EU"  => 0.50, "NA"  => 0.50, "AS"  => 0.50, "RW" => 0.50),
    "smallwood" => Dict("SEF" => 0.50, "SWF" => 0.55, "GEF" => 0.50, "OF" => 0.50,
                        "EU"  => 0.60, "NA"  => 0.60, "AS"  => 0.55, "RW" => 0.55))
nothing # hide

# **Demand.** The own-price elasticities are **below 1**: demand for forest
# products is usually found to be price-inelastic, construction lumber and
# paper most of all, since their cost is a small part of the building or the
# printed product they end up in. An inelastic demand makes prices, rather
# than quantities, absorb a supply shock — which matters for the storm at the
# end.
#
# They are not given product by product, because the products are not separate
# markets. Softwood lumber and panels do the same job in a building and a
# builder swaps one for the other readily; hardwood lumber is a more distant
# substitute for either. Pellets and firewood are the same heat in different
# packaging. Paper competes with nothing else here. Each group gets one curve,
# with its own elasticity, and an elasticity of substitution `σ` saying how its
# members trade off against each other.

η_construction, η_fuel, η_paper = 0.45, 0.35, 0.30
σ_lumber, σ_construction, σ_fuel = 5.0, 1.6, 3.0
lumber = ["sawn_sw", "panel", "plywood"]
fuels  = ["smallwood", "pellets"]
nothing # hide

# Wood fuel also competes with gas, and the cross-price elasticity says how
# much: a dearer gas shifts the whole wood-fuel curve out.

gas_cross_elasticity = 0.35
nothing # hide

# **Mills.** How readily each mill switches between the kinds of wood it can
# use when one of them gets dearer.

σ_mills = Dict("plymill" => 1.8, "panelmill" => 2.5, "pulpmill" => 1.8, "pelletmill" => 3.0)
nothing # hide

# ### Who substitutes for whom
#
# Last, and most consequential: how the varieties of a product from different
# origins replace each other. If they were perfect substitutes everywhere,
# each market would simply buy from the cheapest source and a small cost
# advantage would take everything. If instead they were all imperfect
# substitutes with one elasticity, French oak would be no closer to French pine
# than to Asian hardwood, which is just as wrong.
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
# [`OriginNest`](@ref)s express: the French regions together at a very high
# elasticity, then France and the rest of the EU at a high one, then Europe as
# a whole against the other blocs at a lower one.

σ_blocs = Dict("swr" => 6.0, "hwr" => 6.0, "smallwood" => 7.0, "chips" => 7.0,
               "pulp" => 5.0, "sawn_sw" => 3.5, "sawn_hw" => 3.0, "panel" => 3.0,
               "plywood" => 2.5, "paper" => 4.0, "pellets" => 6.0)
σ_europe = Dict(p => 3σ for (p, σ) in σ_blocs)
σ_france = 30.0

europe(p) = OriginNest(sigma   = σ_europe[p],
                       origins = ["EU"],
                       nests   = [OriginNest(sigma = σ_france, origins = france)])
nothing # hide

# Commodities graded to a standard (chips, roundwood, pulp) substitute readily
# between blocs; manufactured products, where origin carries a reputation, do
# so less. Within Europe every product substitutes three times more easily than
# between blocs, and within France more easily still. Note the direction of the
# inequality the package enforces: a group must hold together *at least as
# tightly* as it holds to the outside, `30 ≥ 3σ ≥ σ`. The reverse would say
# that members of a group are worse substitutes for each other than for
# outsiders, which no cost function can represent.
#
# Why not make the French regions *perfect* substitutes, `sigma = Inf`, and
# France a single market up to transport? Because with mills that have
# constant returns, that would leave the model unable to say *where in France*
# the milling happens — see
# [Perfect substitutes and constant returns](@ref) for why, and for the
# remedies. A finite elasticity, however high, gives each region's product a
# demand of its own, which settles it.

# ## From data to model
#
# Now the third kind of number: every parameter that is neither observed nor
# chosen, worked out so that the base year is an equilibrium.

made     = lookup(production, :quantity, :region, :product)
πˢ       = lookup(producer_prices, :price, :region, :product)   # producer prices
consumed = lookup(consumption, :quantity, :region, :product)
nothing # hide

# ### Every flow, including the one nobody reports
#
# What a region keeps of its own variety is what it made less what it
# exported. With it the flows form a full origin × destination matrix per
# product, the diagonal included:

flows = lookup(trade, :quantity, :product, :from, :to)
exported(r, p) = sum(x for ((q, o, _), x) in flows if q == p && o == r; init = 0.0)
for ((r, p), q) in made
    flows[(p, r, r)] = q - exported(r, p)
end
minimum(values(flows))

# Not negative anywhere, which is the first check a base year has to pass: a
# region cannot export more than it makes.
#
# ### The price at the destination
#
# A variety costs its buyer what its producer was paid, plus the duty if there
# is one, plus the freight to get it there:

delivered(p, o, r) = o == r ? πˢ[(o, p)] :
                     πˢ[(o, p)] * (1 + duty(p, o)) + transport_costs[(p, o, r)]
nothing # hide

# ### The second price: what users pay
#
# Each region buys a mix of varieties, its own and imported ones, at different
# delivered prices, so **two prices** coexist in every region: what its
# producers get for their variety, read from the statistics above, and what its
# users pay for the mix. The second is not a separate observation. It is the
# average delivered price of what the region bought,
#
# ```math
# \pi^c_r = \frac{\sum_o X_{o,r}\, q_{o,r}}{\sum_o X_{o,r}},
# ```
#
# the value of its purchases divided by their volume — which is how a
# statistical office computes a consumer unit value, and also the price the
# calibrated model charges for the mix when it reproduces the base year.

πᶜ = Dict{Tuple{String,String},Float64}()
for r in regions, p in products
    bought = [(o, x) for ((q, o, d), x) in flows if q == p && d == r && x > 0]
    isempty(bought) && continue
    πᶜ[(r, p)] = sum(x * delivered(p, o, r) for (o, x) in bought) / sum(last, bought)
end

DataFrame((region = r, producer = πˢ[(r, "sawn_sw")], consumer = πᶜ[(r, "sawn_sw")])
          for r in regions)

# Softwood lumber, for instance. French mills get about 190 for it while
# French builders pay 210, because a good part of what they buy comes from
# the rest of the EU, dearer once delivered. In Asia-Pacific it is the other
# way round: imports are cheaper than the local product, and users pay less
# than local mills get.
#
# ### Trade shares
#
# A destination's composite is a CES over its origins, and its shares
# ``\delta`` are **not** the observed shares: they are the shares each origin
# would have *if all delivered prices were equal*, so they carry the home bias
# and the trade pattern that prices alone do not explain. Inverting the CES
# demand condition gives them from the flows and the delivered prices, level by
# level of the tree; [`armington_calibration`](@ref) does it for a whole
# product at once.
#
# It also sets one more number per destination, the `scale` of the composite.
# A CES of quantities equals their plain sum only where their prices are equal;
# elsewhere it falls short, by a few per cent for the French roundwood here.
# The scale restores the sum at the base year, so that the composite is
# counted in m³ and priced at ``\pi^c_r``, and the m³ the flows bring in are
# the m³ the mills and the households take out.

armington = [armington_calibration(
                 product = p, sigma = σ_blocs[p], nests = [europe(p)],
                 flows   = Dict((o, r) => x for ((q, o, r), x) in flows if q == p),
                 prices  = Dict((o, r) => delivered(p, o, r)
                                for o in regions, r in regions if haskey(πˢ, (o, p))))
             for p in products]
nothing # hide

# A route that carried nothing in the base year — or too little to be
# reported — gets no share, and stays closed whatever prices do: the
# well-known limit of Armington calibration.
#
# ### Supply
#
# Each forest supply curve runs through the observed harvest at the observed
# producer price, with the chosen slope:

supply = [SupplySpec(product = p, region = r, p0 = πˢ[(r, p)], q0 = made[(r, p)],
                     elasticity = supply_elasticity[p][r])
          for p in keys(supply_elasticity) for r in regions]
nothing # hide

# ### Final demand
#
# A curve over a bundle of products is a [`DemandNest`](@ref): a tree of
# products with one demand curve on the bundle. It has to be a tree rather
# than a matrix of cross-price elasticities, for the reason set out in
# [Substitution in final demand](@ref) — a free matrix is not the gradient of
# any benefit function, and the model would stop being the market it claims to
# represent.
#
# #### Calibrating a bundle
#
# A demand curve over a bundle needs three numbers that no statistic reports:
# the **shares** of its members, and the **price** and **quantity** of the
# bundle itself, which are the `p0` and `q0` of its curve. What the data give
# instead is a consumer price and a quantity for each product.
# [`ces_calibration`](@ref) turns the second into the first, choosing the three
# so that at the base-year prices the model buys exactly the base-year
# quantities. For the wood fuel of `OF`:

ces_calibration(quantities = [consumed[("OF", p)] for p in fuels],
                prices     = [πᶜ[("OF", p)] for p in fuels], sigma = σ_fuel)

# * The **shares** are what each product's share would be if both cost the
#   same. Pellets get 96% of the bundle although households in `OF` buy three
#   times more firewood: pellets are bought little *because* they are dearer,
#   and would be bought far more at an equal price.
# * The **price** is the CES price index of the bundle — the cost of one unit of
#   it at the base-year prices.
# * The **quantity** is the base-year spending divided by that price. It is an
#   index number, not a volume: m³ of firewood and tonnes of pellets would mean
#   nothing added up, while their values at a common price do. That is what a
#   CES composite is.
#
# The prices are **consumer** prices, ``\pi^c``: the household in `OF` buys the
# mix of French and imported pellets, not the pellets of `OF`.
#
# #### Wood in buildings
#
# Construction is a nest inside a nest — softwood lumber, particleboard and
# plywood together at `σ_lumber`, hardwood lumber against that group at
# `σ_construction`. The three in the inner group are what a builder actually
# chooses between for a floor or a wall; appearance-grade hardwood is a
# different kind of purchase. A nested bundle is calibrated from the inside
# out: the inner group on its three products; its price and quantity then make
# it **one member** of the outer bundle, next to hardwood lumber; and the outer
# price and quantity are the reference point of the curve.

function construction_demand(r)
    inner = ces_calibration(quantities = [consumed[(r, p)] for p in lumber],
                            prices     = [πᶜ[(r, p)] for p in lumber], sigma = σ_lumber)
    outer = ces_calibration(quantities = [consumed[(r, "sawn_hw")], inner.quantity],
                            prices     = [πᶜ[(r, "sawn_hw")], inner.price],
                            sigma = σ_construction)
    return DemandSpec(region = r, p0 = outer.price, q0 = outer.quantity,
        elasticity = η_construction,
        product = DemandNest(sigma = σ_construction,
            products = ["sawn_hw"], shares = [outer.shares[1]],
            nests = [DemandNest(sigma = σ_lumber, share = outer.shares[2],
                                products = lumber, shares = inner.shares)]))
end
nothing # hide

# `η_construction` is the own-price elasticity of the whole bundle: how much
# less wood goes into buildings when all of it gets dearer. How the bundle
# splits when only one product gets dearer is the job of the two `σ`.
#
# #### Wood for heat
#
# Wood fuel is a single nest, and the place where the price of gas comes in.
# An exogenous price needs no nest at all — being a constant, it shifts the
# whole curve — so the cross-price elasticity goes into
# [`exogenous_shift`](@ref). At the reference level the shifter is exactly 1,
# and the curve is the calibrated one; the scenario at the end moves it.

function fuel_demand(r, gas = gas_reference)
    fuel = ces_calibration(quantities = [consumed[(r, p)] for p in fuels],
                           prices     = [πᶜ[(r, p)] for p in fuels], sigma = σ_fuel)
    return DemandSpec(region = r, p0 = fuel.price, q0 = fuel.quantity, elasticity = η_fuel,
        shift = exogenous_shift(gas = (level = gas, reference = gas_reference,
                                       elasticity = gas_cross_elasticity)),
        product = DemandNest(sigma = σ_fuel, products = fuels, shares = fuel.shares))
end
nothing # hide

# #### All of final demand
#
# Paper competes with nothing here, so it keeps a plain curve through its own
# consumer price and consumption. The whole demand side is then one curve per
# bundle and region, built by a function of the gas price so that the gas
# scenario at the end can rebuild it at another level:

build_demand(gas = gas_reference) =
    [[construction_demand(r) for r in regions];
     [fuel_demand(r, gas) for r in regions];
     [DemandSpec(product = "paper", region = r, p0 = πᶜ[(r, "paper")],
                 q0 = consumed[(r, "paper")], elasticity = η_paper) for r in regions]]

demand = build_demand()
nothing # hide

# ### The mills
#
# Seven processes, of three kinds.
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
# **Fixed proportions, single input.** The paper mill takes pulp and nothing
# else.
#
# A bundle is calibrated like a demand nest, from the observed mix at the
# consumer prices of the region it was observed in. Its quantities are all m³
# of wood, so here they *do* add up, and the `scale` of the calibration lets
# the bundle keep counting m³ — the `composite` is then the m³ of wood a unit
# of output takes, as the conversion factor says it:

function bundle(process)
    mix = input_mix[input_mix.process .== process, :]
    c = ces_calibration(quantities = mix.quantity,
                        prices = [πᶜ[(r, p)] for (r, p) in zip(mix.region, mix.product)],
                        sigma = σ_mills[process])
    return ces(composite = sum(mix.quantity), products = mix.product, shares = c.shares,
               sigma = σ_mills[process], phi = c.scale)
end

bundle("pelletmill")

# The shares are the same wherever the process runs, so a pellet mill in
# `SWF`, facing other prices, uses another mix — the one those prices make
# cheapest. That is why one region's mix was enough. A region whose observed
# mix departed from what its prices imply would need a process of its own.
#
# The technology, then, is the conversion factors and those bundles:

technology = [
    ("sawmill_sw", [leontief(product = "swr", coeff = 1.0)], ["sawn_sw" => 0.50, "chips" => 0.35]),
    ("sawmill_hw", [leontief(product = "hwr", coeff = 1.0)], ["sawn_hw" => 0.45, "chips" => 0.35]),
    ("plymill",    [bundle("plymill")],                      ["plywood" => 1.0,  "chips" => 0.8]),
    ("panelmill",  [bundle("panelmill")],                    ["panel"   => 1.0]),
    ("pulpmill",   [bundle("pulpmill")],                     ["pulp"    => 1.0]),
    ("papermill",  [leontief(product = "pulp", coeff = 1.10)], ["paper" => 1.0]),
    ("pelletmill", [bundle("pelletmill")],                   ["pellets" => 1.0]),
]
nothing # hide

# #### What the mills cost
#
# The `vacost` of a process is the cost of everything not modelled as an
# input: labour, energy, capital. No statistic gives it in the units the model
# wants, but the base year implies it. A mill that ran and did not expand
# earned no more than its costs, so its cost is its **margin**: what its
# outputs fetched at the producer prices of its region, less what its inputs
# cost at the consumer prices there.
#
# For a CES bundle the cost of the inputs is the price index of the bundle at
# the region's prices — which is what the cheapest mix costs there — per m³ of
# the bundle, hence the division by its scale:

price_index(δ, p, σ) = sum(δ .* p .^ (1 - σ))^(1 / (1 - σ))
input_cost(n::Nest, r) = length(n.products) == 1 ?
    n.composite * πᶜ[(r, only(n.products))] :
    n.composite * price_index(n.shares, [πᶜ[(r, p)] for p in n.products], n.sigma) / n.phi
margin(inputs, outputs, r) = sum(y * πˢ[(r, k)] for (k, y) in outputs) -
                             sum(input_cost(n, r) for n in inputs)

operating(outputs) = [r for r in regions if get(made, (r, first(first(outputs))), 0.0) > 0]
DataFrame((; region = r,
             (Symbol(name) => r in operating(outputs) ? margin(inputs, outputs, r) : missing
              for (name, inputs, outputs) in technology)...)
          for r in regions)

# The margins differ from region to region, as wages and energy do, which is
# why a process takes a `region => cost` dictionary. They are positive
# everywhere, the second check a base year has to pass: a negative one would
# say a mill ran at a loss, and the model would close it.
#
# The margin is taken as all cost. A mill that ran flat out could have earned
# a rent on top, and one year cannot tell the two apart; that is why the
# capacities in the data are those of plants with room to spare.
#
# A process operates where its main output was made — which is what keeps
# pulp out of `SEF` and `GEF` — with the capacities read above:

capacity = Dict(name => Dict(row.region => row.capacity
                             for row in eachrow(capacities) if row.process == name)
                for name in unique(capacities.process))

processes = [Process(name = name, inputs = inputs, outputs = outputs,
                     regions  = operating(outputs),
                     vacost   = Dict(r => margin(inputs, outputs, r) for r in operating(outputs)),
                     capacity = get(capacity, name, Inf))
             for (name, inputs, outputs) in technology]
nothing # hide

# ## Solving

example_market = MarketData(; regions, products, demand, supply, processes,
                              transport_costs, armington, trade_limits, taxes)

res = solve_market(example_market)
nothing # hide

# ## How close does the model land to the data?
#
# Every curve, share and cost above was fitted to one observation, and each
# passes through it exactly: the supply curve of `SWF` through its harvest at
# its producer price, the demand for paper in `OF` through its consumption at
# its consumer price, the trade shares of each destination through its flows.
# But the observations come from different sources, rounded, with small flows
# missing, and they are not quite consistent with each other — the wood a
# region took in does not add up to exactly what its households and its mills
# used. An equilibrium has to satisfy every balance at once, so it settles
# near every observation and on none.
#
# How near is the first thing to look at, every time — *benchmark replication*:

modelled(table, value, keys...) = lookup(table, value, keys...)
relative_gaps(model, observed) = sort([abs(get(model, k, 0.0) / v - 1) for (k, v) in observed])
shipped = Dict(k => x for (k, x) in flows if k[2] != k[3])

comparison = [
    "producer prices" => (modelled(res.prices, :producer_price, :region, :product), πˢ),
    "consumer prices" => (modelled(res.prices, :price, :region, :product), πᶜ),
    "production"      => (modelled(res.production, :quantity, :region, :product), made),
    "consumption"     => (modelled(res.consumption, :quantity, :region, :product), consumed),
    "trade"           => (modelled(res.trade, :quantity, :product, :from, :to), shipped)]

DataFrame((observed = name, entries = length(observed),
           median_gap = gaps[cld(end, 2)], largest_gap = last(gaps))
          for (name, (model, observed)) in comparison
          for gaps in (relative_gaps(model, observed),))

# A tenth of a per cent in the middle, under one per cent at worst for prices
# and quantities, a few per cent for some trade flows: these are gaps of the
# size of the rounding, and the largest sit on the smallest flows, where the
# last digit is a large share. That says the data are consistent up to their
# own precision, and the calibration did its job.
#
# Gaps of ten or twenty per cent would say something else: that the data
# contradict each other — a balance that does not close, a price from another
# year, a consumption survey at odds with the trade figures — or that a
# calibration step went wrong, a unit, a sign, a forgotten duty. Either is
# worth chasing before asking the model anything, because every scenario
# starts from where the model landed, not from the data. See
# [Calibration](@ref) for more on what these gaps mean.

# ## Prices
#
# The two prices of each region, as the model has them.

unstack(res.prices, :region, :product, :price)

#-

unstack(res.prices, :region, :product, :producer_price)

# The `missing` entries are not a failure: `SEF` and `GEF` have no pulp mill,
# so they have no pulp of their own to price. They buy it, and the first table
# gives what they pay.
#
# Notice how close the four French rows are to each other in both tables, and
# how far they are from `NA` or `AS`. That is the nest at work.

# ## What binds, and what it is worth
#
# The pellet plants run below capacity everywhere, so in the base year another
# unit of plant would earn nothing:

res.activity[res.activity.process .== "pelletmill", :]

# The Landes evacuation limit is slack too.

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
# price:

quantity_of(table, r, p) = (rows = table[(table.region .== r) .& (table.product .== p),
                                            :quantity];
                            isempty(rows) ? 0.0 : only(rows))
crossing(direction, r) = sum(res.trade[(res.trade.product .== "smallwood") .&
                                       (res.trade[!, direction] .== r), :quantity]; init = 0.0)

smallwood_use = DataFrame(
    (region      = r,
     harvest     = quantity_of(res.production, r, "smallwood"),
     firewood    = quantity_of(res.consumption, r, "smallwood"),
     net_exports = crossing(:from, r) - crossing(:to, r))
    for r in regions)
smallwood_use.to_the_mills =
    smallwood_use.harvest .- smallwood_use.firewood .- smallwood_use.net_exports
smallwood_use

# The shares are worth reading off. North America burns an eighth of what it
# cuts, and imports wood for its mills on top; Asia-Pacific and the rest of the
# world burn a third or more. The French regions burn between a fifth (`SWF`,
# `GEF`) and two fifths (`OF`, where most of the population is), and `SEF` and
# `GEF` **export most of what is left rather than milling it**, because French
# industry is small next to its neighbours' and the wood is worth more across
# the border.
#
# That is the base year, so the model returns it because it was calibrated on
# it. What the model adds is the mechanism: the same wood is offered to the
# stove and to the mill at one price, and when anything moves, the split moves
# with whichever use outbids the other — as the gas scenario below shows.

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

# The widest gaps are on manufactured products, a few per cent of their price:
# paper made in `OF` and in `SEF` is, to this model, not quite the same paper.
# Roundwood and smallwood stay within a whisker of the freight.
#
# That is the price of `σ_france` being finite, and it is worth paying. With
# `sigma = Inf` both would vanish — no cross-hauling, no gap above freight —
# but where the French mills operate would be left undetermined, as
# [Perfect substitutes and constant returns](@ref) explains. The higher
# `σ_france`, the closer France comes to one market and the closer the model
# comes to that indeterminacy; 30 keeps it well away while leaving the four
# regions nearly interchangeable.

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
# The **gap inside France** is the Landes' own. `SWF` loses nearly three times
# as much as its neighbours, which are near-perfect substitutes for it and
# should share in the glut. What holds the glut back is that the wood cannot
# get out:

storm.limits

# The ceiling that cost nothing in the base year is now full, and its
# congestion rent is the price of a m³ of capacity to move wood out of the
# Landes. What cannot leave has to be used where it lies, and the Landes
# mills take what they can:

storm.activity[storm.activity.region .== "SWF", :]

# The softwood sawmills nearly double their intake, and the panel and pulp
# mills grow too; only the pellet plant, already at its capacity, cannot. It
# is not enough to clear the salvage at the old price, which is why a salvage
# price collapse is a local event: not because the wood is worth less in
# principle, but because it cannot reach anyone who would pay more for it.

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

# Households burn a third more wood, and the pull reaches back up the chain
# to the price of standing timber, 3% dearer everywhere — which is how an
# energy policy becomes a forestry question.
#
# But the two columns do not move together, and the reason is the mills. The
# shifter raises demand for the *bundle*; which half of it grows depends on
# which half can grow. Pellets are made in plants:

gas.activity[gas.activity.process .== "pelletmill", :]

# Every one of them now runs full and earns a rent of some 30 a tonne. They
# could grow only into the room they had spare, so pellet consumption grows
# least in the rest of the world, whose plants were the fullest, and most in
# North America, whose plants had the most room. The rest of the extra demand
# lands on the fuel that needs no plant: a log in a stove.
#
# Three of the features set up in this example meet in that one table — a
# bundle of substitutes in demand, a capacity bound on transformation, and a
# price from outside the model — and the result is not what any of them would
# have given alone.

# ## Where to go from here
#
# * Point `data` at your own base year. The script does not care where the
#   numbers came from, only that they balance; the replication table says
#   whether they do.
# * Change the elasticities in `σ_blocs` and `σ_europe` and watch the storm's
#   gradient widen or flatten, or those in `σ_lumber`, `σ_construction` and
#   `σ_fuel` and watch demand move between products instead of between
#   regions. None of them changes the base year — that is what calibration
#   guarantees — only how the economy moves away from it.
# * Raise `σ_france` and watch France draw together into one market, or lower
#   it and watch the regions drift apart; remove the pellet capacities and see
#   how much of the energy shock moves from firewood to pellets.
# * Put real elasticities in. [Buongiorno (2015)](@ref buongiorno2015) for the
#   demand side, [McDaniel & Balistreri (2003)](@ref mcdaniel2003) for the
#   Armington ones.
# * Give the EU market its own `destination` specification to represent a
#   non-tariff measure applying to one bloc only.
# * Drive a multi-period run from outside: call `solve_market` once per period
#   and update each region's supply `shift` from the standing volume left after
#   the harvest, and each mill's `capacity` from what was built. The engine
#   stays static on purpose — the time loop belongs to whatever uses it, and
#   the calibration above stays where it is for the whole run.

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
