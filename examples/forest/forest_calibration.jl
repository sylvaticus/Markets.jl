# # The forest sector: calibrating with partial data
#
# The [first page](@ref "The forest sector: France in the world") wrote the
# model with the numbers at hand before any calibration effort, and its
# reference run landed some way from them. This page builds **the same model**
# — the same regions, products, mills and elasticities — from the statistics a
# real application would actually have: incomplete, rounded, and not quite
# consistent with each other. The aim is a reference run that sits on the
# data, as near as the data allow.
#
# The general ideas are in [Calibration](@ref) and
# [Calibrating with partial data](@ref); this page is the worked case.
#
# !!! warning "For illustrative purposes only"
#     The data are plausible, but they are not statistics. They are rounded
#     the way published figures are, small flows are missing as they are from
#     trade statistics, and they do not balance exactly — like real data.

using Markets, DataFrames, DelimitedFiles

france  = ["SEF", "SWF", "GEF", "OF"]
blocs   = ["EU", "NA", "AS", "RW"]
regions = [france; blocs]
products = ["swr", "hwr", "smallwood", "chips", "pulp",
            "sawn_sw", "sawn_hw", "panel", "plywood", "paper", "pellets"]
nothing # hide

# ## Three kinds of numbers
#
# Everything the model needs falls into one of three groups, and keeping them
# apart is most of the discipline of calibration.
#
# 1. **Observed**: what happened in the base year — what each region produced,
#    used and consumed, at what prices, what crossed borders, what plant there
#    was. It comes from statistics, and is read from files below.
# 2. **Chosen**: how the economy *responds*, which a single year cannot show —
#    elasticities of supply and demand, and how easily one product or one
#    origin replaces another. They come from the literature or from your own
#    estimation, and are written in this script, where they can be argued
#    about.
# 3. **Calibrated**: everything else — the reference points of the curves,
#    the shares of every CES, the costs of the mills, and here also the flows
#    the statistics do not report. These are never guessed: they are whatever
#    makes the observed year an equilibrium of the chosen behaviour.

# ## The statistics at hand
#
# The data are long-format CSV files, one observation per row, the shape
# statistics come in. `DelimitedFiles`, from the standard library, is enough to
# read them; in a real application you would more likely reach for CSV.jl.

data = pkgdir(Markets, "examples", "forest", "data")

function read_table(file)
    cells, header = readdlm(joinpath(data, file), ',', header = true)
    text(c) = c isa AbstractString ? String(c) : c
    return DataFrame([String(h) => identity.(text.(cells[:, j]))
                      for (j, h) in enumerate(vec(header))])
end

## a column of a table as a dictionary, keyed by the columns named
lookup(table, value, keys...) =
    Dict(Tuple(row[k] for k in keys) => row[value] for row in eachrow(table))
nothing # hide

# **Production**, from forest and industry statistics: the harvest of each
# assortment and the output of each mill, by region.

production = read_table("partial/production.csv")
unstack(production, :region, :product, :quantity)

# **Producer prices**: roadside for logs, mill gate for the rest, by region.
# This is what each region's own product fetched, so there is none where a
# region makes none.

producer_prices = read_table("partial/producer_prices.csv")
unstack(producer_prices, :region, :product, :price)

# **Use**: what each region took in of each product, for its households and
# its mills together, from consumption and industry surveys.

use = read_table("partial/use.csv")
unstack(use, :region, :product, :quantity)

# **Final consumption**, the households' part of it. Firewood needs its own
# source — household energy surveys — because most of it never passes through
# a market that statistics see.

consumption = read_table("partial/consumption.csv")
unstack(consumption, :region, :product, :quantity)

# **Trade**, from customs. Customs see borders, so they report trade between
# the blocs and France's trade **as one country**; they say nothing of which
# French region a log came from or went to, nor of anything moving between
# French regions. Flows under a thousand m³ or tonnes are not reported.

customs = read_table("partial/customs.csv")
first(customs, 8)

# **The mills**: what the EU's industry statistics say each mill that can
# choose its wood actually used, per unit of output, and the pellet plants in
# place, which ran with some room to spare.

input_mix  = read_table("partial/input_mix.csv")
capacities = read_table("partial/capacity.csv")
nothing # hide

# **Freight**, per unit of quantity and per unit of distance, and the
# distances, relative to the Europe–North America crossing:

freight   = read_table("freight.csv")
distances = read_table("distance.csv")

dist = lookup(distances, :distance, :from, :to)
distance(a, b) = a == b ? 0.0 : get(dist, (a, b), get(dist, (b, a), 0.0))
rate = Dict(zip(freight.product, freight.cost))
transport_costs = Dict((p, a, b) => rate[p] * distance(a, b)
                       for p in products, a in regions, b in regions if a != b)
nothing # hide

# **Policies and logistics**, as on the first page: an export duty of a quarter
# of the border price on logs leaving the rest of the world, a ceiling on what
# the Landes can ship, and the base-year price of gas, which wood fuel competes
# with.

log_duty = 0.25
taxes = [Tax(name = "rw_log_duty", on = "trade", product = ["swr", "hwr"], from = "RW",
             rate = log_duty, relative = true)]
duty(p, from) = from == "RW" && p in ("swr", "hwr") ? log_duty : 0.0
trade_limits = [TradeLimit(name = "landes_evacuation", products = ["swr", "smallwood"],
                           from = "SWF", capacity = 18.0)]
gas_reference = 30.0
nothing # hide

# ### Do the data hang together?
#
# Not quite, and it is worth seeing how much before relying on them. Every m³
# produced anywhere is used somewhere, so world production and world use
# should be equal, product by product:

made = lookup(production, :quantity, :region, :product)
used = lookup(use, :quantity, :region, :product)
world(table, p) = sum(q for ((_, q2), q) in table if q2 == p)
DataFrame((product = p, produced = world(made, p), used = world(used, p),
           gap_pct = 100 * (world(used, p) / world(made, p) - 1)) for p in products)

# At most two tenths of a per cent apart: rounding, the usual condition of
# statistics from different sources. Larger gaps would deserve a look before
# going further.

# ## Behaviour: what one year cannot tell
#
# The same elasticities as on the [first page](@ref "The forest sector: France in the world"),
# where each is discussed: supply, demand and its bundles, the mills' choice
# of wood, and the substitution between origins, with the French regions
# near-perfect substitutes for each other.

supply_elasticity = Dict(
    "swr"       => Dict("SEF" => 0.40, "SWF" => 0.45, "GEF" => 0.40, "OF" => 0.40,
                        "EU"  => 0.60, "NA"  => 0.60, "AS"  => 0.55, "RW" => 0.55),
    "hwr"       => Dict("SEF" => 0.35, "SWF" => 0.35, "GEF" => 0.35, "OF" => 0.35,
                        "EU"  => 0.50, "NA"  => 0.50, "AS"  => 0.50, "RW" => 0.50),
    "smallwood" => Dict("SEF" => 0.50, "SWF" => 0.55, "GEF" => 0.50, "OF" => 0.50,
                        "EU"  => 0.60, "NA"  => 0.60, "AS"  => 0.55, "RW" => 0.55))

η_construction, η_fuel, η_paper = 0.45, 0.35, 0.30
σ_lumber, σ_construction, σ_fuel = 5.0, 1.6, 3.0
lumber = ["sawn_sw", "panel", "plywood"]
fuels  = ["smallwood", "pellets"]
gas_cross_elasticity = 0.35

σ_mills = Dict("plymill" => 1.8, "panelmill" => 2.5, "pulpmill" => 1.8, "pelletmill" => 3.0)

σ_blocs = Dict("swr" => 6.0, "hwr" => 6.0, "smallwood" => 7.0, "chips" => 7.0,
               "pulp" => 5.0, "sawn_sw" => 3.5, "sawn_hw" => 3.0, "panel" => 3.0,
               "plywood" => 2.5, "paper" => 4.0, "pellets" => 6.0)
σ_europe = Dict(p => 3σ for (p, σ) in σ_blocs)
σ_france = 30.0
europe(p) = OriginNest(sigma = σ_europe[p], origins = ["EU"],
                       nests = [OriginNest(sigma = σ_france, origins = france)])
nothing # hide

# ## From statistics to model
#
# ### The flows nobody reported
#
# Calibrating trade needs a full matrix of flows between regions, and the
# statistics give part of it. What they give is enough to fill in the rest:
#
# * between the blocs, the flows themselves, from customs;
# * what each bloc keeps of its own production — production less exports;
# * and for the French corner of the matrix, its edges: what each French
#   region produced and used, and what France as a whole bought from and sold
#   to each bloc.
#
# The corner itself — which French region shipped how much to which, and
# which of them the trade with each bloc went through — is filled in by
# [`impute_flows`](@ref): the flows closest to a gravity pattern that add up
# to all those edges. In that small table the French regions stand for
# themselves, and each bloc stands for its trade with France: what it sent
# into France, as an origin, and what it took from France, as a destination.
# Blocs do not trade with each other in it — that trade is already known.

customs_flow = lookup(customs, :quantity, :product, :from, :to)
from_customs(p, o, r) = get(customs_flow, (p, o, r), 0.0)   # unreported: under a thousand

function french_corner(p)
    into_france = Dict(b => from_customs(p, b, "France") for b in blocs)
    out_of_france = Dict(b => from_customs(p, "France", b) for b in blocs)
    return impute_flows(
        production = merge(Dict(r => made[(r, p)] for r in france if haskey(made, (r, p))),
                           into_france),
        use        = merge(Dict(r => used[(r, p)] for r in france if haskey(used, (r, p))),
                           out_of_france),
        distance = distance, decay = 1.4, home_bias = 3.0,
        totals = Dict((blocs, blocs) => 0.0))
end
nothing # hide

# The gravity parameters say how fast trade falls off with distance and how
# much a region leans toward itself. They are guesses, and need not be good
# ones: what the imputed flows must get right is their totals, which they do
# by construction. How little the guess matters is measured in
# [Calibrating with partial data](@ref).
#
# With the corner filled in, the whole matrix:

flows = Dict{Tuple{String,String,String},Float64}()
for p in products
    for ((o, r), x) in french_corner(p)
        x > 0 && (flows[(p, o, r)] = x)
    end
    for o in blocs, r in blocs
        o != r && from_customs(p, o, r) > 0 && (flows[(p, o, r)] = from_customs(p, o, r))
    end
    for b in blocs
        haskey(made, (b, p)) || continue
        exported = sum(from_customs(p, b, r) for r in [blocs; "France"] if r != b)
        flows[(p, b, b)] = made[(b, p)] - exported
    end
end
minimum(values(flows))

# Not negative anywhere, which is the first check a base year has to pass: a
# region cannot export more than it makes. For softwood logs, the French
# corner looks like this:

unstack(DataFrame((from = o, to = r, quantity = x) for ((p, o, r), x) in flows
                  if p == "swr" && (o in france || r in france)), :from, :to, :quantity)

# ### Two prices per region
#
# A variety costs its buyer what its producer was paid, plus the duty if there
# is one, plus the freight to get it there:

πˢ = lookup(producer_prices, :price, :region, :product)   # producer prices
delivered(p, o, r) = o == r ? πˢ[(o, p)] :
                     πˢ[(o, p)] * (1 + duty(p, o)) + transport_costs[(p, o, r)]
nothing # hide

# Each region buys a mix of varieties at different delivered prices, so **two
# prices** coexist in every region: what its producers get, from the
# statistics, and what its users pay for the mix. The second is the average
# delivered price of what the region bought,
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
# way round: imports are cheaper than the local product.
#
# ### Trade shares
#
# A destination's composite is a CES over its origins, and its shares are
# **not** the observed shares: they are the shares each origin would have *if
# all delivered prices were equal*. [`armington_calibration`](@ref) inverts
# the CES demand condition through every level of the origin tree, from the
# flows and the delivered prices. It also sets the `scale` of each composite,
# so that the composite counts the m³ the flows bring in: the m³ the flows
# bring in are then the m³ the mills and the households take out.

armington = [armington_calibration(
                 product = p, sigma = σ_blocs[p], nests = [europe(p)],
                 flows   = Dict((o, r) => x for ((q, o, r), x) in flows if q == p),
                 prices  = Dict((o, r) => delivered(p, o, r)
                                for o in regions, r in regions if haskey(πˢ, (o, p))))
             for p in products]
nothing # hide

# A route that carried nothing in the base year, or too little to be
# reported, gets no share, and stays closed whatever prices do: the
# well-known limit of Armington calibration.
#
# ### Supply
#
# Each forest supply curve runs through the observed harvest at the observed
# producer price:

supply = [SupplySpec(product = p, region = r, p0 = πˢ[(r, p)], q0 = made[(r, p)],
                     elasticity = supply_elasticity[p][r])
          for p in keys(supply_elasticity) for r in regions]
nothing # hide

# ### Demand
#
# Each bundle is calibrated as on the first page, with
# [`ces_calibration`](@ref), but on each region's **consumer** prices and
# consumption rather than on one price for all: the household in `OF` buys the
# mix of French and imported pellets, at what that mix costs in `OF`.

consumed = lookup(consumption, :quantity, :region, :product)

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

function fuel_demand(r, gas = gas_reference)
    fuel = ces_calibration(quantities = [consumed[(r, p)] for p in fuels],
                           prices     = [πᶜ[(r, p)] for p in fuels], sigma = σ_fuel)
    return DemandSpec(region = r, p0 = fuel.price, q0 = fuel.quantity, elasticity = η_fuel,
        shift = exogenous_shift(gas = (level = gas, reference = gas_reference,
                                       elasticity = gas_cross_elasticity)),
        product = DemandNest(sigma = σ_fuel, products = fuels, shares = fuel.shares))
end

build_demand(gas = gas_reference) =
    [[construction_demand(r) for r in regions];
     [fuel_demand(r, gas) for r in regions];
     [DemandSpec(product = "paper", region = r, p0 = πᶜ[(r, "paper")],
                 q0 = consumed[(r, "paper")], elasticity = η_paper) for r in regions]]

demand = build_demand()
nothing # hide

# ### The mills
#
# A mill that can choose its wood is calibrated like a demand bundle, from the
# observed mix at the consumer prices of the region it was observed in. Its
# quantities are all m³ of wood, so here they *do* add up, and the `scale` of
# the calibration lets the bundle keep counting m³: the `composite` is then
# the m³ of wood a unit of output takes, as the conversion factor says it.

function bundle(process)
    mix = input_mix[input_mix.process .== process, :]
    c = ces_calibration(quantities = mix.quantity,
                        prices = [πᶜ[(r, p)] for (r, p) in zip(mix.region, mix.product)],
                        sigma = σ_mills[process])
    return ces(composite = sum(mix.quantity), products = mix.product, shares = c.shares,
               sigma = σ_mills[process], phi = c.scale)
end

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

# The shares are the same wherever the process runs, so a pellet mill in
# `SWF`, facing other prices, uses another mix — the one those prices make
# cheapest. That is why one region's mix was enough.
#
# The `vacost` of a mill — labour, energy, capital — is never observed in the
# units the model wants, but the base year implies it. A mill that ran and did
# not expand earned no more than its costs, so its cost is its **margin**:
# what its outputs fetched at the producer prices of its region, less what its
# inputs cost at the consumer prices there. For a bundle, the cost of the
# inputs is its price index at the region's prices, per m³ of the bundle.

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
# say a mill ran at a loss, and the model would close it. A mill running flat
# out could have earned a rent on top of its costs, which one year cannot tell
# apart from them; that is why the capacities in the data are those of plants
# with room to spare.

capacity = Dict(name => Dict(row.region => row.capacity
                             for row in eachrow(capacities) if row.process == name)
                for name in unique(capacities.process))

processes = [Process(name = name, inputs = inputs, outputs = outputs,
                     regions  = operating(outputs),
                     vacost   = Dict(r => margin(inputs, outputs, r) for r in operating(outputs)),
                     capacity = get(capacity, name, Inf))
             for (name, inputs, outputs) in technology]
nothing # hide

# ## The calibrated reference run

calibrated = MarketData(; regions, products, demand, supply, processes,
                          transport_costs, armington, trade_limits, taxes)
res = solve_market(calibrated)
nothing # hide

# ### How close does it land to the data?
#
# Every curve, share and cost was fitted to one observation, and each passes
# through it exactly. But the observations do not quite agree with each other,
# and an equilibrium has to satisfy every balance at once, so it settles near
# every observation and on none. How near is the first thing to look at, every
# time. Model trade is compared at the level customs report it, with the
# French regions added up:

modelled(table, value, keys...) = lookup(table, value, keys...)
relative_gaps(model, observed) = sort([abs(get(model, k, 0.0) / v - 1) for (k, v) in observed])
at_customs_level(trade) = combine(groupby(transform(trade,
        :from => ByRow(r -> r in france ? "France" : r) => :from,
        :to   => ByRow(r -> r in france ? "France" : r) => :to), [:product, :from, :to]),
    :quantity => sum => :quantity)
reported(trade) = (t = at_customs_level(trade); t[t.from .!= t.to, :])

comparison = [
    "producer prices" => (modelled(res.prices, :producer_price, :region, :product), πˢ),
    "production"      => (modelled(res.production, :quantity, :region, :product), made),
    "consumption"     => (modelled(res.consumption, :quantity, :region, :product), consumed),
    "customs trade"   => (lookup(reported(res.trade), :quantity, :product, :from, :to),
                          lookup(customs, :quantity, :product, :from, :to))]

DataFrame((observed = name, entries = length(observed),
           median_gap = gaps[cld(end, 2)], largest_gap = last(gaps))
          for (name, (model, observed)) in comparison
          for gaps in (relative_gaps(model, observed),))

# Within a fraction of a per cent for prices and quantities, a few per cent at
# worst for the smallest trade flows: gaps of the size of the rounding. The
# data are consistent up to their own precision, and the calibration did its
# job. Gaps of ten or twenty per cent would say the data contradict each other
# — a balance that does not close, a price from another year — or that a step
# went wrong, a unit, a sign, a forgotten duty; either is worth chasing before
# asking the model anything, because every scenario starts from where the
# model landed.
#
# Compare the [first page](@ref "The limits of uncalibrated data"), whose
# reference run missed its own inputs by far more.

# ## The same scenarios
#
# The storm in the Landes and the dearer gas of the first page, unchanged:
# the salvage triples what `SWF` offers at any price, and gas doubles.

reconfigure(d; kwargs...) =
    MarketData(; (f => get(kwargs, f, getfield(d, f)) for f in fieldnames(MarketData))...)
salvage = [s.region == "SWF" && s.product == "swr" ?
               SupplySpec(product = "swr", region = "SWF", p0 = s.p0, q0 = s.q0,
                          elasticity = s.elasticity, shift = 3.0) : s
           for s in supply]
storm = solve_market(reconfigure(calibrated, supply = salvage))
gas   = solve_market(reconfigure(calibrated, demand = build_demand(60.0)))

price_of(r, p, result) = only(result.prices[(result.prices.region .== r) .&
                                            (result.prices.product .== p), :price])
quantity_of(table, r, p) = (rows = table[(table.region .== r) .& (table.product .== p),
                                            :quantity];
                            isempty(rows) ? 0.0 : only(rows))
change(after, before) = 100 * (after / before - 1)

DataFrame((region = r,
           storm_log_price = change(price_of(r, "swr", storm), price_of(r, "swr", res)),
           gas_pellets  = change(quantity_of(gas.consumption, r, "pellets"),
                                 quantity_of(res.consumption, r, "pellets")),
           gas_firewood = change(quantity_of(gas.consumption, r, "smallwood"),
                                 quantity_of(res.consumption, r, "smallwood")))
          for r in regions)

# Percentage changes from the calibrated reference run. They are read the
# same way as on the first page; whether they differ from the first page's,
# and from what perfect data would give, is the subject of the
# [last page](@ref "The forest sector: complete data, and what calibration changes"),
# which reads these results from a file:

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

save_outcomes("partial.csv", ["reference" => res, "storm" => storm, "gas" => gas])
