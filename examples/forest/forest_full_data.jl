# # The forest sector: complete data, and what calibration changes
#
# The [first page](@ref "The forest sector: France in the world") built the
# model from the numbers at hand before any calibration, and the
# [second](@ref "The forest sector: calibrating with partial data") from
# incomplete, rounded statistics. This page builds it from **complete and
# consistent** data — perfect information, which no real application has —
# and then asks the question the three pages exist for: how much do the
# scenarios depend on which of the three you start from?
#
# !!! warning "For illustrative purposes only"
#     The data are plausible, but they are not statistics: complete,
#     consistent data are never observed. They are a base year built to be an
#     equilibrium, so that this page can show what a calibration on perfect
#     information gives.

using Markets, DataFrames, DelimitedFiles

france  = ["SEF", "SWF", "GEF", "OF"]
blocs   = ["EU", "NA", "AS", "RW"]
regions = [france; blocs]
products = ["swr", "hwr", "smallwood", "chips", "pulp",
            "sawn_sw", "sawn_hw", "panel", "plywood", "paper", "pellets"]

data        = pkgdir(Markets, "examples", "forest", "data")
results_dir = pkgdir(Markets, "examples", "forest", "results")

function read_table(path)
    cells, header = readdlm(path, ',', header = true)
    cell(c) = c isa AbstractString ? (isempty(c) ? missing : String(c)) : c
    return DataFrame([String(h) => identity.(cell.(cells[:, j]))
                      for (j, h) in enumerate(vec(header))])
end
lookup(table, value, keys...) =
    Dict(Tuple(row[k] for k in keys) => row[value] for row in eachrow(table))
nothing # hide

# ## Everything there is to know
#
# Complete data differ from the partial ones in one thing above all: **trade**.
# Customs report what crosses a border; complete data give every flow between
# every pair of regions, French regions included. Everything else is the same
# kind of observation, only exact:

production      = read_table(joinpath(data, "full", "production.csv"))
producer_prices = read_table(joinpath(data, "full", "producer_prices.csv"))
consumption     = read_table(joinpath(data, "full", "consumption.csv"))
trade           = read_table(joinpath(data, "full", "trade.csv"))
input_mix       = read_table(joinpath(data, "full", "input_mix.csv"))
capacities      = read_table(joinpath(data, "full", "capacity.csv"))
freight         = read_table(joinpath(data, "freight.csv"))
distances       = read_table(joinpath(data, "distance.csv"))
first(trade[in.(trade.from, Ref(france)) .& in.(trade.to, Ref(france)), :], 6)

# Counted observation by observation, the three pages work from:
#
# | data | first page | second page | this page |
# |:---|---:|---:|---:|
# | harvest and roadside price, by region | 48 | 48 | 48 |
# | output of the mills, by region | — | 62 | 62 |
# | producer prices of the mills' products, by region | — | 62 | 62 |
# | consumer prices | 7, one per product | derived | derived |
# | consumption, by region | 56 | 56 | 56 |
# | use, by region | — | 88 | implied |
# | trade | a gravity rule | 220 customs flows | 595 flows |
# | the mix of wood in the mills | guessed | 14 | 14 |
# | the cost of the mills | guessed | derived | derived |
# | capacities, freight, distances | 47 | 47 | 47 |
#
# The second page is the realistic one; this one is the limit it tends to.

# ## The same calibration, without the gaps
#
# The steps are those of the
# [second page](@ref "The forest sector: calibrating with partial data"),
# explained there, with one fewer: nothing needs imputing. What a region keeps
# of its own production is still not a trade flow, and is still production
# less exports.

dist = lookup(distances, :distance, :from, :to)
distance(a, b) = a == b ? 0.0 : get(dist, (a, b), get(dist, (b, a), 0.0))
rate = Dict(zip(freight.product, freight.cost))
transport_costs = Dict((p, a, b) => rate[p] * distance(a, b)
                       for p in products, a in regions, b in regions if a != b)

log_duty = 0.25
taxes = [Tax(name = "rw_log_duty", on = "trade", product = ["swr", "hwr"], from = "RW",
             rate = log_duty, relative = true)]
duty(p, from) = from == "RW" && p in ("swr", "hwr") ? log_duty : 0.0
evacuation = @limit("landes_evacuation",
                    exports(region = "SWF", product = ["swr", "smallwood"]) <= 18.0)
gas_reference = 30.0

made     = lookup(production, :quantity, :region, :product)
πˢ       = lookup(producer_prices, :price, :region, :product)
consumed = lookup(consumption, :quantity, :region, :product)

flows = lookup(trade, :quantity, :product, :from, :to)
exported(r, p) = sum(x for ((q, o, _), x) in flows if q == p && o == r; init = 0.0)
for ((r, p), q) in made
    flows[(p, r, r)] = q - exported(r, p)
end

delivered(p, o, r) = o == r ? πˢ[(o, p)] :
                     πˢ[(o, p)] * (1 + duty(p, o)) + transport_costs[(p, o, r)]
πᶜ = Dict{Tuple{String,String},Float64}()
for r in regions, p in products
    bought = [(o, x) for ((q, o, d), x) in flows if q == p && d == r && x > 0]
    isempty(bought) && continue
    πᶜ[(r, p)] = sum(x * delivered(p, o, r) for (o, x) in bought) / sum(last, bought)
end
nothing # hide

# The behaviour, chosen as on the other pages:

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

# Trade shares, supply and demand curves:

armington = [armington_calibration(
                 product = p, sigma = σ_blocs[p], nests = [europe(p)],
                 flows   = Dict((o, r) => x for ((q, o, r), x) in flows if q == p),
                 prices  = Dict((o, r) => delivered(p, o, r)
                                for o in regions, r in regions if haskey(πˢ, (o, p))))
             for p in products]

supply = [SupplySpec(product = p, region = r, p0 = πˢ[(r, p)], q0 = made[(r, p)],
                     elasticity = supply_elasticity[p][r])
          for p in keys(supply_elasticity) for r in regions]

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

# and the mills, their bundles from the observed mix and their costs from
# their margins:

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
price_index(δ, p, σ) = sum(δ .* p .^ (1 - σ))^(1 / (1 - σ))
input_cost(n::Nest, r) = length(n.products) == 1 ?
    n.composite * πᶜ[(r, only(n.products))] :
    n.composite * price_index(n.shares, [πᶜ[(r, p)] for p in n.products], n.sigma) / n.phi
margin(inputs, outputs, r) = sum(y * πˢ[(r, k)] for (k, y) in outputs) -
                             sum(input_cost(n, r) for n in inputs)
operating(outputs) = [r for r in regions if get(made, (r, first(first(outputs))), 0.0) > 0]
plant_limits = [@limit("$(row.process) plant $(row.region)",
                       activity(region = row.region, process = row.process) <= row.capacity)
                for row in eachrow(capacities)]
processes = [Process(name = name, inputs = inputs, outputs = outputs,
                     regions  = operating(outputs),
                     vacost   = Dict(r => margin(inputs, outputs, r) for r in operating(outputs)))
             for (name, inputs, outputs) in technology]

complete = MarketData(; regions, products, demand, supply, processes,
                        transport_costs, armington, taxes,
                        limits = [plant_limits; evacuation])
res = solve_market(complete; solver_options = (tol = 1e-10,), tol = 1e-9)
nothing # hide

# The solver's tolerance is tighter than the default, and so is the one below
# which the result tables leave a flow out, only so that the comparison that
# follows is sharp.

# ## The data, returned
#
# Consistent data are an equilibrium of the calibrated model, and a convex
# program has only one: the model must give them back, to the precision of the
# solver. That is *benchmark replication*, the standard test that the data and
# the structure agree:

modelled(table, value, keys...) = lookup(table, value, keys...)
relative_gaps(model, observed) = sort([abs(get(model, k, 0.0) / v - 1) for (k, v) in observed])
shipped = Dict(k => x for (k, x) in flows if k[2] != k[3])

DataFrame((observed = name, entries = length(observed),
           largest_gap = last(relative_gaps(model, observed)))
          for (name, model, observed) in (
              ("producer prices", modelled(res.prices, :producer_price, :region, :product), πˢ),
              ("consumer prices", modelled(res.prices, :price, :region, :product), πᶜ),
              ("production",  modelled(res.production, :quantity, :region, :product), made),
              ("consumption", modelled(res.consumption, :quantity, :region, :product), consumed),
              ("trade",       modelled(res.trade, :quantity, :product, :from, :to), shipped)))

# Every price, quantity and flow comes back. The model knew none of them: it
# was given curves, shares and costs, and found the base year as the point
# where they balance.

# ## Three reference runs, two scenarios
#
# The storm in the Landes and the dearer gas, once more:

reconfigure(d; kwargs...) =
    MarketData(; (f => get(kwargs, f, getfield(d, f)) for f in fieldnames(MarketData))...)
salvage = [s.region == "SWF" && s.product == "swr" ?
               SupplySpec(product = "swr", region = "SWF", p0 = s.p0, q0 = s.q0,
                          elasticity = s.elasticity, shift = 3.0) : s
           for s in supply]
storm = solve_market(reconfigure(complete, supply = salvage))
gas   = solve_market(reconfigure(complete, demand = build_demand(60.0)))
nothing # hide

# The other two pages saved their reference runs and scenarios. Put side by
# side with this one, in the same shape:

function outcome(result, scenario)
    t = outerjoin(select(result.prices, :region, :product, :price, :producer_price),
                  rename(result.production, :quantity => :production),
                  select(result.consumption, :region, :product, :quantity => :consumption),
                  on = [:region, :product])
    return insertcols!(t, 1, :scenario => scenario)
end

cases = ["uncalibrated (first page)"  => read_table(joinpath(results_dir, "uncalibrated.csv")),
         "partial data (second page)" => read_table(joinpath(results_dir, "partial.csv")),
         "complete data (this page)"  => vcat(outcome(res, "reference"), outcome(storm, "storm"),
                                              outcome(gas, "gas"))]
nothing # hide

# For each, how far its reference run sits from the complete data, and what
# the two scenarios do to it:

value(t, scenario, r, p, column) =
    (rows = t[(t.scenario .== scenario) .& (t.region .== r) .& (t.product .== p), column];
     isempty(rows) ? missing : only(rows))
change(t, scenario, r, p, column) =
    100 * (value(t, scenario, r, p, column) / value(t, "reference", r, p, column) - 1)
median_and_largest(gaps) = (round(100 * gaps[cld(end, 2)]; digits = 1),
                            round(100 * last(gaps); digits = 1))
span(xs) = string(round(minimum(xs); digits = 1), " to ", round(maximum(xs); digits = 1))

comparison = DataFrame(
    (case = name,
     producer_prices_gap = median_and_largest(relative_gaps(
         Dict((r, p) => value(t, "reference", r, p, :producer_price) for (r, p) in keys(πˢ)), πˢ)),
     production_gap = median_and_largest(relative_gaps(
         Dict((r, p) => coalesce(value(t, "reference", r, p, :production), 0.0)
              for (r, p) in keys(made)), made)),
     storm_SWF = round(change(t, "storm", "SWF", "swr", :price); digits = 1),
     storm_elsewhere = span(change(t, "storm", r, "swr", :price) for r in regions if r != "SWF"),
     gas_pellets  = span(change(t, "gas", r, "pellets", :consumption) for r in regions),
     gas_firewood = span(change(t, "gas", r, "smallwood", :consumption) for r in regions))
    for (name, t) in cases)

# The gaps are percentages, median and largest; the scenarios are percentage
# changes from each case's own reference run: the log price in `SWF` and in the
# other regions under the storm, and the pellets and firewood consumed under
# the dearer gas.
#
# The partial data give practically the results of perfect information. Their
# reference run sits within a per cent of the complete data — the rounding —
# and the scenarios agree to a tenth or two of a point, although every flow
# touching a French region was imputed rather than observed. Detail that is
# missing costs little as long as the totals around it are right.
#
# The uncalibrated reference run is another matter. Its prices are 5% off in
# the median and a quarter at worst, and the output of some mills is off by
# double — the French ones, whose costs and recipes were guessed. Where a
# scenario runs into those errors, its result changes with them:
#
# * the storm takes the softwood price in `SWF` down by 12% instead of 5%,
#   because the uncalibrated run already ships about 13 of the 18 Mm³ the
#   Landes roads allow, against about 10 in the data, so the limit binds
#   sooner and the salvage has nowhere to go;
# * the dearer gas moves pellets and firewood by different amounts, because
#   the guessed mills have different margins and different room in their
#   plants.
#
# Where the scenario runs into nothing of the kind, the uncalibrated run gets
# it about right: the storm's effect outside `SWF` agrees within a few tenths
# of a point, because constant-elasticity curves give the same percentage
# response wherever they are evaluated.
#
# So calibration matters, but not everywhere equally, and perfect data are not
# what it takes: a reference run in the right place — on the right totals, at
# the right prices, at the right distance from the limits that matter — is
# enough to get the scenarios right.
#
# The general lessons — what a scenario is robust to, what it is not, and
# where to spend the effort on data — are drawn in
# [Calibrating with partial data](@ref).
