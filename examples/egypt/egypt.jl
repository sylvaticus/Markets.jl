# # Egyptian agriculture, 1980
#
# In the late 1970s the government of Egypt and the UNDP drew up a master plan
# for the country's water: how much the Nile could supply, month by month, and
# what it would be worth to invest in more. At its centre was an agricultural
# sector model by Gerald Kutcher, published with the World Bank's sector
# models and kept, in a compact version, as
# [`egypt`](https://www.gams.com/latest/gamslib_ml/libhtml/gamslib_egypt.html)
# in the GAMS model library.
#
# It belongs to the same family as this package: farmers choose what to grow,
# region by region, and the national market sets prices where the surplus of
# consumers and producers together is largest. This page ports it, checks the
# port against a line-by-line transcription of the original, and then asks it
# the question it was built for.
#
# Kutcher, G. P. (1980), *The Agro-Economic Model*, Master Plan for Water
# Resources Development, UNDP-EGY-73/024.

using Markets, DataFrames, DelimitedFiles

# ## The data
#
# The GAMS tables are in `examples/egypt/data`, one CSV file each, in long
# format: one value per row, keyed as in the original.

data = pkgdir(Markets, "examples", "egypt", "data")
function read_table(file)
    cells, header = readdlm(joinpath(data, file), ',', header = true)
    text(c) = c isa AbstractString ? String(c) : c
    return DataFrame([String(h) => identity.(text.(cells[:, j])) for (j, h) in enumerate(vec(header))])
end
## a table as a dictionary from its keys to its value, which is the last column
function table(file)
    t = read_table(file)
    keys = names(t)[1:end-1]
    key(row) = length(keys) == 1 ? row[keys[1]] : Tuple(row[k] for k in keys)
    return Dict(key(row) => float(row[end]) for row in eachrow(t))
end
nothing # hide

# Five regions, named after the Nile barrages that water them, in three
# agricultural zones; the farming calendar starts in October.

farms  = ["u-egypt", "m-egypt", "e-delta", "m-delta", "w-delta"]
zone   = Dict("u-egypt" => "upper", "m-egypt" => "middle", "e-delta" => "delta",
              "m-delta" => "delta", "w-delta" => "delta")
months = ["oct", "nov", "dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep"]
nothing # hide

# The crops: yields by region, per feddan (0.42 ha); land, labour and water
# needed in each month, by zone; fodder left for the animals, by zone; and the
# rest of what a feddan costs.

yields   = table("yield.csv")       # (crop, region) => t/feddan
indprod  = table("indprod.csv")     # (crop, region) => t/feddan of an industrial by-product
feed     = table("feed.csv")        # (crop, zone)   => t/feddan of fodder
land     = table("land.csv")        # (crop, zone, month) => share of the month the field is taken
laborm   = table("laborm.csv")      # (crop, zone, month) => man-days per feddan
water    = table("water.csv")       # (crop, zone, month) => water per feddan
tractor  = table("tractor.csv")     # (crop, zone)   => tractor hours per feddan
cropdat  = table("cropdat.csv")     # (crop, item)   => seed, inputs, protein and starch content
nothing # hide

# The regions: their area (1000 feddans), farmers (thousands), and the protein
# and starch their herds need (1000 t); monthly wages; the area that may be
# under vegetables in each season; upper bounds on cotton, citrus and
# sugarcane; and what it costs to carry straw between three of the regions.

agrreg         = table("agrreg.csv")    # (region, item)
wages          = table("regwagm.csv")   # (region, month) => £E per day
vegetable_area = table("veg.csv")       # (region, season)
upbnds         = table("upbnds.csv")    # (crop, region)
transfer_cost  = Dict(k => v for (k, v) in table("tranc.csv") if v > 0)
nothing # hide

# The national market: base price, consumption and price elasticity of every
# marketed crop, the prices at which it can be imported and exported, and how
# much may be.

demdat       = table("demdat.csv")      # (crop, item)
trade_limits = table("tlimit.csv")      # (crop, "imports"/"exports") => 1000 t
marketed     = unique(first.(collect(keys(demdat))))
nothing # hide

# The scalars and small sets are in the GAMS code itself. Some crops are sold
# under another name (summer and nili maize are both maize); cotton, soybeans
# and sesame also yield vegetable oil; part of what is sold is lost in storage.

reclassified = Dict("s-maize" => "maize", "n-maize" => "maize", "s-tomato" => "s+n-tomato",
                    "n-tomato" => "s+n-tomato", "s-potato" => "potato", "n-potato" => "potato")
oil          = Dict("cotton" => 0.175, "soyabeans" => 0.175, "sesame" => 0.470)
losses       = Dict("wheat" => 0.090, "barley" => 0.114, "horse-bean" => 0.152,
                    "lentils" => 0.188, "w-tomato" => 0.100, "cotton" => 0.540,
                    "s-maize" => 0.049, "s-tomato" => 0.100, "sorghum" => 0.055,
                    "w-onion" => 0.057, "rice" => 0.045, "groundnut" => 0.038,
                    "s-potato" => 0.285, "citrus" => 0.100)
seasons      = ["winter", "summer", "nili"]
vegetables   = Dict("winter" => ["w-veget", "w-tomato", "w-onion"],
                    "summer" => ["s-veget", "s-tomato", "s-potato"],
                    "nili"   => ["n-veget", "n-tomato", "n-potato"])
green        = ["lo-berseem", "sh-berseem", "sugarcane", "s-maize", "n-maize", "sorghum"]
dry          = ["wheat", "barley", "rice"]
transported  = ["wheat", "barley"]      # whose straw can be carried between regions
byproduct_price = Dict("horse-bean" => 21.75, "rice" => 5.75)
working_days    = 20      # per month
nitrogen_price  = 178     # £E per t
phosphate_price = 200     # £E per t
tractor_hour    = 1.5     # £E
concentrate_price = 1000  # £E per t of protein or starch bought in
green_dry_ratio = 0.06723
meat_and_milk   = 749576  # a constant in the objective, the value of livestock output
nothing # hide

# (Three of the losses — summer maize, tomato and potato — are on names that
# are never sold, so in the original they never apply. The port keeps that.)
#
# ## From GAMS to Markets.jl
#
# The original is a linear program in its own terms: cropping activities,
# balances, and an objective that adds up surplus. Each piece has a
# counterpart here.
#
# | `egypt.gms`                              | here                                                      |
# |:-----------------------------------------|:----------------------------------------------------------|
# | `xcrop(r,c)`, area under a crop          | a [`Process`](@ref) per crop and region, run at that level |
# | `yld`, `feed`, losses in storage         | the outputs of the process, net of the losses             |
# | `netcs(c,r)`, inputs other than labour   | its value-added cost                                      |
# | `labbal`, `flab`, `tlab`                 | labour of each month as a product, made by family and hired labour |
# | `landb`, `waterb`, `vegetb`, `xcrop.up`  | [`Limit`](@ref)s on the activities                        |
# | `nutb`, `anut`                           | fodder fed to a herd whose size is fixed, or concentrates bought in |
# | `gfodb`                                  | a limit on the fodder produced and carried                |
# | demand on a grid of ten points           | a linear [`DemandCurve`](@ref) through the same base point |
# | `pe`, `pm`, `tlimit`                     | a flat world demand and supply, limits on the flows       |
# | `trans`, `stran`                         | trade in straw at its transfer cost                       |
#
# Two regions are added to the five: `"egypt"`, the national market where all
# the crops are sold and imports arrive, and `"world"`.

regions  = [farms; "egypt"; "world"]
fodders  = sort(unique(first.(collect(keys(feed)))))
labour   = ["labour $tm" for tm in months]
products = [marketed; ["$c fodder" for c in fodders]; "protein"; "starch"; labour]
nothing # hide

# **Crops.** A crop is grown in a region if it yields something there, a
# product for the market or fodder. Its outputs per feddan are what is sold,
# net of storage losses, and its fodder; its inputs the labour of each month;
# its cost the tractor hours, seed, fertiliser and pesticide, less the value of
# its industrial by-product.

function crop_outputs(c, r)
    out = Pair{String,Float64}[]
    y = get(yields, (c, r), 0.0)
    if y > 0
        sold = get(reclassified, c, c)
        push!(out, sold => y * (1 - get(losses, sold, 0.0)))
        haskey(oil, c) && push!(out, "veg-oil" => y * oil[c])
    end
    f = get(feed, (c, zone[r]), 0.0)
    f > 0 && push!(out, "$c fodder" => f)
    return out
end
crop_cost(c, r) = get(tractor, (c, zone[r]), 0.0) * tractor_hour +
    get(cropdat, (c, "misc"), 0.0) + get(cropdat, (c, "pest"), 0.0) +
    get(cropdat, (c, "n-fer"), 0.0) * nitrogen_price +
    get(cropdat, (c, "p-fer"), 0.0) * phosphate_price -
    get(indprod, (c, r), 0.0) * get(byproduct_price, c, 0.0)
crops = sort(unique(first.(collect(keys(land)))))
grown = [(c, r) for r in farms for c in crops if !isempty(crop_outputs(c, r))]
crop_name(c, r) = "$c in $r"
crop_processes = [Process(name = crop_name(c, r), regions = [r],
                          inputs = [leontief(product = "labour $tm",
                                             coeff = laborm[(c, zone[r], tm)] / working_days)
                                    for tm in months if get(laborm, (c, zone[r], tm), 0.0) > 0],
                          outputs = crop_outputs(c, r),
                          ## the original charges these only where the crop has a yield
                          vacost = get(yields, (c, r), 0.0) > 0 ? crop_cost(c, r) : 0.0)
                  for (c, r) in grown]
crop_processes[1]

# **Labour.** Family labour costs the regional wage; hired labour twice as
# much, and there is no limit to it. Labour is measured in thousands of
# man-months, so a man-day costs a twentieth of it.

labour_processes = [Process(name = "$who labour $tm", inputs = Nest[], outputs = ["labour $tm" => 1.0],
                            regions = farms,
                            vacost = Dict(r => factor * wages[(r, tm)] * working_days for r in farms))
                    for (who, factor) in (("family", 1), ("hired", 2)) for tm in months]
nothing # hide

# **Animals.** The herds of each region need a given amount of protein and
# starch. It comes from the fodder, at its protein and starch content, or from
# concentrates bought at £E 1000 a tonne; what the fodder offers beyond the
# need is wasted. The herd is a process whose activity a limit below holds at
# one.

feeding = [Process(name = "feed $c", inputs = [leontief(product = "$c fodder", coeff = 1.0)],
                   outputs = [n => cropdat[(c, n)] / 100 for n in ("protein", "starch")
                              if get(cropdat, (c, n), 0.0) > 0])
           for c in fodders]
concentrates = [Process(name = "concentrate $n", inputs = Nest[], outputs = [n => 1.0],
                        vacost = concentrate_price, regions = farms) for n in ("protein", "starch")]
herds = [Process(name = "herd $r", regions = [r], outputs = Pair{String,Float64}[],
                 inputs = [leontief(product = n, coeff = agrreg[(r, n)]) for n in ("protein", "starch")])
         for r in farms]
nothing # hide

# **What is not used is thrown away.** The original lets production exceed
# sales, and nutrients exceed the herds' needs; in a market, that is free
# disposal.

disposals = [[disposal(product = c, regions = ["egypt"]) for c in marketed];
             [disposal(product = n, regions = farms) for n in ("protein", "starch")]]
processes = [crop_processes; labour_processes; feeding; concentrates; herds; disposals]
nothing # hide

# **The national market.** Demand is a straight line through the base price and
# consumption, with the base elasticity, exactly as the original computes its
# coefficients. The original never lets the price leave the range from 0.6 to
# 1.666 times the base price, which bounds consumption: a limit below keeps
# that range.

demand = [DemandCurve(product = c, region = "egypt", p0 = demdat[(c, "base-p")],
                     q0 = demdat[(c, "domes-cons")], elasticity = -demdat[(c, "elas")],
                     form = :linear)
          for c in marketed]
demand[1]

# **The world** sells at the import price and buys at the export price, in any
# quantity, for the crops the original lets be imported or exported.

importable = [c for c in marketed if get(trade_limits, (c, "imports"), 0.0) > 0]
exportable = [c for c in marketed if get(trade_limits, (c, "exports"), 0.0) > 0]
supply = [SupplyCurve(product = c, region = "world", p0 = demdat[(c, "import-p")], q0 = 1,
                     elasticity = Inf) for c in importable]
append!(demand, [DemandCurve(product = c, region = "world", p0 = demdat[(c, "export-p")], q0 = 1,
                            elasticity = Inf) for c in exportable])
nothing # hide

# **Routes.** Crops go from the farms to the national market at no cost, as in
# the original, which has a single national balance; imports come in and
# exports go out through it; straw moves between three regions at its
# transfer cost.

transport_costs = Dict{Tuple{String,String,String},Float64}()
for c in marketed, r in farms
    transport_costs[(c, r, "egypt")] = 0.0
end
for c in importable
    transport_costs[(c, "world", "egypt")] = 0.0
end
for c in exportable
    transport_costs[(c, "egypt", "world")] = 0.0
end
for c in transported, ((r, rp), cost) in transfer_cost
    transport_costs[("$c fodder", r, rp)] = cost
end
nothing # hide

# ## The constraints, as limits
#
# Everything else in the original is a constraint on the activities: land in
# each month, water over the year, the area under vegetables in each season,
# the ceilings on cotton, citrus and sugarcane. Each is a family of
# [`Limit`](@ref)s, written as the original writes it.

grows(c, r) = (c, r) in grown
land_limits = @limits "land" [r in farms, tm in months] sum(
    land[(c, zone[r], tm)] * activity(region = r, process = crop_name(c, r))
    for c in crops if grows(c, r) && haskey(land, (c, zone[r], tm))) <= agrreg[(r, "area")]
water_limit = @limit "water" sum(
    water[(c, zone[r], tm)] * activity(region = r, process = crop_name(c, r))
    for (c, r) in grown for tm in months if haskey(water, (c, zone[r], tm))) <= 30000 * 1000
vegetable_limits = @limits "vegetables" [s in seasons, r in farms] sum(
    activity(region = r, process = crop_name(c, r))
    for c in vegetables[s] if grows(c, r)) <= 1.5 * vegetable_area[(r, s)]
area_limits = [@limit("area[$c, $r]", activity(region = r, process = crop_name(c, r)) <= a)
               for ((c, r), a) in upbnds if grows(c, r)]
nothing # hide

# Family labour is limited to 40% of the farmers in any month, and each herd is
# fed exactly once.

family_limits = @limits "family labour" [r in farms, tm in months] activity(
    region = r, process = "family labour $tm") <= 0.4 * agrreg[(r, "farmers")]
herd_limits = @limits "herd" [r in farms] activity(region = r, process = "herd $r") == 1
nothing # hide

# Green fodder must not outweigh dry fodder by more than a given ratio, where
# dry fodder counts the straw carried in and not the straw carried out. A
# limit can name the production of a product, and its imports and exports,
# only where they exist, hence the conditions.

makes(r, p) = any(proc -> (proc.regions isa Colon || r in proc.regions) &&
                          any(first(o) == p for o in proc.outputs), processes)
moves(r, p) = any(k -> k[1] == p && (k[2] == r || k[3] == r), keys(transport_costs))
fodder_limits = @limits "green fodder" [r in farms] green_dry_ratio * sum(
        production(region = r, product = "$c fodder") for c in green if makes(r, "$c fodder")) -
    sum(production(region = r, product = "$c fodder") for c in dry if makes(r, "$c fodder")) -
    sum(imports(region = r, product = "$c fodder") - exports(region = r, product = "$c fodder")
        for c in transported if moves(r, "$c fodder"); init = zero(LimitExpression)) <= 0
nothing # hide

# What can be traded with the world is limited for some crops, and consumption
# stays where the original's grid of prices keeps it.

trade_caps = [@limit("$dir[$c]", flow(product = c, from = from, to = to) <= cap)
              for ((c, dir), cap) in trade_limits if isfinite(cap) && cap > 0
              for (from, to) in (dir == "imports" ? ("world", "egypt") : ("egypt", "world"),)]
slope(c)  = demdat[(c, "base-p")] / demdat[(c, "domes-cons")] / demdat[(c, "elas")]
choke(c)  = demdat[(c, "base-p")] - slope(c) * demdat[(c, "domes-cons")]
at_price(c, factor) = (factor * demdat[(c, "base-p")] - choke(c)) / slope(c)
consumption_range = @limits "consumption" [c in marketed] max(0.0, at_price(c, 1.666)) <=
    consumption(region = "egypt", product = c) <= at_price(c, 0.6)
limits = [land_limits; water_limit; vegetable_limits; area_limits; family_limits; herd_limits;
          fodder_limits; trade_caps; consumption_range]
length(limits)

# ## The equilibrium
#
# The solver is asked for a relative precision of 10⁻⁷ rather than its usual
# 10⁻⁸: the problem is nearly an LP, with many crops tied at the margin, and
# depending on the platform the last digit can cost Ipopt more than it is
# worth.

egypt = MarketData(regions = regions, products = products, demand = demand, supply = supply,
                   processes = processes, limits = limits, transport_costs = transport_costs)
result = solve_market(egypt; solver_options = (tol = 1e-7,))
result.status

# What is grown where, in thousands of feddans:

crop_of = Dict(crop_name(c, r) => c for (c, r) in grown)
areas = filter(row -> haskey(crop_of, row.process), result.activity)
areas.crop = [crop_of[p] for p in areas.process]
area_table = unstack(select(areas, :crop, :region, :level), :region, :level; fill = 0.0)
sort!(area_table, :crop)

# The national market: price, consumption, and what crosses the border.

flow_of(c, from, to) = (rows = result.trade[(result.trade.product .== c) .& (result.trade.from .== from) .&
                                            (result.trade.to .== to), :quantity];
                        isempty(rows) ? 0.0 : only(rows))
national = DataFrame(crop = marketed,
    base_price = [demdat[(c, "base-p")] for c in marketed],
    price = [only(result.prices[(result.prices.region .== "egypt") .& (result.prices.product .== c), :price])
             for c in marketed],
    base_consumption = [demdat[(c, "domes-cons")] for c in marketed],
    consumption = [only(result.consumption[(result.consumption.region .== "egypt") .&
                                           (result.consumption.product .== c), :quantity])
                   for c in marketed],
    imports = [flow_of(c, "world", "egypt") for c in marketed],
    exports = [flow_of(c, "egypt", "world") for c in marketed])

# And what the scarce resources are worth: the rent of each limit that binds,
# in £E per unit. Land is scarce at the turn of the seasons, in November and in
# May, when the crops of one season are still in the ground and those of the
# next want it; family labour at the peaks of work; cotton, citrus and
# sugarcane at their ceilings; and several exports at their caps. Water is
# not: the crops the market asks for need 28.9 of the 30 million units
# available.

binding = filter(row -> row.rent > 1e-6, result.limits)
select(binding, :limit, :value, :rent)

# ## Checking the port
#
# [`egypt_gams.jl`](https://github.com/sylvaticus/Markets.jl/blob/main/examples/egypt/egypt_gams.jl)
# transcribes the original line by line into JuMP: the same variables, the
# same equations, the same objective, on the data loaded above. With
# `grid = false` its demand curves are the exact lines rather than ten points
# on each, which is the only change, and which the port shares.

include(joinpath(pkgdir(Markets), "examples", "egypt", "egypt_gams.jl"))
exact = egypt_gams(grid = false)
(transcription = exact.objective - meat_and_milk, port = result.objective,
 difference = exact.objective - meat_and_milk - result.objective)

# The objectives agree, and so do the prices:

comparison = DataFrame(crop = marketed, transcription = [exact.price[c] for c in marketed],
                       port = national.price)
comparison.difference = comparison.port .- comparison.transcription
comparison

# The original itself solves the grid. Ten points per curve, joined by chords,
# put welfare a little below the curve between the points, so the grid's
# solution sits near the exact one but not on it:

original = egypt_gams(grid = true)
DataFrame(crop = marketed, original = [original.price[c] for c in marketed],
          exact = [exact.price[c] for c in marketed])

# ## The question it was built for
#
# What is water worth, and what happens if there is less of it? Take a fifth of
# it away, solving to the same precision as before:

drier = MarketData(regions = regions, products = products, demand = demand, supply = supply,
                   processes = processes, transport_costs = transport_costs,
                   limits = [l.name == "water" ? @limit("water", sum(
                                 water[(c, zone[r], tm)] * activity(region = r, process = crop_name(c, r))
                                 for (c, r) in grown for tm in months
                                 if haskey(water, (c, zone[r], tm))) <= 0.8 * 30000 * 1000) : l
                             for l in limits])
dry_result = solve_market(drier; solver_options = (tol = 1e-7,))
water_rent(res) = only(res.limits[res.limits.limit .== "water", :rent])
(water_rent = water_rent(result), with_a_fifth_less = water_rent(dry_result),
 welfare_lost = result.objective - dry_result.objective)

# Water now binds, and earns a rent; the transcription, given the same
# water, agrees on it and on the welfare:

exact_drier = egypt_gams(grid = false, water_available = 0.8 * 30000 * 1000)
(transcription = (water_rent = exact_drier.water_rent,
                  objective = exact_drier.objective - meat_and_milk),
 port = (water_rent = water_rent(dry_result), objective = dry_result.objective))

# The crops that give way are the thirsty summer ones, rice above all; the
# winter crops, wheat and berseem, which need least water, take back some of
# the land:

total_area(res) = (a = filter(row -> haskey(crop_of, row.process), res.activity);
                   combine(groupby(transform(a, :process => ByRow(p -> crop_of[p]) => :crop), :crop),
                           :level => sum => :area))
change = outerjoin(rename(total_area(result), :area => :now), rename(total_area(dry_result), :area => :drier),
                   on = :crop)
change.now = coalesce.(change.now, 0.0); change.drier = coalesce.(change.drier, 0.0)
change.difference = change.drier .- change.now
sort!(filter(row -> abs(row.difference) > 1, change), :difference)
