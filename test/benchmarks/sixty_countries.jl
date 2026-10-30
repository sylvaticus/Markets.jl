# Part of Markets.jl package [https://github.com/sylvaticus/Markets.jl] - MIT Licence

# ----------------------------------------------------------------------------
#  Benchmark: sixty countries and an ad valorem tariff on every route
#
#  The forest example at scale: its 11 products and 7 processes, its demand
#  groups and its Armington trade, in 60 countries on 6 continents, each
#  continent a group of origins as Europe is in the example. Sizes and
#  positions are drawn at random, with a fixed seed, and scaled to the
#  example's world totals.
#
#  It reports how the solver copes with the untaxed economy, then with a 20%
#  ad valorem tariff on all 38,940 routes, which `solve_market` finds by
#  repeated solves: how many, how fast they converge, and what they cost.
#  `--cold` adds the same loop with every solve started from scratch, for
#  comparison. It is not part of `Pkg.test()`: it takes a few minutes.
#
#      julia --project=docs test/benchmarks/sixty_countries.jl [--cold]
# ----------------------------------------------------------------------------

using Markets, DataFrames, Random, Printf, Logging

# each solve of the tariff loop logs its gap, iterations and status
ENV["JULIA_DEBUG"] = "Markets"

rng = MersenneTwister(2026)
n_continents, per_continent = 6, 10
countries = [@sprintf("C%02d", i) for i in 1:n_continents*per_continent]
continent_of = Dict(c => (i - 1) ÷ per_continent + 1 for (i, c) in enumerate(countries))
continents = [[c for c in countries if continent_of[c] == k] for k in 1:n_continents]

# positions: continents on a circle about 1 apart, countries scattered around
# their centre, so that distances span the forest example's range (0.07–1.35)
centre = [(cos(2π * k / n_continents), sin(2π * k / n_continents)) .* 0.6 for k in 1:n_continents]
pos = Dict(c => centre[continent_of[c]] .+ 0.12 .* (randn(rng), randn(rng)) for c in countries)
distance(a, b) = a == b ? 0.0 : hypot((pos[a] .- pos[b])...)

# sizes: where the forests are and where the people are are related, not equal
forest_w = exp.(randn(rng, length(countries)))
people_w = 0.5 .* forest_w ./ sum(forest_w) .+ 0.5 .* (w = exp.(randn(rng, length(countries))); w ./ sum(w))
forest_w ./= sum(forest_w)
fw = Dict(zip(countries, forest_w)); pw = Dict(zip(countries, people_w))

products = ["swr", "hwr", "smallwood", "chips", "pulp",
            "sawn_sw", "sawn_hw", "panel", "plywood", "paper", "pellets"]

# world totals of the forest example, shared out
supply_total = Dict("swr" => 1302.2, "hwr" => 617.6, "smallwood" => 1172.3)
supply_price = Dict("swr" => 100.0, "hwr" => 118.0, "smallwood" => 67.0)
supply = [SupplyCurve(product = p, region = c,
                      p0 = supply_price[p] * (1 + 0.12 * randn(rng)),
                      q0 = supply_total[p] * fw[c] * exp(0.3 * randn(rng)),
                      elasticity = 0.4 + 0.2 * rand(rng))
          for p in ("swr", "hwr", "smallwood") for c in countries]

ref_price = Dict("sawn_sw" => 205.0, "sawn_hw" => 280.0, "panel" => 240.0, "plywood" => 285.0,
                 "paper" => 935.0, "smallwood" => 70.0, "pellets" => 285.0)
ref_total = Dict("sawn_sw" => 309.6, "sawn_hw" => 95.9, "panel" => 236.2, "plywood" => 143.5,
                 "paper" => 335.4, "smallwood" => 350.5, "pellets" => 70.5)
q(p, c) = ref_total[p] * pw[c]

σ_lumber, σ_construction, σ_fuel = 5.0, 1.6, 3.0
lumber = ["sawn_sw", "panel", "plywood"]
fuels = ["smallwood", "pellets"]
function construction_demand(r)
    inner = ces_calibration(quantities = [q(p, r) for p in lumber],
                            prices = [ref_price[p] for p in lumber], sigma = σ_lumber)
    outer = ces_calibration(quantities = [q("sawn_hw", r), inner.quantity],
                            prices = [ref_price["sawn_hw"], inner.price], sigma = σ_construction)
    DemandCurve(region = r, p0 = outer.price, q0 = outer.quantity, elasticity = 0.45,
        product = DemandNest(sigma = σ_construction, products = ["sawn_hw"], shares = [outer.shares[1]],
            nests = [DemandNest(sigma = σ_lumber, share = outer.shares[2],
                                products = lumber, shares = inner.shares)]))
end
function fuel_demand(r)
    fuel = ces_calibration(quantities = [q(p, r) for p in fuels],
                           prices = [ref_price[p] for p in fuels], sigma = σ_fuel)
    DemandCurve(region = r, p0 = fuel.price, q0 = fuel.quantity, elasticity = 0.35,
                product = DemandNest(sigma = σ_fuel, products = fuels, shares = fuel.shares))
end
demand = [[construction_demand(r) for r in countries];
          [fuel_demand(r) for r in countries];
          [DemandCurve(product = "paper", region = r, p0 = ref_price["paper"],
                       q0 = q("paper", r), elasticity = 0.30) for r in countries]]

pulp_regions = [c for c in countries if rand(rng) < 0.6]
processes = [
    Process(name = "sawmill_sw", inputs = [leontief(product = "swr", coeff = 1.0)],
            outputs = ["sawn_sw" => 0.50, "chips" => 0.35], vacost = 40),
    Process(name = "sawmill_hw", inputs = [leontief(product = "hwr", coeff = 1.0)],
            outputs = ["sawn_hw" => 0.45, "chips" => 0.35], vacost = 45),
    Process(name = "plymill",
            inputs = [ces(composite = 2.0, products = ["swr", "hwr"], shares = [0.45, 0.55], sigma = 1.8)],
            outputs = ["plywood" => 1.0, "chips" => 0.8], vacost = 180),
    Process(name = "panelmill",
            inputs = [ces(composite = 1.30, products = ["chips", "smallwood", "swr", "hwr"],
                          shares = [0.40, 0.25, 0.23, 0.12], sigma = 2.5)],
            outputs = ["panel" => 1.0], vacost = 120),
    Process(name = "pulpmill", regions = pulp_regions,
            inputs = [ces(composite = 4.0, products = ["swr", "hwr", "chips", "smallwood"],
                          shares = [0.35, 0.22, 0.23, 0.20], sigma = 1.8)],
            outputs = ["pulp" => 1.0], vacost = 300),
    Process(name = "papermill", inputs = [leontief(product = "pulp", coeff = 1.10)],
            outputs = ["paper" => 1.0], vacost = 200),
    Process(name = "pelletmill",
            inputs = [ces(composite = 2.2, products = ["chips", "smallwood", "swr", "hwr"],
                          shares = [0.40, 0.40, 0.12, 0.08], sigma = 3.0)],
            outputs = ["pellets" => 1.0], vacost = 90),
]

freight = Dict("swr" => 12.0, "hwr" => 12.0, "smallwood" => 13.0, "chips" => 14.0,
               "pulp" => 35.0, "sawn_sw" => 30.0, "sawn_hw" => 32.0, "panel" => 35.0,
               "plywood" => 38.0, "paper" => 55.0, "pellets" => 28.0)
transport_costs = Dict((p, a, b) => freight[p] * distance(a, b)
                       for p in products, a in countries, b in countries if a != b)

σ_blocs = Dict("swr" => 6.0, "hwr" => 6.0, "smallwood" => 7.0, "chips" => 7.0,
               "pulp" => 5.0, "sawn_sw" => 3.5, "sawn_hw" => 3.0, "panel" => 3.0,
               "plywood" => 2.5, "paper" => 4.0, "pellets" => 6.0)
size_of = Dict(c => sum(s.q0 for s in supply if s.region == c) for c in countries)
gravity = Dict((o, r) => size_of[o] * exp(-1.4 * distance(o, r)) * (o == r ? 3.0 : 1.0)
               for o in countries, r in countries)
offers(p) = p == "pulp" ? pulp_regions : countries
armington = [Armington(product = p, sigma = σ_blocs[p],
                       shares = Dict(k => v for (k, v) in gravity if k[1] in offers(p)),
                       nests = [OriginNest(sigma = 3σ_blocs[p], origins = filter(in(offers(p)), cs))
                                for cs in continents])
             for p in products]

market(taxes) = MarketData(regions = countries, products = products, demand = demand,
                           supply = supply, processes = processes,
                           transport_costs = transport_costs, armington = armington, taxes = taxes)
tariff20 = [Tax(name = "tariff", on = "trade", rate = 0.2, relative = true)]


println("economy: ", repr(market(tariff20)), ", ", length(transport_costs), " routes")

# A first solve compiles Markets, JuMP and Ipopt, which takes longer than the
# solve itself: a small economy built from the same pieces — a demand group, a
# CES process, Armington trade and an ad valorem tariff — pays for it here, so
# that the timings below are of solving alone and can be compared.
warmup = MarketData(regions = ["A", "B"], products = ["w", "g", "h"],
    demand = [DemandCurve(product = DemandNest(sigma = 2.0, products = ["g", "h"]),
                          region = r, p0 = 100, q0 = 10, elasticity = 0.5) for r in ("A", "B")],
    supply = [SupplyCurve(product = p, region = r, p0 = 30, q0 = r == "A" ? 40 : 5,
                          elasticity = 1.0) for p in ("w", "h") for r in ("A", "B")],
    processes = [Process(name = "mill", outputs = ["g" => 1.0], vacost = 10,
                         inputs = [ces(composite = 1.5, products = ["w", "h"],
                                       shares = [0.7, 0.3], sigma = 2.0)])],
    transport_costs = Dict((p, o, r) => 5.0 for p in ("w", "g", "h")
                           for (o, r) in (("A", "B"), ("B", "A"))),
    armington = [Armington(product = p, sigma = 4.0,
                           shares = Dict((o, r) => o == r ? 3.0 : 1.0
                                         for o in ("A", "B"), r in ("A", "B")))
                 for p in ("w", "g", "h")],
    taxes = tariff20)
with_logger(NullLogger()) do
    solve_market(warmup)
end

GC.gc()
wall = @elapsed free = solve_market(market(Tax[]))
@printf("\nno tax: %s, %d variables, %d Ipopt iterations, Ipopt %.1f s, wall %.1f s\n",
        free.status, Markets.JuMP.num_variables(free.model), free.iterations, free.solve_time, wall)

println("\n20% ad valorem on every route:")
GC.gc()
wall = @elapsed taxed = solve_market(market(tariff20))
@printf("%s: %d solves, %d Ipopt iterations, Ipopt %.1f s, wall %.1f s\n",
        taxed.status, taxed.solves, taxed.iterations, taxed.solve_time, wall)

if "--cold" in ARGS
    println("\nthe same, every solve from scratch:")
    GC.gc()
    wall = @elapsed cold = solve_market(market(tariff20);
                                        solver_options = (warm_start_init_point = "no", mu_init = 0.1))
    @printf("%s: %d solves, %d Ipopt iterations, Ipopt %.1f s, wall %.1f s; objective %.8g vs %.8g\n",
            cold.status, cold.solves, cold.iterations, cold.solve_time, wall,
            cold.objective, taxed.objective)
end

# what the tariff did: revenue, against 20% of the origin value of what was
# shipped, and the trade it cut
pp = Dict((r, p) => x for (r, p, x) in zip(taxed.prices.region, taxed.prices.product,
                                           taxed.prices.producer_price))
implied = sum(0.2 * pp[(row.from, row.product)] * row.quantity for row in eachrow(taxed.trade)
              if row.from != row.to)
@printf("\nrevenue %.6g, 20%% of the origin value of the trade tables %.6g\n",
        only(taxed.taxes.revenue), implied)
traded(res) = sum(row.quantity for row in eachrow(res.trade) if row.from != row.to; init = 0.0)
@printf("cross-border trade %.1f without the tariff, %.1f with it (%.1f%%)\n",
        traded(free), traded(taxed), 100 * (traded(taxed) / traded(free) - 1))
@printf("welfare %.6g without the tariff, %.6g with it: a loss of %.4g\n",
        free.welfare, taxed.welfare, free.welfare - taxed.welfare)
