# # Two regions and a dairy
#
# The smallest model that is still a market: two regions, two products, one
# transformation and one trade route. Everything else in this package is more
# of these same four things.
#
# The **Pasture** has the cows; the **City** has the people. Milk is either
# drunk as it is or turned into cheese, in either place, and either milk or
# cheese can be carried between the two at a cost. Four questions follow, and
# none of them has an obvious answer: how much of the milk is drunk and how
# much goes to the dairy, where the dairy stands, what each region pays, and
# what travels. Nobody decides. The answer is whatever clears both markets at
# once.

using Markets, DataFrames

# ## The economy
#
# Two regions and two products. Quantities are tonnes, prices euros per tonne.

regions  = [:Pasture, :City]
products = [:milk, :cheese]
nothing # hide

# Milk comes from the herds, on a curve through one observed point: last year
# the Pasture sold 2000 at €400 the tonne. The City keeps a few herds at its
# edges, where land is dear.

supply = [SupplySpec(product = :milk, region = :Pasture, p0 = 400, q0 = 2000, elasticity = 0.6),
          SupplySpec(product = :milk, region = :City,    p0 = 450, q0 =  300, elasticity = 0.6)]
nothing # hide

# Now the part that makes this more than an arithmetic exercise: **milk is both
# a final product and an input.** People drink it, and dairies turn it into
# cheese. Each use has its own demand, and neither is told how much milk it may
# have.

demand = [DemandSpec(product = :milk,   region = :Pasture, p0 =  400, q0 = 200, elasticity = 0.3),
          DemandSpec(product = :milk,   region = :City,    p0 =  450, q0 = 600, elasticity = 0.3),
          DemandSpec(product = :cheese, region = :Pasture, p0 = 5300, q0 =  30, elasticity = 0.5),
          DemandSpec(product = :cheese, region = :City,    p0 = 5350, q0 = 120, elasticity = 0.5)]
nothing # hide

# One transformation, and it may run in either region. Cheesemaking is a fixed
# recipe — ten tonnes of milk to the tonne of cheese, the rest leaving as whey
# and water — and €1300 covers everything the model does not name: the rennet,
# the cellar, the cheesemaker.

processes = [Process(name    = :dairy,
                     inputs  = [leontief(product = :milk, coeff = 10.0)],
                     outputs = [:cheese => 1.0],
                     vacost  = 1300)]
nothing # hide

# Both products can be carried either way, at €45 the tonne. Writing a route
# down is what makes a product tradable — anything with no route stays where it
# was made — and routes are directional, hence four entries for two products.
# Note that the rate is the *same* for both: nothing in the data says that one
# is cheaper to move than the other.

transport_costs = Dict((:milk,   :Pasture, :City) => 45.0, (:milk,   :City, :Pasture) => 45.0,
                       (:cheese, :Pasture, :City) => 45.0, (:cheese, :City, :Pasture) => 45.0)
nothing # hide

# ## Solving

economy = MarketData(; regions, products, demand, supply, processes, transport_costs)
result  = solve_market(economy)
nothing # hide

# ## What comes out
#
# Prices first. They are not assumptions: they are the shadow values of the
# material balances — what one more tonne of each product would be worth in
# each place.

unstack(result.prices, :region, :product, :price)

# They come back within a euro or two of the points the curves were anchored
# on, because the data were mutually consistent to begin with. A consistent
# benchmark reproduces itself, which is all that [Calibration](@ref) means.
#
# Then what is made, what is used, and what moves:

unstack(result.production, :region, :product, :quantity)

#-

unstack(result.consumption, :region, :product, :quantity)

#-

result.trade

# Three things in those tables were never put in.
#
# **The milk splits itself.** Of the 2301 tonnes the herds give, 801 are drunk
# and 1500 go to the dairy. Nothing in the data said how to divide them: the
# two uses bid against each other for the same milk, and the price that settles
# the contest is the one that makes every market balance at once.
#
# **The City makes no cheese at all** — its entry in the production table is
# `missing`, meaning the dairy never runs there. Bringing in a tonne of cheese
# costs €45 of freight; bringing in the ten tonnes of milk behind it costs
# €450. So every cheese eaten in the City was made in the Pasture, and the
# model was not told that, it worked it out.
#
# **Milk still crosses, but only the milk that is drunk.** The City's own herds
# cover part of its thirst and the rest is carried in. Not one tonne of milk is
# carried in order to be made into cheese.

# ## Two things worth checking by hand
#
# A model is worth no more than the conditions it satisfies, and in one this
# small each can be checked on a single line.
#
# **The dairy earns no profit.** Where it runs, the cheese price must cover
# exactly ten tonnes of milk plus the €1300 — any more and another dairy would
# open, any less and this one would shut:

price(r, p) = only(result.prices[(result.prices.region .== r) .&
                                 (result.prices.product .== p), :price])
(cheese = price(:Pasture, :cheese), milk_plus_making = 10 * price(:Pasture, :milk) + 1300)

# **Carrying earns no profit either.** On a route that is used, the price gap
# between its two ends is exactly the freight — any wider and somebody would
# ship more until it closed. Both routes carry something here, so both gaps
# are €45:

(milk   = price(:City, :milk)   - price(:Pasture, :milk),
 cheese = price(:City, :cheese) - price(:Pasture, :cheese))

# ## Asking it something
#
# The tables above are not the point — you could have measured those. The point
# is the question you cannot measure: *what if something changed?*
#
# Suppose the City takes to drinking milk, and wants a fifth more of it at
# every price. That is a **shifter** on the demand curve; the observed point
# the curve was fitted to stays exactly where it was.

thirsty = [d.product == :milk && d.region == :City ?
               DemandSpec(product = :milk, region = :City, p0 = d.p0, q0 = d.q0,
                          elasticity = d.elasticity, shift = 1.2) : d
           for d in demand]
thirst = solve_market(MarketData(; regions, products, demand = thirsty, supply,
                                   processes, transport_costs))

after(r, p) = only(thirst.prices[(thirst.prices.region .== r) .&
                                 (thirst.prices.product .== p), :price])
DataFrame((region = r, product = p, before = round(price(r, p); digits = 1),
           after  = round(after(r, p); digits = 1),
           change_pct = round(100 * (after(r, p) / price(r, p) - 1); digits = 1))
          for r in regions for p in products)

# Nobody's taste for *cheese* changed, and yet cheese is dearer — in the
# Pasture as much as in the City. The two uses of milk compete for the same
# herds, so a City that drinks more leaves less to be made into cheese, and
# everyone who eats cheese pays for it. A model with one curve per product and
# no transformation tying them together could not have produced that at all: it
# is the dual use of milk, and nothing else, that carries the effect across.
#
# Cheese also rises by much less, in percentage terms, than milk does — 4.2%
# against 5.6%. A good part of what a tonne of cheese costs is the making, and
# a thirstier City does not make a cheesemaker dearer.
#
# Behind the prices, the quantities tell the same story twice: the City drinks
# 602 tonnes of milk before and 711 after, while the dairy's output falls from
# 150 tonnes of cheese to 147. The extra glasses came out of the cheese.
#
# ## Where to go next
#
# * [Using the module](@ref) for the rest of the data types: joint products,
#   substitutable inputs, imperfect substitution between origins, bounds and
#   taxes.
# * [The forest sector: France in the world](@ref) for the same four ideas at
#   full size — eight regions, eleven products, and a storm.

#src ---------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__                                          #src
    println("\nPRICES");        show(unstack(result.prices, :region, :product, :price))         #src
    println("\n\nPRODUCTION");  show(unstack(result.production, :region, :product, :quantity))  #src
    println("\n\nCONSUMPTION"); show(unstack(result.consumption, :region, :product, :quantity)) #src
    println("\n\nTRADE");       show(result.trade)                             #src
    println()                                                                  #src
end                                                                            #src
