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

regions  = ["Pasture", "City"]
products = ["milk", "cheese"]
nothing # hide

# Milk comes from the herds, on a curve through one observed point: last year
# the Pasture sold 2000 at €400 the tonne. The City keeps a few herds at its
# edges, where land is dear.

supply = [SupplySpec(product = "milk", region = "Pasture", p0 = 400, q0 = 2000, elasticity = 0.6),
          SupplySpec(product = "milk", region = "City",    p0 = 450, q0 =  300, elasticity = 0.6)]
nothing # hide

# Now the part that makes this more than an arithmetic exercise: **milk is both
# a final product and an input.** People drink it, and dairies turn it into
# cheese. Each use has its own demand, and neither is told how much milk it may
# have.

demand = [DemandSpec(product = "milk",   region = "Pasture", p0 =  400, q0 = 200, elasticity = 0.3),
          DemandSpec(product = "milk",   region = "City",    p0 =  450, q0 = 600, elasticity = 0.3),
          DemandSpec(product = "cheese", region = "Pasture", p0 = 5300, q0 =  30, elasticity = 0.5),
          DemandSpec(product = "cheese", region = "City",    p0 = 5350, q0 = 120, elasticity = 0.5)]
nothing # hide

# One transformation, and it may run in either region. Cheesemaking is a fixed
# recipe — ten tonnes of milk to the tonne of cheese, the rest leaving as whey
# and water — and €1300 covers everything the model does not name: the rennet,
# the cellar, the cheesemaker.

processes = [Process(name    = "dairy",
                     inputs  = [leontief(product = "milk", coeff = 10.0)],
                     outputs = ["cheese" => 1.0],
                     vacost  = 1300)]
nothing # hide

# Both products can be carried either way, at €45 the tonne. Writing a route
# down is what makes a product tradable — anything with no route stays where it
# was made — and routes are directional, hence four entries for two products.
# Note that the rate is the *same* for both: nothing in the data says that one
# is cheaper to move than the other.

transport_costs = Dict(("milk",   "Pasture", "City") => 45.0, ("milk",   "City", "Pasture") => 45.0,
                       ("cheese", "Pasture", "City") => 45.0, ("cheese", "City", "Pasture") => 45.0)
nothing # hide

# ## Looking before solving
#
# The pieces go into one `MarketData`. Typing its name gives back only how much
# it holds; printing it gives a summary of the economy that was actually
# written down: what each product does, how the processes chain, where goods
# may go and at what cost. It is the cheapest check there is that the data say
# what was meant.

economy = MarketData(; regions, products, demand, supply, processes, transport_costs)
#-
print(economy)

# The mathematical program it becomes can be read too, before anything is
# solved: [`jump_model`](@ref) returns the JuMP model that `solve_market` would
# hand to the solver, with every variable and constraint named after the region
# and product it belongs to, and `latex_formulation` typesets it.

latex_formulation(jump_model(economy))

# !!! warning "Check the signs against `print`"
#     JuMP 1.32 drops the minus sign of a negative coefficient written in
#     scientific notation when it typesets a model
#     ([JuMP.jl#4268](https://github.com/jump-dev/JuMP.jl/issues/4268)). The
#     two milk terms at the start of the objective above are such
#     coefficients: both are negative, as `print(jump_model(economy))` shows.
#     Until the bug is fixed, read the signs from `print`.
#
# ### Reading the program
#
# There are four kinds of variable, and each is named after where it lives:
#
# * `D[region,product]` is **final demand**, the quantity consumed:
#   `D[City,milk]` is the milk the City drinks.
# * `S[region,product]` is **primary supply**: `S[Pasture,milk]` is the milk
#   the Pasture's herds give.
# * `z[region,process]` is the **activity level** of a process in a region:
#   `z[Pasture,dairy]` is how much the Pasture's dairy works, and
#   `z[City,dairy]` how much the City's does. The first index is the location
#   and the second names the process; neither is an input or an output. What
#   the dairy takes and gives is in its recipe: each unit of activity takes 10
#   tonnes of milk and yields 1 tonne of cheese, so here one unit of activity
#   is one tonne of cheese.
# * `T[product,from,to]` is the quantity **shipped** on a route:
#   `T[milk,Pasture,City]` is the milk carried from the Pasture to the City.
#   There is one for each route declared in `transport_costs`, so four here.
#
# The constraints are of two kinds. The first is a **material balance**,
# `balance[region,product]`, for every product in every region: what is
# supplied, made and imported equals what is consumed, used up and exported.
# With ``y_{k,p}`` the yield of product ``p`` per unit of activity of process
# ``k``, and ``u_{k,p}`` the quantity of it the process uses,
#
# ```math
# S_{r,p} + \sum_k y_{k,p}\, z_{r,k} + \sum_o T_{p,o,r}
#   \;=\; D_{r,p} + \sum_k u_{k,p}\, z_{r,k} + \sum_d T_{p,r,d}
# ```
#
# which is where the coefficients of `z` come from: −10 in the milk balance of
# its region, +1 in the cheese balance. The dual value of each balance is the
# price of that product in that region. The second kind are **bounds**: every
# quantity is non-negative, and each `D` is kept a millionth of its reference
# quantity away from zero, because the demand term below cannot be
# differentiated at zero. That floor never binds at a solution.
#
# ### The objective
#
# What is maximised is the **net social surplus**: what consumers would pay
# for what they get, less what it costs to produce, make and carry. Written
# with its parameters rather than its numbers, it is
#
# ```math
# \max_{D,\,S,\,z,\,T}\;
#   \sum_{r,p} \frac{a_{r,p}}{1 - 1/\eta_{r,p}}\, D_{r,p}^{\,1 - 1/\eta_{r,p}}
#   \;-\; \sum_{r,p} \frac{b_{r,p}}{1 + 1/\varepsilon_{r,p}}\, S_{r,p}^{\,1 + 1/\varepsilon_{r,p}}
#   \;-\; \sum_{r,k} c_k\, z_{r,k}
#   \;-\; \sum_{p,o,d} \tau_{p,o,d}\, T_{p,o,d}
# ```
#
# Each of the four sums is one block of the objective printed above.
#
# * **Gross consumer benefit.** Each demand curve, with reference point
#   ``(q_0, p_0)``, elasticity ``\eta`` and shifter ``m``, has the inverse
#   ``P^D(D) = a\, D^{-1/\eta}`` with ``a = p_0 (m q_0)^{1/\eta}``, and the
#   term is its antiderivative. Its slope at any ``D`` is the price consumers
#   would pay for one more unit. For milk in the Pasture, ``p_0 = 400``,
#   ``q_0 = 200`` and ``\eta = 0.3`` give the exponent
#   ``1 - 1/0.3 = -2.33`` and the coefficient ``-8.02 \times 10^9``. Both are
#   negative whenever ``\eta < 1``, and the term still rises with ``D``. For
#   cheese, ``\eta = 0.5`` makes the exponent exactly −1.
# * **Cost of primary supply.** Each supply curve has the inverse
#   ``P^S(S) = b\, S^{1/\varepsilon}`` with ``b = p_0 (m q_0)^{-1/\varepsilon}``,
#   and the term is the area under it: what the herds cost to keep, unit by
#   unit. With ``\varepsilon = 0.6`` the exponent is ``1 + 1/0.6 = 2.67``.
# * **Cost of making**, ``c_k`` per unit of activity: the `vacost` of €1300
#   for each tonne the dairy makes, in either region.
# * **Cost of carrying**, ``\tau`` per unit shipped: the €45 of each route.
#
# Why maximising this gives the market is in the conditions for its optimum.
# With ``\pi_{r,p}`` the dual of each balance, they say that consumers buy until
# the price they would pay for one more unit, ``P^D(D)``, falls to
# ``\pi_{r,p}``. Herds produce until their marginal cost ``P^S(S)`` rises to
# it. A dairy runs only where
# ``\pi_{r,\text{cheese}} = 10\, \pi_{r,\text{milk}} + 1300``. A route is used
# only where the price gap across it is exactly ``\tau``. Those
# are the conditions of a competitive market, and the last two are checked by
# hand further down. The general program, with taxes, bounds and imperfect
# substitutes, is on the [Code implementation](@ref) page.

# ## Solving

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
(cheese = price("Pasture", "cheese"), milk_plus_making = 10 * price("Pasture", "milk") + 1300)

# **Carrying earns no profit either.** On a route that is used, the price gap
# between its two ends is exactly the freight — any wider and somebody would
# ship more until it closed. Both routes carry something here, so both gaps
# are €45:

(milk   = price("City", "milk")   - price("Pasture", "milk"),
 cheese = price("City", "cheese") - price("Pasture", "cheese"))

# ## Asking it something
#
# The tables above are not the point — you could have measured those. The point
# is the question you cannot measure: *what if something changed?*
#
# Suppose the City takes to drinking milk, and wants a fifth more of it at
# every price. That is a **shifter** on the demand curve; the observed point
# the curve was fitted to stays exactly where it was.

thirsty = [d.product == "milk" && d.region == "City" ?
               DemandSpec(product = "milk", region = "City", p0 = d.p0, q0 = d.q0,
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
