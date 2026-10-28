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

supply = [SupplyCurve(product = "milk", region = "Pasture", p0 = 400, q0 = 2000, elasticity = 0.6),
          SupplyCurve(product = "milk", region = "City",    p0 = 450, q0 =  300, elasticity = 0.6)]
nothing # hide

# Now the part that makes this more than an arithmetic exercise: **milk is both
# a final product and an input.** People drink it, and dairies turn it into
# cheese. Each use has its own demand, and neither is told how much milk it may
# have.

demand = [DemandCurve(product = "milk",   region = "Pasture", p0 =  400, q0 = 200, elasticity = 0.3),
          DemandCurve(product = "milk",   region = "City",    p0 =  450, q0 = 600, elasticity = 0.3),
          DemandCurve(product = "cheese", region = "Pasture", p0 = 5300, q0 =  30, elasticity = 0.5),
          DemandCurve(product = "cheese", region = "City",    p0 = 5350, q0 = 120, elasticity = 0.5)]
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
# ### Constraints
#
# The constraints are of two kinds. The first is a **material balance**,
# `balance[region,product]`, for every product in every region: what is
# supplied, made and imported equals what is consumed, used up and exported.
#
# ```math
# S_{r,p} + \sum_k y_{k,p}\, z_{r,k} + \sum_o T_{p,o,r}
#   \;=\; D_{r,p} + \sum_k u_{k,p}\, z_{r,k} + \sum_d T_{p,r,d}
# ```
#
# ``k`` runs over the processes, and this economy has only one, `dairy`. For
# each process ``k`` and product ``p``:
#
# * ``y_{k,p}`` is how much of ``p`` one unit of activity of ``k``
#   **produces**: its yield, from `outputs`;
# * ``u_{k,p}`` is how much of ``p`` one unit of activity of ``k`` **uses
#   up**: its input requirement, from `inputs`.
#
# For the dairy they are
#
# | product | ``y_{\text{dairy},p}`` (out) | ``u_{\text{dairy},p}`` (in) |
# |:--------|----:|----:|
# | milk    | 0   | 10  |
# | cheese  | 1   | 0   |
#
# so in region ``r``, ``\sum_k y_{k,p}\, z_{r,k}`` is everything the processes
# there make of ``p``, and ``\sum_k u_{k,p}\, z_{r,k}`` everything they use up
# of it. With one process each sum has a single term that is not zero, and the
# two balances of a region read
#
# ```
# milk:    S + 0·z + imports = D + 10·z + exports
# cheese:  0 + 1·z + imports = D +  0·z + exports
# ```
#
# which is where the coefficients of `z` in the printed model come from: −10 in
# the milk balance of its region, +1 in the cheese balance.
#
# What is transformed is in neither `S` nor `D`. `S` is **primary supply**
# only, what enters from outside the modelled processes: the milk of the
# herds. `D` is **final demand** only, what leaves the economy to be consumed:
# the milk people drink and the cheese they eat. The milk the dairy takes is
# the ``u\, z`` term, and the cheese it makes the ``y\, z`` one.
#
# That the balance is this simple comes from the dairy's recipe being fixed —
# a [`leontief`](@ref) input, ten tonnes of milk whatever the prices. Each
# ``u_{k,p}`` is then a constant, and what a process uses is that constant
# times its activity. A process whose inputs substitute for each other chooses
# its mix as prices move, and the quantity it uses of each becomes a variable
# of its own; the [forest example](@ref "The forest sector: France in the world")
# has those.
#
# The dual value of each balance is the price of that product in that region.
#
# The second kind of constraint is the **bounds**. `S`, `z` and `T` are only
# required to be non-negative. Each `D` is instead kept a little above zero, by
# a floor of a millionth of its own reference quantity — its `q0`, times its
# shift, which is 1 here:
#
# | variable            | `q0` | floor  |
# |:--------------------|-----:|-------:|
# | `D[Pasture,milk]`   | 200  | 0.0002 |
# | `D[City,milk]`      | 600  | 0.0006 |
# | `D[Pasture,cheese]` |  30  | 0.00003 |
# | `D[City,cheese]`    | 120  | 0.00012 |
#
# The printed model shows them as `0.00019999999999999998` and the like, which
# is only floating point: a millionth of 200 has no exact binary form.
#
# The floor is there because, with ``\eta < 1``, the demand term of the
# objective below raises `D` to a negative power — ``D^{-2.33}`` for milk,
# ``D^{-1}`` for cheese — which is infinite at zero, and so is its slope. The
# solver evaluates both at every step it takes, so `D` must never get there,
# not even on the way to the answer. At a solution the floor never binds: as
# consumption shrinks towards zero its price grows without limit, and nobody
# buys so little. `S` needs no floor, because it is raised to a power above 1
# (``S^{2.67}``), smooth at zero; nor do `z` and `T`, which enter only
# linearly.
#
# The floor is per variable because what counts as close to zero depends on
# the quantity: a millionth of its own reference is negligible for each,
# whether a product trades in tens of tonnes or in millions. A single absolute
# floor would have to suit the smallest product, and would then be far smaller
# than needed for the largest. The `qfloor` keyword of `jump_model` and
# `solve_market` replaces them with one value for every variable.
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
#   ``P^D(D) = a\, D^{-1/\eta}``, and the term is its antiderivative. Its
#   slope at any ``D`` is the price consumers would pay for one more unit.
# * **Cost of primary supply.** Each supply curve, with elasticity
#   ``\varepsilon``, has the inverse ``P^S(S) = b\, S^{1/\varepsilon}``, and
#   the term is the area under it: what the herds cost to keep, unit by unit.
# * **Cost of making**, ``c_k`` per unit of activity: the `vacost` of €1300
#   for each tonne the dairy makes, in either region.
# * **Cost of carrying**, ``\tau`` per unit shipped: the €45 of each route.
#
# The scales ``a`` and ``b`` are what put each curve through its reference
# point. The curve must give the price ``p_0`` at the quantity ``m q_0``, the
# reference quantity as shifted, so
#
# ```math
# a \left(m q_0\right)^{-1/\eta} = p_0
#   \;\;\Rightarrow\;\; a = p_0 \left(m q_0\right)^{1/\eta},
# \qquad
# b \left(m q_0\right)^{1/\varepsilon} = p_0
#   \;\;\Rightarrow\;\; b = p_0 \left(m q_0\right)^{-1/\varepsilon}
# ```
#
# and the coefficient the printed objective shows in front of each power is
# that scale divided by the exponent. For milk in the Pasture, ``p_0 = 400``,
# ``q_0 = 200``, ``m = 1`` and ``\eta = 0.3`` give
# ``a = 400 \times 200^{1/0.3} = 1.871 \times 10^{10}`` and the exponent
# ``1 - 1/0.3 = -2.33``, so the coefficient is
# ``1.871 \times 10^{10} / (-2.33) = -8.02 \times 10^9``. Every curve of the
# economy, worked out the same way:
#
# | term                | ``p_0`` | ``q_0`` | elasticity | scale                       | exponent | coefficient              |
# |:--------------------|------:|-----:|-----:|----------------------------:|-------:|-------------------------:|
# | `D[Pasture,milk]`   |  400  |  200 | 0.3  | ``a = 1.871 \times 10^{10}`` | −2.33  | ``-8.02 \times 10^{9}``  |
# | `D[City,milk]`      |  450  |  600 | 0.3  | ``a = 8.198 \times 10^{11}`` | −2.33  | ``-3.51 \times 10^{11}`` |
# | `D[Pasture,cheese]` | 5300  |   30 | 0.5  | ``a = 4.77 \times 10^{6}``   | −1     | ``-4.77 \times 10^{6}``  |
# | `D[City,cheese]`    | 5350  |  120 | 0.5  | ``a = 7.704 \times 10^{7}``  | −1     | ``-7.704 \times 10^{7}`` |
# | `S[Pasture,milk]`   |  400  | 2000 | 0.6  | ``b = 0.00126``              | 2.67   | ``0.000472``             |
# | `S[City,milk]`      |  450  |  300 | 0.6  | ``b = 0.0335``               | 2.67   | ``0.0126``               |
#
# The demand exponent and coefficient are both negative whenever ``\eta < 1``,
# and the term still rises with ``D``; for cheese, ``\eta = 0.5`` makes the
# exponent exactly −1. The sizes of the scales say nothing in themselves: ``a``
# is a price times a quantity raised to ``1/\eta``, and is large because
# ``1/\eta`` is.
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
               DemandCurve(d; shift = 1.2) : d
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
