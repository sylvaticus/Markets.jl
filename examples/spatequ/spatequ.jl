# # Takayama and Judge: two products, three regions
#
# The spatial equilibrium of [Takayama and Judge
# (1964)](https://doi.org/10.2307/1236473) is where models like this one
# started: regions with their own supply and demand, a cost of moving goods
# between them, and an equilibrium found by maximising the surplus net of
# transport. Its numerical example is
# [`spatequ`](https://www.gams.com/latest/gamslib_ml/libhtml/gamslib_spatequ.html)
# in the GAMS model library, contributed by Phan Sy Hieu, and this page ports
# it.
#
# It is small enough to check by hand, which is the point: a port is only
# worth something once it is shown to reproduce what it ports.

using Markets, DataFrames

# ## The data
#
# Three regions, two commodities, and straight-line supply and demand:
# quantity as a function of price, `Q = α + β·P`.

regions  = ["Reg1", "Reg2", "Reg3"]
products = ["Com1", "Com2"]
nothing # hide

# The GAMS tables give the intercepts, and a 2×2 matrix of price coefficients
# per region, with cross-price terms.

alpha_d = Dict(("Reg1", "Com1") => 200, ("Reg1", "Com2") => 300,
               ("Reg2", "Com1") => 100, ("Reg2", "Com2") => 200,
               ("Reg3", "Com1") => 160, ("Reg3", "Com2") => 250)
beta_d  = Dict(("Reg1", "Com1") => [-10, 1], ("Reg1", "Com2") => [1, -10],
               ("Reg2", "Com1") => [-5, 1],  ("Reg2", "Com2") => [1, -20],
               ("Reg3", "Com1") => [-8, 1],  ("Reg3", "Com2") => [1, -10])
alpha_s = Dict((r, c) => (c == "Com1" ? -50 : -60) for r in regions, c in products)
beta_s  = Dict(("Reg1", "Com1") => [10, 0.5], ("Reg1", "Com2") => [0.5, 15],
               ("Reg2", "Com1") => [20, 0.5], ("Reg2", "Com2") => [0.5, 25],
               ("Reg3", "Com1") => [10, 0.5], ("Reg3", "Com2") => [0.5, 15])
nothing # hide

# **What the GAMS model actually solves.** Its demand equation is
#
# ```
# DEM(r,c)..  AlphaD(r,c) + sum(cc, BetaD(r,c,cc)*P(r,c)) =e= Qd(r,c);
# ```
#
# and it multiplies every coefficient by `P(r,c)`, the product's *own* price,
# not by `P(r,cc)`; the supply equation does the same. As written, then, the
# cross-price terms are folded into the own-price slope, the two commodities
# are independent markets, and the slope of each is the sum of its row:
# −10 + 1 = −9 for Com1 in Reg1. That is what its LP and its complementarity
# version solve, and what is ported here. (A demand system with cross-price
# terms between two products would be a [`DemandNest`](@ref) in this
# package, whose substitution is a CES rather than a free linear matrix.)

slope(row) = sum(row)
demand = [DemandCurve(product = c, region = r, alpha = alpha_d[(r, c)],
                     beta = slope(beta_d[(r, c)]))
          for r in regions for c in products]
supply = [SupplyCurve(product = c, region = r, alpha = alpha_s[(r, c)],
                     beta = slope(beta_s[(r, c)]))
          for r in regions for c in products]
demand[1]

# Transport costs per unit, the same each way:

cost = Dict(("Reg1", "Reg2") => (2, 3), ("Reg1", "Reg3") => (2, 3),
            ("Reg2", "Reg3") => (1, 2))
transport_costs = Dict{Tuple{String,String,String},Float64}()
for ((a, b), (c1, c2)) in cost, (from, to) in ((a, b), (b, a))
    transport_costs[("Com1", from, to)] = c1
    transport_costs[("Com2", from, to)] = c2
end
nothing # hide

# ## The equilibrium

market = MarketData(regions = regions, products = products, demand = demand,
                    supply = supply, transport_costs = transport_costs)
result = solve_market(market)
result.prices

# Reg2 is where both commodities are cheapest: its supply responds most to the
# price, and its demand is the smallest. Both move out of it to the two other
# regions, whose prices sit exactly one transport cost above.

result.trade

# ## Checking the port
#
# The equilibrium of the GAMS model is characterised by its complementarity
# version: every quantity on its curve, every region in balance, no price gap
# wider than the cost of closing it, and an exact gap wherever goods move.
# Each condition can be checked on the result directly.

price_of(r, c) = only(result.prices[(result.prices.region .== r) .& (result.prices.product .== c), :price])
quantity(df, r, c) = (rows = df[(df.region .== r) .& (df.product .== c), :quantity];
                      isempty(rows) ? 0.0 : only(rows))
shipped(c, from, to) = (rows = result.trade[(result.trade.product .== c) .& (result.trade.from .== from) .&
                                            (result.trade.to .== to), :quantity];
                        isempty(rows) ? 0.0 : only(rows))
checks = DataFrame(region = String[], product = String[], price = Float64[],
                   demand_gap = Float64[], supply_gap = Float64[], balance_gap = Float64[])
for r in regions, c in products
    p = price_of(r, c)
    qd, qs = quantity(result.consumption, r, c), quantity(result.production, r, c)
    inflow  = sum(shipped(c, o, r) for o in regions if o != r)
    outflow = sum(shipped(c, r, o) for o in regions if o != r)
    push!(checks, (r, c, p,
                   qd - (alpha_d[(r, c)] + slope(beta_d[(r, c)]) * p),
                   qs - (alpha_s[(r, c)] + slope(beta_s[(r, c)]) * p),
                   qs + inflow - outflow - qd))
end
checks

# The gaps are at the solver's tolerance. And between any two regions, the
# price difference never exceeds the transport cost, and equals it on every
# route used:

arbitrage = DataFrame(product = String[], from = String[], to = String[],
                      price_gap = Float64[], cost = Float64[], shipped = Float64[])
for ((c, from, to), τ) in sort(collect(transport_costs))
    push!(arbitrage, (c, from, to, price_of(to, c) - price_of(from, c), τ, shipped(c, from, to)))
end
arbitrage

# Those are the conditions the GAMS complementarity model `P2R3_MCP` imposes,
# so this is its solution. What took four model statements there is one
# [`MarketData`](@ref) here, and the same economy can now be given processes,
# limits, taxes or imperfect substitution between origins without being
# written again.
