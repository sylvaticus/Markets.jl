# A line-by-line transcription in JuMP of the GAMS model library's egypt.gms
# (Kutcher, 1980), used to check the port in egypt.jl: same sets, same
# variables, same equations, same objective. It reads the data that egypt.jl
# has already loaded.
#
# `grid = true` is the original: the demand curves are replaced by ten points
# each, and consumption is a convex combination of them, which makes the model
# an LP. `grid = false` keeps the same curves, but exact: the welfare of a
# consumption level is the integral of the linear demand, as in the port.

using Markets: JuMP, Ipopt

function egypt_gams(; grid::Bool, water_available = 30000 * 1000)
    m = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_silent(m)
    crops = unique(first.(collect(keys(yields))) ∪ first.(collect(keys(feed))) ∪
                   first.(collect(keys(land))))
    fodder_crops = unique(first.(collect(keys(feed))))
    ## the adjusted yield matrix yld(cn, c, r)
    yld(cn, c, r) = (y = get(yields, (c, r), 0.0);
                     cn == "veg-oil" && haskey(oil, c) ? y * oil[c] :
                     get(reclassified, c, c) == cn ? y : 0.0)
    netcs(c, r) = get(yields, (c, r), 0.0) > 0 ?
        get(tractor, (c, zone[r]), 0.0) * tractor_hour + get(cropdat, (c, "misc"), 0.0) +
        get(cropdat, (c, "pest"), 0.0) + get(cropdat, (c, "n-fer"), 0.0) * nitrogen_price +
        get(cropdat, (c, "p-fer"), 0.0) * phosphate_price -
        get(indprod, (c, r), 0.0) * get(byproduct_price, c, 0.0) : 0.0
    JuMP.@variables(m, begin
        xcrop[farms, crops] >= 0
        imports[marketed] >= 0
        exports[marketed] >= 0
        sales[marketed] >= 0
        cons[marketed]                     # national consumption, Σ qs·natq on the grid
        anut[["protein", "starch"], farms] >= 0
        trans[c in transported, r in farms, rp in farms; haskey(transfer_cost, (r, rp))] >= 0
        fodder[farms, fodder_crops] >= 0
        tlab[farms, months] >= 0
        0 <= flab[r in farms, months] <= 0.4 * agrreg[(r, "farmers")]
    end)
    for (c, r) in keys(upbnds)
        JuMP.set_upper_bound(xcrop[r, c], upbnds[(c, r)])
    end
    for c in marketed
        cap = get(trade_limits, (c, "imports"), 0.0)
        isinf(cap) || JuMP.set_upper_bound(imports[c], cap)
        cap = get(trade_limits, (c, "exports"), 0.0)
        isinf(cap) || JuMP.set_upper_bound(exports[c], cap)
    end
    ## demand: p = α + β·q, through the base point with the base elasticity
    β = Dict(c => demdat[(c, "base-p")] / demdat[(c, "domes-cons")] / demdat[(c, "elas")]
             for c in marketed)
    α = Dict(c => demdat[(c, "base-p")] - β[c] * demdat[(c, "domes-cons")] for c in marketed)
    qmin = Dict(c => (1.666 * demdat[(c, "base-p")] - α[c]) / β[c] for c in marketed)
    qmax = Dict(c => (0.6 * demdat[(c, "base-p")] - α[c]) / β[c] for c in marketed)
    if grid
        g = 1:10
        qs = Dict((c, k) => qmin[c] + (qmax[c] - qmin[c]) / 9 * (k - 1) for c in marketed, k in g)
        JuMP.@variable(m, natq[marketed, g] >= 0)
        JuMP.@constraint(m, conv[c in marketed], sum(natq[c, k] for k in g) == 1)
        JuMP.@constraint(m, [c in marketed], cons[c] == sum(qs[(c, k)] * natq[c, k] for k in g))
        welfare = JuMP.@expression(m, sum(natq[c, k] * (α[c] * qs[(c, k)] + 0.5 * β[c] * qs[(c, k)]^2)
                                          for c in marketed, k in g))
    else
        for c in marketed
            JuMP.set_lower_bound(cons[c], qmin[c])
            JuMP.set_upper_bound(cons[c], qmax[c])
        end
        welfare = JuMP.@expression(m, sum(α[c] * cons[c] + 0.5 * β[c] * cons[c]^2 for c in marketed))
    end
    JuMP.@constraints(m, begin
        comb[cn in marketed], sales[cn] <= sum(xcrop[r, c] * yld(cn, c, r) for r in farms, c in crops)
        demb[cn in marketed], exports[cn] + cons[cn] - sales[cn] * (1 - get(losses, cn, 0.0)) -
                              imports[cn] == 0
        landb[r in farms, tm in months],
            sum(xcrop[r, c] * get(land, (c, zone[r], tm), 0.0) for c in crops) <= agrreg[(r, "area")]
        labbal[r in farms, tm in months],
            sum(xcrop[r, c] * get(laborm, (c, zone[r], tm), 0.0) for c in crops) / working_days <=
            tlab[r, tm] + flab[r, tm]
        waterb, sum(xcrop[r, c] * get(water, (c, zone[r], tm), 0.0)
                    for r in farms, c in crops, tm in months) <= water_available
        vegetb[s in seasons, r in farms],
            sum(xcrop[r, c] for c in vegetables[s]) <= vegetable_area[(r, s)] * 1.5
        fodb[r in farms, cf in fodder_crops],
            fodder[r, cf] == xcrop[r, cf] * get(feed, (cf, zone[r]), 0.0) +
                (cf in transported ? sum(trans[cf, rp, r] - trans[cf, r, rp]
                                         for rp in farms if haskey(transfer_cost, (r, rp)); init = 0.0) : 0.0)
        nutb[nt in ["protein", "starch"], r in farms],
            sum(fodder[r, cf] * get(cropdat, (cf, nt), 0.0) * 0.01 for cf in fodder_crops) +
            anut[nt, r] >= agrreg[(r, nt)]
        gfodb[r in farms], sum(fodder[r, c] for c in green if c in fodder_crops) * green_dry_ratio <=
                           sum(fodder[r, c] for c in dry if c in fodder_crops)
    end)
    JuMP.@objective(m, Max, meat_and_milk + welfare +
        sum(exports[c] * demdat[(c, "export-p")] - imports[c] * demdat[(c, "import-p")]
            for c in marketed) -
        sum(trans[c, r, rp] * transfer_cost[(r, rp)] for (c, r, rp) in eachindex(trans)) -
        sum(xcrop[r, c] * netcs(c, r) for r in farms, c in crops) -
        sum(anut[nt, r] * concentrate_price for nt in ["protein", "starch"], r in farms) -
        sum((flab[r, tm] + 2 * tlab[r, tm]) * wages[(r, tm)] * working_days
            for r in farms, tm in months))
    JuMP.optimize!(m)
    return (status = JuMP.termination_status(m), objective = JuMP.objective_value(m),
            price = Dict(c => JuMP.shadow_price(demb[c]) for c in marketed),
            consumption = Dict(c => JuMP.value(cons[c]) for c in marketed),
            area = Dict((r, c) => JuMP.value(xcrop[r, c]) for r in farms, c in crops),
            water_rent = JuMP.shadow_price(waterb))
end
