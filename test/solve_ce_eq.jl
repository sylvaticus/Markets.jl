using Pkg
cd(@__DIR__)
Pkg.activate(".")
using DataFrames, OdsIO, AxisArrays

# *** Getting input data and wrangling them ***
inputData       = "data/data.ods"
sets            = ods_read(inputData;sheetName="sets",retType="DataFrame")
const products  = collect(skipmissing(sets.products))
const regions   = collect(skipmissing(sets.regions))
(np, nr)        = length(products), length(regions)

elasticities_table = ods_read(inputData;sheetName="elasticities",retType="DataFrame")
contantTerms_table = ods_read(inputData;sheetName="constantTerms",retType="DataFrame")
init_table         = ods_read(inputData;sheetName="init",retType="DataFrame")

ϵ_d = AxisArray(zeros(np,nr,np,nr), p=products, r=regions, in_p = products, in_r=regions)
ϵ_s = AxisArray(zeros(np,nr,np,nr), p=products, r=regions, in_p = products, in_r=regions)

const_d = AxisArray(zeros(np,nr), p=products, r=regions)
const_s = AxisArray(zeros(np,nr), p=products, r=regions)
d0      = AxisArray(zeros(np,nr), p=products, r=regions)
s0      = AxisArray(zeros(np,nr), p=products, r=regions)
pr0     = AxisArray(zeros(np,nr), p=products, r=regions)

for r in eachrow(elasticities_table)
    ϵ_d[p = r.p, r = r.r, in_p = r.in_p, in_r = r.in_r] = r.d_value
    ϵ_s[p = r.p, r = r.r, in_p = r.in_p, in_r = r.in_r] = r.s_value
end
for r in eachrow(contantTerms_table)
    const_d[p = r.p, r = r.r] = r.const_d
    const_s[p = r.p, r = r.r] = r.const_s
end
for r in eachrow(init_table)
    d0[p = r.p, r = r.r]  = r.d0
    s0[p = r.p, r = r.r]  = r.s0
    pr0[p = r.p, r = r.r] = r.pr0
end

out = solveEquilibrium(const_d,ϵ_d,d0,const_s,ϵ_s,s0,pr0)


println("\n- Objective value (total costs): ", out.objective_value)
println("\n- Optimal prices:\n")
[println("($(products[p]),$(regions[r])): $(out.optPr[p,r])") for p in 1:np, r in 1:nr]
println("\n- Optimal demand quantities:\n")
[println("($(products[p]),$(regions[r])): $(out.optD[p,r])") for p in 1:np, r in 1:nr]
println("\n- Optimal supply quantities:\n")
[println("($(products[p]),$(regions[r])): $(out.optS[p,r])") for p in 1:np, r in 1:nr]

@test isapprox(out.optD, [141.555   152.322;115.118   100.0;306.072   267.472;487.378   645.038; 33.8543   26.5953], atol=0.001)