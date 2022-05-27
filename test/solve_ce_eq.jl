#using Pkg
cd(@__DIR__)
#Pkg.activate(".")
using DataFrames, OdsIO

# *** Getting input data and wrangling them ***
inputData       = "data/data.ods"
sets            = ods_read(inputData;sheetName="sets",retType="DataFrame")
const products  = collect(skipmissing(sets.products))
const regions   = collect(skipmissing(sets.regions))
(np, nr)        = length(products), length(regions)

elasticities_table = ods_read(inputData;sheetName="elasticities",retType="DataFrame")
constantTerms_table = ods_read(inputData;sheetName="constantTerms",retType="DataFrame")
init_table         = ods_read(inputData;sheetName="init",retType="DataFrame")

elasticities_dict  = Dict( [(r.p,r.r,r.in_p,r.in_r) => (r.d_value,r.s_value) for r in eachrow(elasticities_table)])
constantTerms_dict = Dict( [(r.p,r.r)               => (r.const_d,r.const_s) for r in eachrow(constantTerms_table)])
init_dict          = Dict( [(r.p,r.r)               => (r.d0,r.s0,r.pr0)     for r in eachrow(init_table)])

ϵ_d = zeros(np,nr,np,nr)
ϵ_s = zeros(np,nr,np,nr)

const_d = zeros(np,nr)
const_s = zeros(np,nr)
d0      = zeros(np,nr)
s0      = zeros(np,nr)
pr0     = zeros(np,nr)

for (ip,p) in enumerate(products), (ir,r) in enumerate(regions)
    for (ip2,p2) in enumerate(products), (ir2,r2) in enumerate(regions)
        ϵ_d[ip,ir,ip2,ir2] = get(elasticities_dict,(p,r,p2,r2),[0,0])[1]
        ϵ_s[ip,ir,ip2,ir2] = get(elasticities_dict,(p,r,p2,r2),[0,0])[2]
    end
    const_d[ip,ir] = get(constantTerms_dict,(p,r),[0,0])[1]
    const_s[ip,ir] = get(constantTerms_dict,(p,r),[0,0])[2]
    d0[ip,ir]      = get(init_dict,(p,r),[1,1,1])[1]
    s0[ip,ir]      = get(init_dict,(p,r),[1,1,1])[2]
    pr0[ip,ir]     = get(init_dict,(p,r),[1,1,1])[3]
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