using Pkg
cd(@__DIR__)
Pkg.activate(".")
using DataFrames, OdsIO

# *** Getting input data and wrangling them ***
inputData       = "data/data.ods"
sets            = ods_read(inputData;sheetName="sets",retType="DataFrame")
const products  = collect(skipmissing(sets.products))
const regions   = collect(skipmissing(sets.regions))
(np, nr)        = length(products), length(regions)

priceElasticities_table  = ods_read(inputData;sheetName="priceElasticities",retType="DataFrame")
tradeElasticities_table  = ods_read(inputData;sheetName="tradeElasticities",retType="DataFrame")
shareParameters_table          = ods_read(inputData;sheetName="shareParameters",retType="DataFrame")
constantTerms_table = ods_read(inputData;sheetName="constantTerms",retType="DataFrame")
init_table          = ods_read(inputData;sheetName="init",retType="DataFrame")
init2_table          = ods_read(inputData;sheetName="init2",retType="DataFrame")



priceElasticities_dict  = Dict( [(r.p,r.r,r.pIO) => (r.esI,r.esO,r.edI,r.edO,) for r in eachrow(priceElasticities_table)])
tradeElasticities_dict  = Dict( [(r.p) => (r.ets,r.etd) for r in eachrow(tradeElasticities_table)])
shareParameters_dict = Dict( [(r.p,r.r,r.r2)               => (r.a,r.b) for r in eachrow(shareParameters_table)])
constantTerms_dict = Dict( [(r.p,r.r)               => (r.const_s,r.const_d) for r in eachrow(constantTerms_table)])
init_dict          = Dict( [(r.p,r.r)               => (r.sc0,r.dc0,r.pcs0,r.pcd0)     for r in eachrow(init_table)])
init2_dict         = Dict( [(r.p,r.r,r.r2)               => (r.s0,r.d0,r.pa0)     for r in eachrow(init2_table)])

ϵsI     = zeros(np,nr,np)
ϵdI     = zeros(np,nr,np)
ϵsO     = zeros(np,nr,np)
ϵdO     = zeros(np,nr,np)
a       = zeros(np,nr,nr)
b       = zeros(np,nr,nr)
ϵts     = zeros(np)
ϵtd     = zeros(np)
const_s = zeros(np,nr)
const_d = zeros(np,nr)
sc0     = zeros(np,nr)
dc0     = zeros(np,nr)
s0      = zeros(np,nr,nr)
d0      = zeros(np,nr,nr)
pcs0    = zeros(np,nr)
pcd0    = zeros(np,nr)
pa0     = zeros(np,nr,nr)

for (ip,p) in enumerate(products)
    ϵts[ip] = get(tradeElasticities_dict,(p),[0,0])[1]
    ϵtd[ip] = get(tradeElasticities_dict,(p),[0,0])[1]
    for (ir,r) in enumerate(regions)
        for (ip2,p2) in enumerate(products)
            ϵsI[ip,ir,ip2] = get(priceElasticities_dict,(p,r,p2),[0,0,0,0])[1]
            ϵsO[ip,ir,ip2] = get(priceElasticities_dict,(p,r,p2),[0,0,0,0])[2]
            ϵdI[ip,ir,ip2] = get(priceElasticities_dict,(p,r,p2),[0,0,0,0])[3]
            ϵdO[ip,ir,ip2] = get(priceElasticities_dict,(p,r,p2),[0,0,0,0])[4]
        end
        for (ir2,r2) in enumerate(regions)
            s0[ip,ir,ir2]  = get(init2_dict,(p,r,r2),[0,0,0])[1]
            d0[ip,ir,ir2]  = get(init2_dict,(p,r,r2),[0,0,0])[2]
            pa0[ip,ir,ir2] = get(init_dict,(p,r),[1,1,1])[3]
            a[ip,ir,ir2]   = get(shareParameters_dict,(p,r,r2),[0,0])[1]
            b[ip,ir,ir2]   = get(shareParameters_dict,(p,r,r2),[0,0])[2]
        end
        const_s[ip,ir] = get(constantTerms_dict,(p,r),[0,0])[1]
        const_d[ip,ir] = get(constantTerms_dict,(p,r),[0,0])[2]
        sc0[ip,ir]      = get(init_dict,(p,r),[1,1,1,1])[1]
        dc0[ip,ir]      = get(init_dict,(p,r),[1,1,1,1])[2]
        pcs0[ip,ir]     = get(init_dict,(p,r),[1,1,1,1])[3]
        pcd0[ip,ir]     = get(init_dict,(p,r),[1,1,1,1,])[4]

    end
end

out = solveEquilibrium(const_s,ϵsI,ϵsO,a,ϵts,sc0,s0,const_d,ϵdI,ϵdO,b,ϵtd,dc0,d0,pcs0,pcd0,pa0)

out2 = runTestMktEq2(const_s,ϵ_s,a,ϵts,sc0,s0,const_d,ϵ_d,b,ϵtd,dc0,d0,pcs0,pcd0,pa0)



println("\n- Objective value (total costs): ", out.objective_value)
#=
println("\n- Optimal prices:\n")
[println("($(products[p]),$(regions[r])): $(out.optPr[p,r])") for p in 1:np, r in 1:nr]
println("\n- Optimal demand quantities:\n")
[println("($(products[p]),$(regions[r])): $(out.optD[p,r])") for p in 1:np, r in 1:nr]
println("\n- Optimal supply quantities:\n")
[println("($(products[p]),$(regions[r])): $(out.optS[p,r])") for p in 1:np, r in 1:nr]
=#
@test out.objective_value > -Inf # :-)

#@test isapprox(out.optD, [141.555   152.322;115.118   100.0;306.072   267.472;487.378   645.038; 33.8543   26.5953], atol=0.001)




# Market without trade

out = solveEquilibriumWithoutTrade()


println("\n- Objective value (total costs): ", out.objective_value)
println("\n- Optimal prices:\n")
[println("($(products[p]),$(regions[r])): $(out.optPr[p,r])") for p in 1:np, r in 1:nr]
println("\n- Optimal demand quantities:\n")
[println("($(products[p]),$(regions[r])): $(out.optD[p,r])") for p in 1:np, r in 1:nr]
println("\n- Optimal supply quantities:\n")
[println("($(products[p]),$(regions[r])): $(out.optS[p,r])") for p in 1:np, r in 1:nr]

@test isapprox(out.optD, [141.555   152.322;115.118   100.0;306.072   267.472;487.378   645.038; 33.8543   26.5953], atol=0.001)