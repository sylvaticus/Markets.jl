export solveEquilibrium

using JuMP, Ipopt

"""
    solveEquilibrium

Solve an economic equilibrium by passing elasticities parameters of CD-form supply/demand curves.
"""
function solveEquilibrium(const_s,ϵ_s,a,ϵts,sc0,s0,const_d,ϵ_d,b,ϵtd,dc0,d0,pcs0,pcd0,pl0;functionalForm="constantElasticity",debug=false)

    (np, nr) = size(ϵ_d)

    # *** Model declaration ***
    #m = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>1000, "tol"=>0.001, "linear_solver"=>"mumps", "print_level" =>1))
    m = Model(optimizer_with_attributes(Ipopt.Optimizer))

    # *** Parameters declarations (useful to avoid model reconstruction if these changes)  ***
    #@NLparameter(m, _ϵ_d[p = 1:np, r = 1:nr, in_p = 1:np, in_r = 1:nr ] == ϵ_d[p,r,in_p,in_r])

    # *** Variables declaration ***
    @variable(m, sc[p = 1:np, r = 1:nr] >= 0,  start=sc0[p,r])   # Supply, composite quantity
    @variable(m, dc[p = 1:np, r = 1:nr] >= 0,  start=dc0[p,r])   # Demand, composite quantity
    @variable(m, pcs[p = 1:np, r = 1:nr] >= 0,  start=pcs0[p,r])   # Supply, composite price
    @variable(m, pcd[p = 1:np, r = 1:nr] >= 0,  start=pcd0[p,r])   # Supply, composite price
    @variable(m, s[p = 1:np, r = 1:nr, rto = 1:nr] >= 0,  start=s0[p,r,rto])      # Supply, from region r to region rto
    @variable(m, d[p = 1:np, rfrom = 1:nr, r = 1:nr] >= 0,  start=d0[p,rfrom,r])  # Demand, from region rfrom to region r
    @variable(m, pl[p = 1:np, r = 1:nr] >= 0,  start=pl0[p,r])    # Local market prices


    # *** Constraints declaration (and definition) ***
    @NLconstraint(m, compSupply[p in 1:np, r in 1:nr],  sc[p,r]  == const_s[p,r] * prod(pcs[in_p,r]^ ϵ_s[p,r,in_p] for in_p in 1:np if ϵ_s[p,r,in_p] != 0.0) )
    @NLconstraint(m, compDemand[p in 1:np, r in 1:nr],  dc[p,r]  == const_d[p,r] * prod(pcd[in_p,r]^ ϵ_d[p,r,in_p] for in_p in 1:np if ϵ_d[p,r,in_p] != 0.0) )
    @NLconstraint(m, compSupplyAggregation[p in 1:np, r in 1:nr],  sc[p,r]  ==  sum(a[p,r,rto]   * s[p,r,rto]^   ϵts[p] for rto   in 1:nr if a[p,r,rto]   !=0 )^(1/ϵts[p]) )
    @NLconstraint(m, compDemandAggregation[p in 1:np, r in 1:nr],  dc[p,r]  ==  sum(b[p,rfrom,r] * d[p,rfrom,r]^ ϵtd[p] for rfrom in 1:nr if b[p,rfrom,r] !=0 )^(1/ϵtd[p]) )

    @NLconstraint(m, supplyByDest[p in 1:np, r in 1:nr, rto in 1:nr],   s[p,r,rto]   == a[p,r,rto]   * (pcs[p,r]/pl[p,rto])   ^ (1/(1-ϵts[p])) )
    @NLconstraint(m, demandByOrig[p in 1:np, rfrom in 1:nr, r in 1:nr], d[p,rfrom,r] == b[p,rfrom,r] * (pcd[p,r]/pl[p,rfrom]) ^ (1/(1-ϵtd[p])) )

    @NLconstraint(m, supplyMonetaryBalance[p in 1:np, r in 1:nr],  pcs[p,r] * sc[p,r]  ==  sum(pl[p,rto] * s[p,r,rto]     for rto in 1:nr   ) )
    @NLconstraint(m, demandMonetaryBalance[p in 1:np, r in 1:nr],  pcd[p,r] * dc[p,r]  ==  sum(pl[p,rfrom] * d[p,rfrom,r] for rfrom in 1:nr ) )
    @NLconstraint(m, physicalBalance[p in 1:np, r in 1:nr],  sum(s[p,r,rto] for rto in 1:nr) ==  sum(d[p,rfrom,r] for rfrom in 1:nr ) )



    # *** Objective function ***
    @NLobjective(m, Min,
        sum((d[p,rfrom,rto] - s[p,rfrom,rto])^2 for p in 1:np, rfrom in 1:nr, rto in 1:nr)
    )

    # *** Print human-readable version of the model ***
    #print(m)

    # *** Solve the model  ***
    optimize!(m)
    status = termination_status(m)

    # *** Printing outputs (optimal quantities/prices) ***
    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT  || status == MOI.ITERATION_LIMIT) && has_values(m)
      return (solved=true, model=m, status=status, objective_value=objective_value(m), optD=value.(d), optS  = value.(s), optPl = value.(pl))
    else
      return (solved=false, model=m, status=status, objective_value=objective_value(m), optD=zeros(np,nr,nr), optS  = zeros(np,nr,nr), optPl = zeros(np,nr))
    end
end