export solveEquilibrium, solveEquilibriumWithoutTrade, runTestMktEq, runTestMktEq2

using JuMP
using Ipopt
using AmplNLWriter, Ipopt_jll

"""
    solveEquilibrium

Solve an economic equilibrium by passing elasticities parameters of CD-form supply/demand curves.
"""
function solveEquilibrium(const_s,ϵsI,ϵsO,a,ϵts,sc0,s0,const_d,ϵdI,ϵdO,b,ϵtd,dc0,d0,pcs0,pcd0,pa0;functionalForm="constantElasticity",debug=false)

    (np, nr, _) = size(ϵdI)

    # *** Model declaration ***
    #m = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>1000, "tol"=>0.001, "linear_solver"=>"mumps", "print_level" =>1))
    m = Model(optimizer_with_attributes(Ipopt.Optimizer))
    #m = Model(() -> AmplNLWriter.Optimizer("ipopt",["print_level=0"]))
    #m = Model(() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
    set_optimizer_attribute(m, "tol", 0.01)

    #parCheck(nzpar,expr) = (nzpar == 0 ? 0.0 : return expr)
    #register(m, :parCheck, 2, parCheck, autodiff=true)

    # *** Parameters declarations (useful to avoid model reconstruction if these changes)  ***
    #@NLparameter(m, _ϵ_d[p = 1:np, r = 1:nr, in_p = 1:np, in_r = 1:nr ] == ϵ_d[p,r,in_p,in_r])

    # *** Variables declaration ***
    @variable(m, sc[p = 1:np, r = 1:nr] >= 0,  start=sc0[p,r])   # Supply, composite quantity
    @variable(m, dc[p = 1:np, r = 1:nr] >= 0,  start=dc0[p,r])   # Demand, composite quantity
    @variable(m, pcs[p = 1:np, r = 1:nr] >= 0,  start=pcs0[p,r])   # Supply, composite price
    @variable(m, pcd[p = 1:np, r = 1:nr] >= 0,  start=pcd0[p,r])   # Supply, composite price
    @variable(m, s[p = 1:np, r = 1:nr, rto = 1:nr] >= 0,  start=s0[p,r,rto])       # Supply, from region r to region rto
    @variable(m, d[p = 1:np, rfrom = 1:nr, r = 1:nr] >= 0,  start=d0[p,rfrom,r])   # Demand, from region rfrom to region r
    @variable(m, pa[p = 1:np, rin = 1:nr, rout=1:nr] >= 0,  start=pa0[p,rin,rout]) # Acrtual price of product p coming from region rfrom and going to region rto


    # *** Constraints declaration (and definition) ***
    @NLconstraint(m, compSupply[p in 1:np, r in 1:nr],  sc[p,r]  == const_s[p,r] * prod(pcd[pI,r]^ ϵsI[p,r,pI] for pI in 1:np if ϵsI[p,r,pI] != 0.0) * prod(pcs[pO,r]^ ϵsO[p,r,pO] for pO in 1:np if ϵsO[p,r,pO] != 0.0) ) # Note we use two elasticities because pcd[p,r] and pcs[p,r] are not generally the same
    @NLconstraint(m, compDemand[p in 1:np, r in 1:nr],  dc[p,r]  == const_d[p,r] * prod(pcd[pI,r]^ ϵdI[p,r,pI] for pI in 1:np if ϵdI[p,r,pI] != 0.0) * prod(pcs[pO,r]^ ϵdO[p,r,pO] for pO in 1:np if ϵdO[p,r,pO] != 0.0))
    @NLconstraint(m, compSupplyAggregation[p in 1:np, r in 1:nr],  sc[p,r]  ==  sum(a[p,r,rto]   * s[p,r,rto]^   ϵts[p] for rto   in 1:nr if a[p,r,rto]   !=0 )^(1/ϵts[p]) )
    @NLconstraint(m, compDemandAggregation[p in 1:np, r in 1:nr],  dc[p,r]  ==  sum(b[p,rfrom,r] * d[p,rfrom,r]^ ϵtd[p] for rfrom in 1:nr if b[p,rfrom,r] !=0 )^(1/ϵtd[p]) )
    @NLconstraint(m, supplyByDest[p in 1:np, r in 1:nr, rto in 1:nr],   s[p,r,rto]   == ifelse(a[p,r,rto]   == 0, 0.0, a[p,r,rto]   * (pcs[p,r]/pa[p,r,rto])   ^ (1/(1-ϵts[p])) ) )
    @NLconstraint(m, demandByOrig[p in 1:np, rfrom in 1:nr, r in 1:nr], d[p,rfrom,r] == ifelse(b[p,rfrom,r] == 0, 0.0, b[p,rfrom,r] * (pcd[p,r]/pa[p,rfrom,r]) ^ (1/(1-ϵtd[p])) ))
    @NLconstraint(m, supplyMonetaryBalance[p in 1:np, r in 1:nr],  pcs[p,r] * sc[p,r]        ==  sum(pa[p,r,rto]   * s[p,r,rto]   for rto in 1:nr   ) )
    @NLconstraint(m, demandMonetaryBalance[p in 1:np, r in 1:nr],  pcd[p,r] * dc[p,r]        ==  sum(pa[p,rfrom,r] * d[p,rfrom,r] for rfrom in 1:nr ) )
    @NLconstraint(m, physicalBalance[p in 1:np, r in 1:nr],  sum(s[p,r,rto] for rto in 1:nr) ==  sum(d[p,rfrom,r] for rfrom in 1:nr ) )



    # *** Objective function ***
    @NLobjective(m, Min,
        sum((d[p,rfrom,rto] - s[p,rfrom,rto])^2 for p in 1:np, rfrom in 1:nr, rto in 1:nr)
    )

    # *** Print human-readable version of the model ***
    print(m)

    # *** Solve the model  ***
    optimize!(m)
    status = termination_status(m)
  
    println(status)
    # *** Printing outputs (optimal quantities/prices) ***
    #if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT  || status == MOI.ITERATION_LIMIT) && has_values(m)
    #  return (solved=true, model=m, status=status, objective_value=objective_value(m), optD=value.(d), optS  = value.(s), optPl = value.(pl))
    #else
    #  return (solved=false, model=m, status=status, objective_value=objective_value(m), optD=zeros(np,nr,nr), optS  = zeros(np,nr,nr), optPl = zeros(np,nr))
    #end
end

function runTestMktEq(par=4)

  # Qd = -da p + db
  da = par
  db = 20
  # Qs = sa p - sb
  sa = 8
  sb = 4
  m = Model(() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
  @variable(m, q >=0)
  @NLconstraint(m,testContraint,  q^2 == ifelse(da<2,14,20000))
  @NLconstraint(m,test2Contraint,  q <= 2)
  @objective(m,Max, -(1/da)*(1/2)*q^2 + (db/da)*q -((1/sa)*(1/2)*q^2+(sb/sa)*q))

  print(m)

  optimize!(m)
  status = termination_status(m)

  println(status)
  #println("Objective value: ", objective_value(m))
  #println("Q: ", value.(q))

end


function runTestMktEq2(const_s,ϵ_s,a,ϵts,sc0,s0,const_d,ϵ_d,b,ϵtd,dc0,d0,pcs0,pcd0,pa0;functionalForm="constantElasticity",debug=false)

  (np, nr, _) = size(ϵ_d)

  # *** Model declaration ***
  #m = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>1000, "tol"=>0.001, "linear_solver"=>"mumps", "print_level" =>1))
  #m = Model(optimizer_with_attributes(Ipopt.Optimizer))
  #m = Model(() -> AmplNLWriter.Optimizer("ipopt",["print_level=0"]))
  m = Model(() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe))

  #parCheck(nzpar,expr) = (nzpar == 0 ? 0.0 : return expr)
  #register(m, :parCheck, 2, parCheck, autodiff=true)

  # *** Parameters declarations (useful to avoid model reconstruction if these changes)  ***
  #@NLparameter(m, _ϵ_d[p = 1:np, r = 1:nr, in_p = 1:np, in_r = 1:nr ] == ϵ_d[p,r,in_p,in_r])

  # *** Variables declaration ***
  @variable(m, sc[p = 1:np, r = 1:nr] >= 0,  start=sc0[p,r])   # Supply, composite quantity
  @variable(m, dc[p = 1:np, r = 1:nr] >= 0,  start=dc0[p,r])   # Demand, composite quantity
  @variable(m, pcs[p = 1:np, r = 1:nr] >= 0,  start=pcs0[p,r])   # Supply, composite price
  @variable(m, pcd[p = 1:np, r = 1:nr] >= 0,  start=pcd0[p,r])   # Supply, composite price
  @variable(m, s[p = 1:np, r = 1:nr, rto = 1:nr] >= 0,  start=s0[p,r,rto])       # Supply, from region r to region rto
  @variable(m, d[p = 1:np, rfrom = 1:nr, r = 1:nr] >= 0,  start=d0[p,rfrom,r])   # Demand, from region rfrom to region r
  @variable(m, pa[p = 1:np, rin = 1:nr, rout=1:nr] >= 0,  start=pa0[p,rin,rout]) # Acrtual price of product p coming from region rfrom and going to region rto


  # *** Constraints declaration (and definition) ***
  @NLconstraint(m, compSupply[p in 1:np, r in 1:nr],  sc[p,r]  == const_s[p,r] * prod(pcs[in_p,r]^ ϵ_s[p,r,in_p] for in_p in 1:np if ϵ_s[p,r,in_p] != 0.0) )
  #@NLconstraint(m, compDemand[p in 1:np, r in 1:nr],  dc[p,r]  == const_d[p,r] * prod(pcd[in_p,r]^ ϵ_d[p,r,in_p] for in_p in 1:np if ϵ_d[p,r,in_p] != 0.0) )
  #@NLconstraint(m, compSupplyAggregation[p in 1:np, r in 1:nr],  sc[p,r]  ==  sum(a[p,r,rto]   * s[p,r,rto]^   ϵts[p] for rto   in 1:nr if a[p,r,rto]   !=0 )^(1/ϵts[p]) )
  #@NLconstraint(m, compDemandAggregation[p in 1:np, r in 1:nr],  dc[p,r]  ==  sum(b[p,rfrom,r] * d[p,rfrom,r]^ ϵtd[p] for rfrom in 1:nr if b[p,rfrom,r] !=0 )^(1/ϵtd[p]) )
  #@NLconstraint(m, supplyByDest[p in 1:np, r in 1:nr, rto in 1:nr],   s[p,r,rto]   == ifelse(a[p,r,rto]   == 0, 0.0, a[p,r,rto]   * (pcs[p,r]/pa[p,r,rto])   ^ (1/(1-ϵts[p])) ) )
  #@NLconstraint(m, demandByOrig[p in 1:np, rfrom in 1:nr, r in 1:nr], d[p,rfrom,r] == ifelse(b[p,rfrom,r] == 0, 0.0, b[p,rfrom,r] * (pcd[p,r]/pa[p,rfrom,r]) ^ (1/(1-ϵtd[p])) ))
  #@NLconstraint(m, supplyMonetaryBalance[p in 1:np, r in 1:nr],  pcs[p,r] * sc[p,r]        ==  sum(pa[p,r,rto]   * s[p,r,rto]   for rto in 1:nr   ) )
  #@NLconstraint(m, demandMonetaryBalance[p in 1:np, r in 1:nr],  pcd[p,r] * dc[p,r]        ==  sum(pa[p,rfrom,r] * d[p,rfrom,r] for rfrom in 1:nr ) )
  #@NLconstraint(m, physicalBalance[p in 1:np, r in 1:nr],  sum(s[p,r,rto] for rto in 1:nr) ==  sum(d[p,rfrom,r] for rfrom in 1:nr ) )



  # *** Objective function ***
  @NLobjective(m, Min,
      sum((d[p,rfrom,rto] - s[p,rfrom,rto])^2 for p in 1:np, rfrom in 1:nr, rto in 1:nr)
  )

  # *** Print human-readable version of the model ***
  print(m)

  # *** Solve the model  ***
  optimize!(m)
  status = termination_status(m)

  println(status)
  # *** Printing outputs (optimal quantities/prices) ***
  #if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT  || status == MOI.ITERATION_LIMIT) && has_values(m)
  #  return (solved=true, model=m, status=status, objective_value=objective_value(m), optD=value.(d), optS  = value.(s), optPl = value.(pl))
  #else
  #  return (solved=false, model=m, status=status, objective_value=objective_value(m), optD=zeros(np,nr,nr), optS  = zeros(np,nr,nr), optPl = zeros(np,nr))
  #end
end



function solveEquilibriumWithoutTrade(
  const_d = [50.0 100.0; 100.0 10.0; 200.0 200.0; 200.0 200.0; 40.0 10.0],
  ϵ_d = [-0.8 0.4; 0.3 0.15; 0.0 0.0; 0.0 0.0; 0.0 0.0;;; 0.3 0.15; -0.8 0.4; 0.0 0.0; 0.0 0.0; 0.0 0.0;;; 0.6 0.3; 0.0 0.0; -0.8 0.4; 0.4 0.2; 0.4 0.2;;; 0.0 0.0; 0.6 0.15; 0.4 0.2; -0.8 0.4; 0.4 0.2;;; 0.3 0.15; 0.3 0.15; 0.4 0.2; 0.4 0.2; -0.8 0.4;;;; 0.4 -0.8; 0.15 0.3; 0.0 0.0; 0.0 0.0; 0.0 0.0;;; 0.15 0.3; 0.4 -0.8; 0.0 0.0; 0.0 0.0; 0.0 0.0;;; 0.3 0.6; 0.0 0.0; 0.4 -0.8; 0.2 0.4; 0.2 0.4;;; 0.0 0.0; 0.15 0.6; 0.2 0.4; 0.4 -0.8; 0.2 0.4;;; 0.15 0.3; 0.15 0.3; 0.2 0.4; 0.2 0.4; 0.4 -0.8],
  d0 = [72.0 100.0; 100.0 100.0; 50.0 50.0; 50.0 50.0; 10.0 10.0],
  const_s = [100.0 100.0; 100.0 100.0; 200.0 200.0; 200.0 200.0; 10.0 10.0],
  ϵ_s = [0.8 0.0; 0.0 0.0; -0.6 -0.3; 0.0 0.0; -0.6 -0.3;;; 0.0 0.0; 0.8 0.0; 0.0 0.0; -0.6 -0.3; -0.6 -0.3;;; 0.0 0.0; 0.0 0.0; 0.8 0.0; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.8 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.8 0.0;;;; 0.0 0.8; 0.0 0.0; -0.3 -0.6; 0.0 0.0; -0.3 -0.6;;; 0.0 0.0; 0.0 0.0; 0.0 0.0; -0.3 -0.6; -0.3 -0.6;;; 0.0 0.0; 0.0 0.0; 0.0 0.8; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.8; 0.0 0.0;;; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.8],
  s0 = [72.0 100.0; 100.0 100.0; 50.0 50.0; 50.0 50.0; 10.0 10.0],
  pr0 = [0.64 10.0; 10.0 10.0; 30.0 30.0; 30.0 30.0; 100.0 100.0];
  )


  (np, nr) = size(ϵ_d)

  # *** Model declaration ***
  #m = Model(optimizer_with_attributes(Ipopt.Optimizer))
  m = Model(() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe))


  # *** Variables declaration ***
  @variable(m, d[p = 1:np, r = 1:nr] >= 0,  start=d0[p,r])   # Demand quantity
  @variable(m, s[p = 1:np, r = 1:nr] >= 0,  start=s0[p,r])   # Supply quantity
  @variable(m, pr[p = 1:np, r = 1:nr] >= 0,  start=pr0[p,r]) # Price

  # *** Constraints declaration (and definition) ***
  @NLconstraint(m, demand[p in 1:np, r in 1:nr],  d[p,r]  == const_d[p,r] * prod(pr[in_p,in_r]^ ϵ_d[p,r,in_p,in_r] for in_p in 1:np, in_r in 1:nr if ϵ_d[p,r,in_p,in_r] != 0.0  ) )
  @NLconstraint(m, supply[p in 1:np, r in 1:nr],  s[p,r]  == const_s[p,r] * prod(pr[in_p,in_r]^ ϵ_s[p,r,in_p,in_r] for in_p in 1:np, in_r in 1:nr if ϵ_s[p,r,in_p,in_r] != 0.0) )

  # *** Objective function ***
  @NLobjective(m, Min,
      sum((d[p,r] - s[p,r])^2 for p in 1:np, r in 1:nr)
  )

  # *** Print human-readable version of the model ***
  #print(m)

  # *** Solve the model  ***
  optimize!(m)
  status = termination_status(m)

  # *** Printing outputs (optimal quantities/prices) ***
  if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT  || status == MOI.ITERATION_LIMIT) && has_values(m)
    return (solved=true, model=m, status=status, objective_value=objective_value(m), optD=value.(d), optS  = value.(s), optPr = value.(pr))
  else
    return (solved=false, model=m, status=status, objective_value=objective_value(m), optD=zeros(np,nr), optS  = zeros(np,nr), optPr = zeros(np,nr))
  end
end