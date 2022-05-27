export solveEquilibrium

using JuMP, Ipopt

function solveEquilibrium(const_d,ϵ_d,d0,const_s,ϵ_s,s0,pr0;functionalForm="constantElasticity",debug=false)

    (np, nr) = size(ϵ_d)

    # *** Model declaration ***
    #m = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>1000, "tol"=>0.001, "linear_solver"=>"mumps", "print_level" =>1))
    m = Model(optimizer_with_attributes(Ipopt.Optimizer))

    # *** Parameters declarations (useful to avoid model reconstruction if these changes)  ***
    #@NLparameter(m, _ϵ_d[p = 1:np, r = 1:nr, in_p = 1:np, in_r = 1:nr ] == ϵ_d[p,r,in_p,in_r])
    #@NLparameter(m, _ϵ_s[p = 1:np, r = 1:nr, in_p = 1:np, in_r = 1:nr ] == ϵ_s[p,r,in_p,in_r])
    #@NLparameter(m, _const_d[p = 1:np, r = 1:nr] == const_d[p,r])
    #@NLparameter(m, _const_s[p = 1:np, r = 1:nr] == const_s[p,r])

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