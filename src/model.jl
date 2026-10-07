# ----------------------------------------------------------------------------
#  Model construction
# ----------------------------------------------------------------------------

# small positive floor so CES powers q^ρ stay differentiable away from 0
const QFLOOR = 1e-4

"""
$(TYPEDSIGNATURES)

Build and solve the equilibrium for the economy `d` and return a
[`Results`](@ref).

A warning is emitted if the solver does not report an optimal (or locally
optimal) solution. `optimizer` is any JuMP solver able to handle nonlinear
constraints and to return duals; `silent` suppresses its output.
"""
function solve_market(d::MarketData; optimizer = Ipopt.Optimizer, silent::Bool = true)
    m = Model(optimizer)
    silent && set_silent(m)
    # Ipopt prints its licence banner even when silent unless "sb" is set
    silent && optimizer === Ipopt.Optimizer && set_attribute(m, "sb", "yes")

    demand_of = Dict((s.region, s.product) => s for s in d.demand)
    supply_of = Dict((s.region, s.product) => s for s in d.supply)

    # --- decision variables -------------------------------------------------
    # final demand quantities
    D = Dict{Tuple{Symbol,Symbol},VariableRef}()
    for s in d.demand
        D[(s.region, s.product)] = @variable(m, lower_bound = QFLOOR, start = s.q0)
    end
    # primary supply quantities
    S = Dict{Tuple{Symbol,Symbol},VariableRef}()
    for s in d.supply
        S[(s.region, s.product)] = @variable(m, lower_bound = QFLOOR, start = s.q0)
    end
    # process activity levels
    z = Dict{Tuple{Symbol,Symbol},VariableRef}()
    for proc in d.processes, r in regions_of(proc, d)
        z[(r, proc.name)] = @variable(m, lower_bound = 0.0, start = 1.0)
    end
    # trade flows  T[(product, from, to)]
    T = Dict{Tuple{Symbol,Symbol,Symbol},VariableRef}()
    for ((p, from, to), _) in d.transport
        T[(p, from, to)] = @variable(m, lower_bound = 0.0, start = 0.0)
    end

    # --- assemble per-(region,product) sources & uses -----------------------
    produced = Dict{Tuple{Symbol,Symbol},Any}()   # process outputs + primary supply
    used     = Dict{Tuple{Symbol,Symbol},Any}()   # process inputs + final demand
    addto!(dict, key, term) = (dict[key] = haskey(dict, key) ? dict[key] + term : term)

    for proc in d.processes, r in regions_of(proc, d)
        zk = z[(r, proc.name)]
        for (p, yield) in proc.outputs
            addto!(produced, (r, p), yield * zk)
        end
        for nest in proc.inputs
            if length(nest.products) == 1            # Leontief
                addto!(used, (r, nest.products[1]), nest.composite * zk)
            else                                     # CES
                qs = VariableRef[]
                for p in nest.products
                    q = @variable(m, lower_bound = QFLOOR, start = nest.composite)
                    push!(qs, q)
                    addto!(used, (r, p), q)
                end
                # Calibrated CES: shares are *value* shares (Σδ = 1) and enter
                # the quantity aggregator as δ^(1/σ).  This makes the implied
                # composite-input price equal the input price when prices are
                # uniform (and total physical input = composite·z), so the dual
                # prices reproduce the genuine marginal cost of production.
                ρ = (nest.sigma - 1) / nest.sigma
                aggr = @expression(m, sum(nest.shares[i]^(1 / nest.sigma) * qs[i]^ρ
                                          for i in eachindex(qs)))
                @constraint(m, nest.phi * aggr^(1 / ρ) >= nest.composite * zk)
            end
        end
    end

    # --- material balance per (region, product)  → duals are the prices -----
    balance = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    for r in d.regions, p in d.products
        src = get(produced, (r, p), 0.0)
        use = get(used, (r, p), 0.0)
        haskey(S, (r, p)) && (src += S[(r, p)])
        haskey(D, (r, p)) && (use += D[(r, p)])
        # imports into r minus exports from r
        imp = @expression(m, sum(T[(p, o, r)] for o in d.regions if haskey(T, (p, o, r)); init = 0.0))
        exp = @expression(m, sum(T[(p, r, o)] for o in d.regions if haskey(T, (p, r, o)); init = 0.0))
        balance[(r, p)] = @constraint(m, src + imp - use - exp == 0.0)
    end

    # --- objective: net social surplus --------------------------------------
    #   + gross consumer benefit  Σ ∫₀ᴰ P(d)dd
    #   − resource (supply) cost   Σ ∫₀ˢ P(s)ds
    #   − conversion (value-added) cost
    #   − transport cost
    cons_benefit = @expression(m, sum(
        let s = demand_of[k], a = s.p0 * s.q0^(1 / s.elasticity), e = 1 - 1 / s.elasticity
            (a / e) * D[k]^e
        end for k in keys(D)))
    supply_cost = @expression(m, sum(
        let s = supply_of[k], b = s.p0 * s.q0^(-1 / s.elasticity), e = 1 + 1 / s.elasticity
            (b / e) * S[k]^e
        end for k in keys(S)))
    conv_cost = @expression(m, sum(proc.vacost * z[(r, proc.name)]
                                   for proc in d.processes for r in regions_of(proc, d)))
    trans_cost = @expression(m, sum(d.transport[key] * T[key] for key in keys(T); init = 0.0))

    @objective(m, Max, cons_benefit - supply_cost - conv_cost - trans_cost)

    optimize!(m)
    st = termination_status(m)
    (st == MOI.LOCALLY_SOLVED || st == MOI.OPTIMAL) ||
        @warn "solver returned status $st — results may be unreliable"

    return build_results(d, m, D, S, z, T, balance)
end
