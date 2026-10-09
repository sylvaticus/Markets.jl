# ----------------------------------------------------------------------------
#  Model construction
# ----------------------------------------------------------------------------

# small positive floor so CES powers q^ρ stay differentiable away from 0
const QFLOOR = 1e-4

"""
$(TYPEDSIGNATURES)

The products of `d` that have a finite Armington elasticity, i.e. whose
regional varieties are imperfect substitutes, as a `product => `[`Armington`](@ref)
dictionary.  Products absent from it are homogeneous.
"""
function armington_of(d::MarketData)
    arm = Dict{Symbol,Armington}()
    for a in d.armington
        isfinite(a.sigma) || continue        # σ = Inf ⇒ homogeneous: nothing to do
        a.product in d.products ||
            throw(ArgumentError("Armington specification for unknown product $(a.product)"))
        haskey(arm, a.product) &&
            throw(ArgumentError("duplicate Armington specification for $(a.product)"))
        arm[a.product] = a
    end
    return arm
end

"""
$(TYPEDSIGNATURES)

The origins that destination `r` can buy product `a.product` from — itself,
plus every origin with a transport route into `r` — together with their
normalised value shares.  Origins without a positive share are left out.
"""
function armington_origins(d::MarketData, a::Armington, r::Symbol)
    p = a.product
    available = [r; [o for o in d.regions if o != r && haskey(d.transport, (p, o, r))]]
    if isempty(a.shares)
        return available, fill(1 / length(available), length(available))
    end
    origins = Symbol[]
    shares  = Float64[]
    for o in available
        s = get(a.shares, (o, r), 0.0)
        s > 0 || continue
        push!(origins, o)
        push!(shares, s)
    end
    isempty(origins) && throw(ArgumentError(
        "no origin with a positive Armington share supplies $p in region $r"))
    shares ./= sum(shares)
    return origins, shares
end

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
    # trade flows of the homogeneous products  T[(product, from, to)]
    arm = armington_of(d)
    T = Dict{Tuple{Symbol,Symbol,Symbol},VariableRef}()
    for ((p, from, to), _) in d.transport
        haskey(arm, p) && continue       # Armington products use X instead
        T[(p, from, to)] = @variable(m, lower_bound = 0.0, start = 0.0)
    end
    # Armington products: quantity of the variety of each origin used in each
    # destination, X[(product, origin, destination)], and the composite that
    # the destination actually uses, A[(product, destination)]
    X = Dict{Tuple{Symbol,Symbol,Symbol},VariableRef}()
    A = Dict{Tuple{Symbol,Symbol},VariableRef}()
    origins_of = Dict{Tuple{Symbol,Symbol},Tuple{Vector{Symbol},Vector{Float64}}}()
    for (p, a) in arm, r in d.regions
        origins, shares = armington_origins(d, a, r)
        origins_of[(p, r)] = (origins, shares)
        # a rough but useful starting point: the local reference consumption
        a0 = get(demand_of, (r, p), nothing)
        start = a0 === nothing ? 1.0 : a0.q0
        A[(p, r)] = @variable(m, lower_bound = QFLOOR, start = start)
        for (i, o) in pairs(origins)
            X[(p, o, r)] = @variable(m, lower_bound = QFLOOR, start = shares[i] * start)
        end
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

    # --- balances per (region, product)  → duals are the prices -------------
    # `balance` is what a region's users face (its dual is the price they pay),
    # `origin_balance` is what its producers sell (its dual is the price they
    # get).  For a homogeneous product the two coincide and only `balance` is
    # built; for an Armington product they differ by the composition of the
    # CES bundle, and both are built.
    balance        = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    origin_balance = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    for r in d.regions, p in d.products
        src = get(produced, (r, p), 0.0)
        use = get(used, (r, p), 0.0)
        haskey(S, (r, p)) && (src += S[(r, p)])
        haskey(D, (r, p)) && (use += D[(r, p)])

        if haskey(arm, p)
            # everything region r produces goes to one of the destinations that
            # buy its variety (itself included)
            ships = @expression(m, sum(X[(p, r, dst)] for dst in d.regions
                                       if haskey(X, (p, r, dst)); init = 0.0))
            origin_balance[(r, p)] = @constraint(m, src - ships == 0.0)
            # the composite covers the local uses
            balance[(r, p)] = @constraint(m, A[(p, r)] - use == 0.0)
            # ... and is the CES aggregate of the varieties bought
            origins, shares = origins_of[(p, r)]
            if length(origins) == 1
                @constraint(m, X[(p, origins[1], r)] >= A[(p, r)])
            else
                σ = arm[p].sigma
                ρ = (σ - 1) / σ
                aggr = @expression(m, sum(shares[i]^(1 / σ) * X[(p, origins[i], r)]^ρ
                                          for i in eachindex(origins)))
                @constraint(m, aggr^(1 / ρ) >= A[(p, r)])
            end
        else
            # homogeneous product: imports into r minus exports from r
            imp = @expression(m, sum(T[(p, o, r)] for o in d.regions if haskey(T, (p, o, r)); init = 0.0))
            exp = @expression(m, sum(T[(p, r, o)] for o in d.regions if haskey(T, (p, r, o)); init = 0.0))
            balance[(r, p)] = @constraint(m, src + imp - use - exp == 0.0)
        end
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
    trans_cost = @expression(m,
        sum(d.transport[key] * T[key] for key in keys(T); init = 0.0) +
        sum(d.transport[key] * X[key] for key in keys(X) if key[2] != key[3]; init = 0.0))

    @objective(m, Max, cons_benefit - supply_cost - conv_cost - trans_cost)

    optimize!(m)
    st = termination_status(m)
    (st == MOI.LOCALLY_SOLVED || st == MOI.OPTIMAL) ||
        @warn "solver returned status $st — results may be unreliable"

    return build_results(d, m, D, S, z, T, X, balance, origin_balance)
end
