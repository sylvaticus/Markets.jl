# ----------------------------------------------------------------------------
#  Validation
#
#  The engine builds one balance per (region, product) of the economy, so a
#  name that is not in `regions` or `products` gets variables and objective
#  terms but no balance: a mistyped input becomes a free input, a mistyped
#  consumed product becomes an unbounded one, and the solver reports success.
#  Everything that could happen that way is checked here instead, together
#  with every other way an economy can be inconsistent with itself, so that
#  `validate` alone is the whole verdict on a specification.
# ----------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Check that an economy refers only to names it declares, that nothing in it is
declared twice, and that every part of it can be built as written.  Throws an
`ArgumentError` naming the offending entry on the first problem found; returns
`nothing` when the economy is well formed.

[`solve_market`](@ref) and [`jump_model`](@ref) call it before building
anything, so it rarely has to be called directly.  The checks are:

* `regions` and `products` are non-empty and free of duplicates;
* every product and region named anywhere — in a curve, a process, a transport
  route, an [`Armington`](@ref) specification, a [`Limit`](@ref) or a
  [`Tax`](@ref) — is one of those declared;
* no two [`SupplyCurve`](@ref)s and no two [`DemandCurve`](@ref)s describe the
  same product in the same region, and no two [`Process`](@ref)es share a name;
* a transport route goes between two different regions;
* the [`Armington`](@ref) specifications are consistent: at most one per
  product and destination, and a product imperfectly substitutable everywhere
  or nowhere;
* a selector of a [`Tax`](@ref) or a [`Limit`](@ref) names something the
  economy declares **and** selects at least one quantity of it, no two limits
  share a name, and the members of one demand group carry one ad valorem rate;
* no process produces out of nothing at no cost, which would make the economy
  unbounded.

It also **warns**, without failing, about an [`Armington`](@ref) share given to
an origin that cannot produce the product, or to a pair with no route: the
engine drops such a share, so the pattern of trade differs from the one written
down.
"""
function validate(d::MarketData)
    fail(msg) = throw(ArgumentError(msg))

    # --- the two index sets everything else is checked against --------------
    isempty(d.regions)  && fail("the economy declares no regions")
    isempty(d.products) && fail("the economy declares no products")
    allunique(d.regions)  || fail("these regions are declared more than once: " *
                                  "$(join(duplicates(d.regions), ", "))")
    allunique(d.products) || fail("these products are declared more than once: " *
                                  "$(join(duplicates(d.products), ", "))")
    regions, products = Set(d.regions), Set(d.products)

    check_product(p, what) = p in products ||
        fail("$what names the product $p, which the economy does not declare")
    check_region(r, what) = r in regions ||
        fail("$what names the region $r, which the economy does not declare")
    check_selector(field, set, what, kind) =
        field isa Colon || for x in each(field)
            x in set || fail("$what selects the $kind $x, which the economy does " *
                             "not declare")
        end

    # --- curves -------------------------------------------------------------
    seen_demand = Set{Tuple{String,String}}()
    for s in d.demand
        check_region(s.region, "a demand curve")
        for p in leaves(s.product)
            check_product(p, "the demand curve in $(s.region)")
            (p, s.region) in seen_demand &&
                fail("$p has more than one demand curve in region $(s.region)")
            push!(seen_demand, (p, s.region))
        end
    end
    seen_supply = Set{Tuple{String,String}}()
    for s in d.supply
        check_region(s.region, "a supply curve")
        check_product(s.product, "the supply curve in $(s.region)")
        (s.product, s.region) in seen_supply &&
            fail("$(s.product) has more than one supply curve in region $(s.region)")
        push!(seen_supply, (s.product, s.region))
    end

    # --- processes ----------------------------------------------------------
    names = [p.name for p in d.processes]
    allunique(names) || fail("these processes are declared more than once: " *
                             "$(join(duplicates(names), ", "))")
    for proc in d.processes
        proc.regions isa Colon || for r in regions_of(proc, d)
            check_region(r, "process $(proc.name)")
        end
        for nest in proc.inputs, p in nest.products
            check_product(p, "an input of process $(proc.name)")
        end
        for (p, _) in proc.outputs
            check_product(p, "an output of process $(proc.name)")
        end
        if proc.vacost isa Dict
            for r in keys(proc.vacost)
                check_region(r, "the value-added cost of process $(proc.name)")
            end
            # a region left out would run the process for free, which is
            # never what was meant
            missed = [r for r in regions_of(proc, d) if !haskey(proc.vacost, r)]
            isempty(missed) || fail("process $(proc.name) operates in " *
                "$(join(missed, ", ")) but its value-added cost gives none there")
        end
        # a process with nothing to buy and nothing to pay is a free source of
        # its outputs, and the economy around it has no finite optimum unless a
        # limit holds its activity back — a ceiling on it, not a floor under it
        held = any((isfinite(l.upper) && w > 0 || isfinite(l.lower) && w < 0) &&
                   q.kind == :activity && matches(q.selection.process, proc.name)
                   for l in d.limits for (q, w) in l.expression.terms)
        isempty(proc.inputs) && maximum_vacost(proc) <= 0 && !held &&
            fail("process $(proc.name) has no inputs, no value-added cost and no " *
                 "limit on its activity, so it would produce without limit; give it " *
                 "an input, a cost or a limit")
    end

    # --- trade --------------------------------------------------------------
    for (p, o, r) in keys(d.transport_costs)
        check_product(p, "a transport route")
        check_region(o, "a transport route of $p")
        check_region(r, "a transport route of $p")
        o === r && fail("the transport route of $p from $o to $r goes nowhere")
    end

    # --- Armington ----------------------------------------------------------
    for a in d.armington
        check_product(a.product, "an Armington specification")
        check_selector(a.destination, regions, "the Armington specification of " *
                       "$(a.product)", "region")
        for (o, r) in keys(a.shares)
            check_region(o, "the Armington shares of $(a.product)")
            check_region(r, "the Armington shares of $(a.product)")
        end
        for o in leaves(a.nests)
            check_region(o, "an origin group of $(a.product)")
        end
        a.scale isa Dict && for r in keys(a.scale)
            check_region(r, "the Armington scale of $(a.product)")
        end
    end
    specs = armington_of(d) # one specification per destination, Armington everywhere or nowhere
    # A share on an origin that cannot produce the product, or on a pair with
    # no route, falls on a variety the model will not have, and is dropped. A
    # share matrix written for every pair is a reasonable way to give the
    # shares, so this is not an error; but the pattern of trade would then
    # differ silently from the one written down, so it is said.
    let can = producible(d), uses = consumed(d), dropped = Dict{String,Vector{String}}()
        for ((p, r), a) in specs
            (p, r) in uses || continue
            for ((o, dest), s) in a.shares
                (dest == r && s > 0) || continue
                why = !((p, o) in can) ? "$o produces no $p" :
                      o != r && !haskey(d.transport_costs, (p, o, r)) ? "no route from $o to $r" :
                      nothing
                why === nothing || push!(get!(dropped, p, String[]), "$o → $r ($why)")
            end
        end
        for p in sort!(collect(keys(dropped)))
            @warn "Armington shares of $p fall on varieties the model will not have, " *
                  "and are dropped: $(join(sort!(dropped[p]), ", "))"
        end
    end

    # --- limits and taxes ---------------------------------------------------
    limit_names = [l.name for l in d.limits]
    allunique(limit_names) || fail("these limits are declared more than once: " *
                                   "$(join(duplicates(limit_names), ", "))")
    fields = Dict(:region => (regions, "region"), :from => (regions, "region"),
                  :to => (regions, "region"), :product => (products, "product"),
                  :process => (Set(names), "process"))
    for limit in d.limits, (q, _) in limit.expression.terms, (f, value) in pairs(q.selection)
        set, kind = fields[f]
        check_selector(value, set, "limit $(limit.name)", kind)
    end
    # a selector may name only declared things and still select no quantity
    # the model has: a supply the economy does not curve, a flow with no route
    isempty(d.limits) || let cat = model_catalogue(d)
        foreach(limit -> limit_terms(limit, d, cat), d.limits)
    end
    for t in d.taxes
        what = "tax $(t.name)"
        used, unused = tax_selectors(t.on)
        for (field, value) in ((:product, t.product), (:region, t.region),
                               (:process, t.process), (:from, t.from), (:to, t.to))
            value isa Colon && continue
            field in unused && fail("$what falls on $(t.on), which ignores the " *
                                    "`$field` selector; it would apply more widely " *
                                    "than it looks. Use $(join(used, ", ")) instead")
        end
        check_selector(t.product, products, what, "product")
        check_selector(t.region, regions, what, "region")
        check_selector(t.process, Set(names), what, "process")
        check_selector(t.from, regions, what, "region")
        check_selector(t.to, regions, what, "region")
        bites = t.on == "supply"  ? any(matches(t.product, s.product) && matches(t.region, s.region)
                                        for s in d.supply) :
                t.on == "demand"  ? any(matches(t.product, p) && matches(t.region, s.region)
                                        for s in d.demand for p in leaves(s.product)) :
                t.on == "process" ? any(matches(t.process, proc.name) && matches(t.region, r)
                                        for proc in d.processes for r in regions_of(proc, d)) :
                                    any(matches(t.product, p) && matches(t.from, o) && matches(t.to, r)
                                        for (p, o, r) in keys(d.transport_costs))
        bites || fail("$what falls on $(t.on) but selects no $(t.on == "trade" ? "route" :
                      t.on == "process" ? "process" : "$(t.on) curve") of this economy")
    end
    # several relative rates on one base must not add up to -1 or below, which
    # `tax_rates` checks when asked for the rates of that base: ask for every
    # base the model taxes (the demand ones are asked for just below)
    if !isempty(d.taxes)
        foreach(s -> tax_rates(d.taxes, "supply"; product = s.product, region = s.region),
                d.supply)
        for proc in d.processes, r in regions_of(proc, d)
            tax_rates(d.taxes, "process"; process = proc.name, region = r)
        end
        for (p, o, r) in keys(d.transport_costs)
            tax_rates(d.taxes, "trade"; product = p, from = o, to = r)
        end
    end
    # the members of one demand group share one curve, hence one ad valorem rate
    foreach(s -> demand_relative_tax(d, s), d.demand)
    return nothing
end

"The entries of `xs` that occur more than once, in order of first occurrence."
duplicates(xs) = unique(x for x in xs if count(==(x), xs) > 1)

"The largest value-added cost a process pays anywhere."
maximum_vacost(p::Process) =
    p.vacost isa Dict ? (isempty(p.vacost) ? 0.0 : maximum(values(p.vacost))) : p.vacost

"""
$(TYPEDSIGNATURES)

Which selector fields of a [`Tax`](@ref) the base `on` uses, and which it
ignores.  A tax that sets a selector its base ignores would apply more widely
than it appears to, so [`validate`](@ref) rejects it.
"""
function tax_selectors(on::AbstractString)
    on == "supply"  && return ((:product, :region), (:process, :from, :to))
    on == "demand"  && return ((:product, :region), (:process, :from, :to))
    on == "process" && return ((:process, :region), (:product, :from, :to))
    return ((:product, :from, :to), (:region, :process))          # "trade"
end
