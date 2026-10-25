# ----------------------------------------------------------------------------
#  Validation
#
#  The engine builds one balance per (region, product) of the economy, so a
#  name that is not in `regions` or `products` gets variables and objective
#  terms but no balance: a mistyped input becomes a free input, a mistyped
#  consumed product becomes an unbounded one, and the solver reports success.
#  Everything that could happen that way is checked here instead.
# ----------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Check that an economy refers only to names it declares, and that nothing in it
is declared twice.  Throws an `ArgumentError` naming the offending entry on the
first problem found; returns `nothing` when the economy is well formed.

[`solve_market`](@ref) calls it before building anything, so it rarely has to
be called directly.  The checks are:

* `regions` and `products` are non-empty and free of duplicates;
* every product and region named anywhere — in a curve, a process, a transport
  route, an [`Armington`](@ref) specification, a [`TradeLimit`](@ref) or a
  [`Tax`](@ref) — is one of those declared;
* no two [`SupplySpec`](@ref)s and no two [`DemandSpec`](@ref)s describe the
  same product in the same region, and no two [`Process`](@ref)es share a name;
* a transport route goes between two different regions;
* a selector of a tax or a trade limit matches something;
* no process produces out of nothing at no cost, which would make the economy
  unbounded.
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
        field == "all" || for x in (field isa AbstractString ? [field] : field)
            x in set || fail("$what selects the $kind $x, which the economy does " *
                             "not declare")
        end

    # --- curves -------------------------------------------------------------
    seen_demand = Set{Tuple{String,String}}()
    for s in d.demand
        check_region(s.region, "a demand curve")
        for p in nest_products(s.product)
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
        proc.regions == "all" || for r in regions_of(proc, d)
            check_region(r, "process $(proc.name)")
        end
        for nest in proc.inputs, p in nest.products
            check_product(p, "an input of process $(proc.name)")
        end
        for (p, _) in proc.outputs
            check_product(p, "an output of process $(proc.name)")
        end
        proc.capacity isa Dict && for r in keys(proc.capacity)
            check_region(r, "the capacity of process $(proc.name)")
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
        # its outputs, and the economy around it has no finite optimum
        isempty(proc.inputs) && maximum_vacost(proc) <= 0 &&
            !isfinite(maximum_capacity(proc)) &&
            fail("process $(proc.name) has no inputs, no value-added cost and no " *
                 "capacity, so it would produce without limit; give it an input, a " *
                 "cost or a capacity")
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
        for o in collect_origins!(String[], a.nests)
            check_region(o, "an origin group of $(a.product)")
        end
        a.scale isa Dict && for r in keys(a.scale)
            check_region(r, "the Armington scale of $(a.product)")
        end
    end

    # --- limits and taxes ---------------------------------------------------
    for limit in d.trade_limits
        check_selector(limit.products, products, "trade limit $(limit.name)", "product")
        check_selector(limit.from, regions, "trade limit $(limit.name)", "region")
        check_selector(limit.to, regions, "trade limit $(limit.name)", "region")
        for p in keys(limit.weights)
            check_product(p, "the weights of trade limit $(limit.name)")
        end
    end
    for t in d.taxes
        what = "tax $(t.name)"
        used, unused = tax_selectors(t.on)
        for (field, value) in ((:product, t.product), (:region, t.region),
                               (:process, t.process), (:from, t.from), (:to, t.to))
            value == "all" && continue
            field in unused && fail("$what falls on $(t.on), which ignores the " *
                                    "`$field` selector; it would apply more widely " *
                                    "than it looks. Use $(join(used, ", ")) instead")
        end
        check_selector(t.product, products, what, "product")
        check_selector(t.region, regions, what, "region")
        check_selector(t.process, Set(names), what, "process")
        check_selector(t.from, regions, what, "region")
        check_selector(t.to, regions, what, "region")
    end
    return nothing
end

"The entries of `xs` that occur more than once, in order of first occurrence."
duplicates(xs) = unique(x for x in xs if count(==(x), xs) > 1)

"The largest value-added cost a process pays anywhere."
maximum_vacost(p::Process) =
    p.vacost isa Dict ? (isempty(p.vacost) ? 0.0 : maximum(values(p.vacost))) : p.vacost

"The largest capacity a process may run at anywhere, `Inf` if unbounded."
maximum_capacity(p::Process) =
    p.capacity isa Dict ? (isempty(p.capacity) ? Inf : maximum(values(p.capacity))) :
                          p.capacity

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
