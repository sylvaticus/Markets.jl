# ----------------------------------------------------------------------------
#  Displaying an economy
#
#  A `MarketData` holds its curves, recipes and routes as vectors and a
#  dictionary, which the default display prints as one unreadable line. What
#  is shown instead is what the economy *is*: which products are supplied,
#  consumed, made and traded, how the processes chain into one another, where
#  goods may move and at what cost, and which taxes and limits bite.
# ----------------------------------------------------------------------------

# A number as briefly as it can be read: integers without a decimal point.
fmt(x::Real) = isinf(x) ? (x > 0 ? "∞" : "-∞") :
               isinteger(x) && abs(x) < 1e15 ? string(Int(x)) :
               string(round(x; sigdigits = 4))

plural(n, one, many = one * "s") = string(n, " ", n == 1 ? one : many)

# A selector field as text: its names, or `all`.
listing(sel) = sel isa AbstractString ? sel : join(sel, ", ")

# `coef name`, leaving out a coefficient of 1.
amount(coef, name) = coef == 1 ? name : string(fmt(coef), " ", name)

pad(s, n)  = s * " "^max(0, n - textwidth(s))
lpad_(s, n) = " "^max(0, n - textwidth(s)) * s

# `s` cut to `width` columns, with an ellipsis where it was cut.
function clip(s::AbstractString, width::Integer)
    textwidth(s) <= width && return s
    out, used = IOBuffer(), 0
    for c in s
        used + textwidth(c) > width - 1 && break
        print(out, c)
        used += textwidth(c)
    end
    return String(take!(out)) * "…"
end

# `rows` as aligned columns, indented by two spaces. Columns listed in `right`
# are right-aligned, which is what numbers want.
function columns(rows::AbstractVector{<:AbstractVector{<:AbstractString}}; right = ())
    isempty(rows) && return String[]
    widths = [maximum(textwidth(r[j]) for r in rows) for j in eachindex(first(rows))]
    cell(j, c) = j in right ? lpad_(c, widths[j]) : pad(c, widths[j])
    return [rstrip("  " * join([cell(j, c) for (j, c) in enumerate(r)], "  ")) for r in rows]
end

nest_label(n::Nest) = length(n.products) == 1 ? amount(n.composite, only(n.products)) :
    amount(n.composite, "CES(" * join(n.products, ", ") * "; σ " * fmt(n.sigma) * ")")

group_label(n::DemandNest) =
    "{" * join(vcat(n.products, group_label.(n.nests)), ", ") * "} σ " * fmt(n.sigma)
group_label(n::OriginNest) =
    "{" * join(vcat(n.origins, group_label.(n.nests)), ", ") * "} σ " * fmt(n.sigma)

# The stage of each process in `d`: 1 for a process whose inputs are all primary
# products, and one more than the latest of its inputs otherwise, an input being
# available from the earliest process that makes it. A process caught in a cycle
# that nothing primary feeds has no stage, and gets `typemax(Int)`.
function stages(d::MarketData)
    never = typemax(Int)
    makers = Dict{String,Vector{Process}}()
    for proc in d.processes, (p, _) in proc.outputs
        push!(get!(makers, p, Process[]), proc)
    end
    supplied = Set(s.product for s in d.supply)
    available(p) = p in supplied || !haskey(makers, p) ? 0 : never
    product_stage = Dict{String,Int}()
    for proc in d.processes, n in proc.inputs, p in n.products
        product_stage[p] = available(p)
    end
    stage = Dict(proc.name => never for proc in d.processes)
    # relax until nothing moves: at most one pass per process in a chain
    for _ in 0:length(d.processes)
        changed = false
        for proc in d.processes
            latest = maximum((product_stage[p] for n in proc.inputs for p in n.products);
                             init = 0)
            s = latest == never ? never : latest + 1
            s < stage[proc.name] && (stage[proc.name] = s; changed = true)
            for (p, _) in proc.outputs
                haskey(product_stage, p) && s < product_stage[p] &&
                    (product_stage[p] = s; changed = true)
            end
        end
        changed || break
    end
    return stage
end

# The routes of `product`, in the order its regions are declared, with a pair of
# routes that cost the same both ways shown once, as `A⇄B`.
function route_labels(d::MarketData, product::String)
    costs = Dict((o, r) => c for ((p, o, r), c) in d.transport_costs if p == product)
    rank(r) = something(findfirst(==(r), d.regions), length(d.regions) + 1)
    labels, done = String[], Set{Tuple{String,String}}()
    for (o, r) in sort(collect(keys(costs)); by = k -> (rank(k[1]), rank(k[2]), k))
        (o, r) in done && continue
        c = costs[(o, r)]
        if get(costs, (r, o), nothing) == c
            push!(labels, string(o, "⇄", r, " ", fmt(c)))
            push!(done, (r, o))
        else
            push!(labels, string(o, "→", r, " ", fmt(c)))
        end
        push!(done, (o, r))
    end
    return labels, costs
end

# The `from` and `to` selectors of a tax as text: `A → B`, `from A` or `to B`,
# and empty when both select every region.
function direction(from, to)
    from isa Colon && to isa Colon && return ""
    from isa Colon && return "to " * listing(to)
    to isa Colon && return "from " * listing(from)
    return string(listing(from), " → ", listing(to))
end

const TAX_BASE = Dict("supply" => "of the producer price", "demand" => "of the market price",
                      "process" => "of the value-added cost", "trade" => "of the origin price")

function tax_row(t::Tax)
    kind = t.rate < 0 ? "subsidy" : "tax"
    rate = t.relative ? string(fmt(100 * abs(t.rate)), "% ", TAX_BASE[t.on]) :
                        string(fmt(abs(t.rate)), " per unit")
    scope = if t.on == "trade"
        join(filter(!isempty, [t.product isa Colon ? "" : listing(t.product),
                               direction(t.from, t.to)]), " ")
    else
        what = t.on == "process" ? t.process : t.product
        join(filter(!isempty, [what isa Colon ? "" : listing(what),
                               t.region isa Colon ? "" : "in " * listing(t.region)]), " ")
    end
    return [t.name, string(kind, " on ", t.on), rate, isempty(scope) ? "everywhere" : scope]
end

# The economy in one line: what the REPL shows when it is typed, what `repr`
# gives, and what stands for it inside a vector, a dictionary or another
# structure. `print` gives the summary below.
Base.show(io::IO, d::MarketData) = print(io, "MarketData(",
    plural(length(d.regions), "region"), ", ", plural(length(d.products), "product"), ", ",
    plural(length(d.processes), "process", "processes"), ")")

"""
$(TYPEDSIGNATURES)

A summary of the economy `d`: its products and what happens to each, the
processes in the order they chain, the routes and what they cost, and the
taxes, subsidies and limits. Shown by `print(d)` and `println(d)`, hence by
`string(d)` and `"\$d"`, while typing `d` in the REPL, `repr(d)` and a container
holding it give only its one-line count; printed to a terminal, lines wider
than it are cut, and a route list too long for one line is summarised.
"""
function Base.print(io::IO, d::MarketData)
    limit = get(io, :limit, false) || io isa Base.TTY
    width = limit ? displaysize(io)[2] : typemax(Int)
    line(s) = println(io, clip(s, width))
    regions_in(rs) = rs isa Colon || Set(rs) == Set(d.regions) ? "every region" : listing(rs)
    by_region(dict) = sort(collect(dict);
        by = ((r, _),) -> (something(findfirst(==(r), d.regions), length(d.regions) + 1), r))
    heading(s) = (println(io); printstyled(io, clip(s, width); bold = true); println(io))

    printstyled(io, "MarketData"; bold = true)
    println(io, clip(": " * join([plural(length(d.regions), "region"),
                            plural(length(d.products), "product"),
                            plural(length(d.processes), "process", "processes"),
                            plural(length(d.transport_costs), "trade route"),
                            plural(length(d.taxes), "tax", "taxes")], ", "),
                     width - textwidth("MarketData")))
    line("  regions: " * join(d.regions, ", "))

    # --- products -----------------------------------------------------------
    count_of(pairs) = (c = Dict{String,Int}(); for k in pairs; c[k] = get(c, k, 0) + 1; end; c)
    supplied = count_of(s.product for s in d.supply)
    consumed = count_of(p for s in d.demand for p in leaves(s.product))
    made     = count_of(p for proc in d.processes for (p, _) in proc.outputs)
    used     = count_of(p for proc in d.processes for n in proc.inputs for p in unique(n.products))
    routes   = count_of(k[1] for k in keys(d.transport_costs))
    varieties = Dict{String,Vector{String}}()
    for a in d.armington
        push!(get!(varieties, a.product, String[]), fmt(a.sigma))
    end
    shown(c, p) = haskey(c, p) ? string(c[p]) : "·"
    heading("Products")
    header = ["", "supplied in", "consumed in", "made by", "used by", "routes"]
    isempty(varieties) || push!(header, "Armington σ")
    rows = [header]
    for p in d.products
        row = [p, shown(supplied, p), shown(consumed, p), shown(made, p), shown(used, p),
               shown(routes, p)]
        isempty(varieties) || push!(row, join(unique(get(varieties, p, ["·"])), ", "))
        push!(rows, row)
    end
    foreach(line, columns(rows; right = 2:7))
    line("  (regions supplying and consuming it, processes making and using it)")

    groups = Dict{String,Vector{String}}()
    for s in d.demand
        s.product isa DemandNest || continue
        push!(get!(groups, group_label(s.product), String[]), s.region)
    end
    if !isempty(groups)
        heading("Demand groups")
        foreach(line, columns([[g, "in " * regions_in(r)] for (g, r) in sort(collect(groups))]))
    end

    # --- transformations ----------------------------------------------------
    if !isempty(d.processes)
        stage = stages(d)
        order = sort(d.processes; by = p -> (stage[p.name], findfirst(==(p), d.processes)))
        multistage = length(unique(values(stage))) > 1
        heading(multistage ? "Transformations, by stage" : "Transformations")
        inputs = [join(nest_label.(p.inputs), " + ") for p in order]
        reach = maximum(textwidth, inputs)
        rows, notes = Vector{String}[], String[]
        for (p, inp) in zip(order, inputs)
            note = String[]
            if p.vacost isa Dict
                # one per region would crowd the chart: the range says enough
                lo, hi = isempty(p.vacost) ? (0.0, 0.0) : extrema(values(p.vacost))
                hi == 0 || push!(note, "cost " * (lo == hi ? fmt(lo) : fmt(lo) * "–" * fmt(hi)))
            else
                p.vacost == 0 || push!(note, "cost " * fmt(p.vacost))
            end
            p.regions isa Colon || push!(note, "in " * regions_in(p.regions))
            push!(notes, join(note, " · "))
            out = isempty(p.outputs) ? "(disposed of)" :
                  join((amount(y, q) for (q, y) in p.outputs), " + ")
            row = [p.name, string(inp, " ", "─"^(reach - textwidth(inp) + 1), "→ ", out)]
            s = stage[p.name]
            multistage && pushfirst!(row, s == typemax(Int) ? "?" : string(s))
            push!(rows, row)
        end
        chart = columns(rows; right = multistage ? (1:1) : ())
        # the notes form a last column where the terminal is wide enough for
        # them, and go on a line of their own, under the inputs, where it is not
        at = maximum(textwidth, chart) + 2
        under = 2 + sum(maximum(textwidth(r[j]) for r in rows) + 2 for j in 1:length(rows[1]) - 1)
        for (c, n) in zip(chart, notes)
            if isempty(n)
                line(c)
            elseif at + textwidth(n) <= width
                line(pad(c, at) * n)
            else
                line(c)
                line(" "^under * n)
            end
        end
    end

    # --- trade --------------------------------------------------------------
    traded = [p for p in d.products if haskey(routes, p)]
    if !isempty(traded)
        heading("Trade routes and their cost per unit")
        name_width = maximum(textwidth, traded)
        for p in traded
            labels, costs = route_labels(d, p)
            full = string("  ", pad(p, name_width), "  ", join(labels, ", "))
            if textwidth(full) > width
                ends = unique(vcat(first.(keys(costs)), last.(keys(costs))))
                n = length(ends)
                among = length(costs) == n * (n - 1) ?
                    string("every route among ", plural(n, "region")) :
                    string(plural(length(costs), "route"), " among ", plural(n, "region"))
                full = string("  ", pad(p, name_width), "  ", among, ", costing ",
                              fmt(minimum(values(costs))), "–", fmt(maximum(values(costs))))
            end
            line(full)
        end
    end
    if !isempty(d.armington)
        heading("Imperfect substitutes (Armington)")
        foreach(line, columns([[a.product, "σ " * fmt(a.sigma),
                                a.destination isa Colon ? "everywhere" : "in " * listing(a.destination),
                                join(filter(!isempty,
                                     [isempty(a.nests) ? "" : "groups " * join(group_label.(a.nests), ", "),
                                      isempty(a.shares) ? "" : "calibrated shares"]), " · ")]
                               for a in d.armington]))
    end

    if !isempty(d.limits)
        heading("Limits")
        foreach(line, columns([[l.name,
                                l.lower == l.upper ? string("= ", fmt(l.upper)) :
                                    join(filter(!isempty,
                                         [isfinite(l.lower) ? "≥ " * fmt(l.lower) : "",
                                          isfinite(l.upper) ? "≤ " * fmt(l.upper) : ""]), ", ") *
                                    (isfinite(l.penalty) ? ", or " * fmt(l.penalty) * " a unit past it" : ""),
                                sprint(show, l.expression)]
                               for l in d.limits]))
    end

    # --- policy -------------------------------------------------------------
    if !isempty(d.taxes)
        heading("Taxes and subsidies")
        foreach(line, columns(tax_row.(d.taxes)))
    end
    return nothing
end
