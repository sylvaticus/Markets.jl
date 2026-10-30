# Part of Markets.jl package [https://github.com/sylvaticus/Markets.jl] - MIT Licence
# validated 20261030 Antonello Lobianco - still to double check well the two @limit and @limits macros TODO

# ----------------------------------------------------------------------------
#  Limits: linear constraints on the quantities of the model
# ----------------------------------------------------------------------------

"""
$(TYPEDEF)

Quantities of the model of one kind — supply, consumption, a trade flow, … —
summed over the regions, products and processes their `selection` selects.
Built by the term builders, such as `Markets.supply`; the building block of a
[`LimitExpression`](@ref).
"""
struct QuantityRef
    kind::Symbol
    selection::NamedTuple
end
Base.:(==)(a::QuantityRef, b::QuantityRef) = a.kind == b.kind && a.selection == b.selection
Base.hash(q::QuantityRef, h::UInt) = hash(q.selection, hash(q.kind, h))

"""
$(TYPEDEF)

A linear expression in the quantities of the model: a weighted sum of
[`QuantityRef`](@ref)s plus a constant. It is what a [`Limit`](@ref) bounds,
and it is built by adding, subtracting and scaling the term builders —
`Markets.supply`, `Markets.consumption`, `Markets.production`,
`Markets.activity`, `Markets.inputs`, `Markets.flow`, `Markets.exports` and
`Markets.imports` — by numbers. Nothing else is defined on it, so whatever can
be built is linear: the product of two quantities, a quantity in a
denominator, a power or any other function of one is refused.
"""
struct LimitExpression
    terms::Vector{Pair{QuantityRef,Float64}}
    constant::Float64
end
LimitExpression(q::QuantityRef) = LimitExpression([q => 1.0], 0.0)

# --- the term builders --------------------------------------------------------
# Every selector is a name, a vector of names or `:` for every one, as in the
# rest of the package. They are not exported: their names are those a user would give to
# his own variables, and `@limit` finds them anyway.

ref(kind; selection...) = LimitExpression(QuantityRef(kind, map(selector, NamedTuple(selection))))

"""
    Markets.supply(; region = :, product = :)

Primary supply of the selected products in the selected regions, summed.
"""
supply(; region = (:), product = (:)) = ref(:supply; region, product)

"""
    Markets.consumption(; region = :, product = :)

Final demand of the selected products in the selected regions, summed. For an
[`Armington`](@ref) product it is measured, as the demand curves measure it, in
units of the region's composite.
"""
consumption(; region = (:), product = (:)) = ref(:consumption; region, product)

"""
    Markets.production(; region = :, product = :)

Primary supply plus the output of every process, of the selected products in
the selected regions, summed.
"""
production(; region = (:), product = (:)) = ref(:production; region, product)

"""
    Markets.activity(; region = :, process = :)

Activity level of the selected processes in the selected regions, summed.
"""
activity(; region = (:), process = (:)) = ref(:activity; region, process)

"""
    Markets.inputs(; region = :, process = :, product = :)

What the selected processes draw of the selected products in the selected
regions, summed — for an [`Armington`](@ref) product, in units of the region's
composite.
"""
inputs(; region = (:), process = (:), product = (:)) = ref(:inputs; region, process, product)

"""
    Markets.flow(; product = :, from = :, to = :)

Shipments of the selected products between the selected regions, summed. A
region's deliveries to itself are not flows.
"""
flow(; product = (:), from = (:), to = (:)) = ref(:flow; product, from, to)

"""
    Markets.exports(; region = :, product = :)

Shipments of the selected products leaving the selected regions, summed. For a
group of regions, only what leaves the group counts.
"""
exports(; region = (:), product = (:)) = ref(:exports; region, product)

"""
    Markets.imports(; region = :, product = :)

Shipments of the selected products entering the selected regions, summed. For
a group of regions, only what enters the group counts.
"""
imports(; region = (:), product = (:)) = ref(:imports; region, product)

const QUANTITY_BUILDERS = (:supply, :consumption, :production, :activity, :inputs,
                           :flow, :exports, :imports)

# --- arithmetic: only what keeps an expression linear --------------------------

function combine_terms(a::Vector{Pair{QuantityRef,Float64}}, b::Vector{Pair{QuantityRef,Float64}})
    out = copy(a)
    for (q, w) in b
        i = findfirst(t -> first(t) == q, out)
        i === nothing ? push!(out, q => w) : (out[i] = q => last(out[i]) + w)
    end
    return out
end
scaled(e::LimitExpression, a::Real) =
    LimitExpression([q => a * w for (q, w) in e.terms], a * e.constant)

Base.:+(a::LimitExpression, b::LimitExpression) =
    LimitExpression(combine_terms(a.terms, b.terms), a.constant + b.constant)
Base.:+(a::LimitExpression, b::Real) = LimitExpression(a.terms, a.constant + b)
Base.:+(a::Real, b::LimitExpression) = b + a
Base.:+(a::LimitExpression) = a
Base.:-(a::LimitExpression) = scaled(a, -1)
Base.:-(a::LimitExpression, b::LimitExpression) = a + scaled(b, -1)
Base.:-(a::LimitExpression, b::Real) = a + (-b)
Base.:-(a::Real, b::LimitExpression) = a + scaled(b, -1)
Base.:*(a::Real, b::LimitExpression) = scaled(b, a)
Base.:*(a::LimitExpression, b::Real) = scaled(a, b)
Base.:/(a::LimitExpression, b::Real) = scaled(a, 1 / b)
Base.zero(::Type{LimitExpression}) = LimitExpression(Pair{QuantityRef,Float64}[], 0.0)
Base.zero(::LimitExpression) = zero(LimitExpression)

# The tempting mistakes, refused with a reason rather than a MethodError.
not_linear(what) = throw(ArgumentError(
    "a limit must be linear in the quantities of the model, and $what is not"))
Base.:*(::LimitExpression, ::LimitExpression) = not_linear("the product of two quantities")
Base.:/(::LimitExpression, ::LimitExpression) = not_linear("the ratio of two quantities")
Base.:/(::Real, ::LimitExpression) = not_linear("a quantity in a denominator")
Base.:^(::LimitExpression, ::Integer) = not_linear("a power of a quantity")
Base.:^(::LimitExpression, ::Real) = not_linear("a power of a quantity")
Base.inv(::LimitExpression) = not_linear("the inverse of a quantity")
for f in (:abs, :sqrt, :exp, :log)
    @eval Base.$f(::LimitExpression) = not_linear(string($(string(f)), " of a quantity"))
end

as_expression(e::LimitExpression) = e
as_expression(x::Real) = LimitExpression(Pair{QuantityRef,Float64}[], float(x))
as_expression(x) = throw(ArgumentError(
    "a limit bounds a linear expression in the quantities of the model, got $(repr(x))"))

"""
$(TYPEDEF)

A linear constraint on the quantities of the model:

    lower ≤ expression ≤ upper

where `expression` is a [`LimitExpression`](@ref) — supply, consumption,
production, activity, process inputs and trade flows, weighted by any numbers
and summed. It can cap what several products draw from a shared resource, set
a quota or a minimum, tie one quantity to a share of another, or fix a sum.
With `lower == upper` it is an equality.

The limit is part of the economy: it is checked by [`validate`](@ref), shown in
its summary, and reported in `Results.limits` with its value and its `rent` —
the gain in the objective from relaxing the bound that binds by one unit,
positive for an upper bound, negative for a lower one: what a permit to one more
unit would be worth, or the subsidy that would get agents to meet the minimum on
their own.

Write it with [`@limit`](@ref), or build the expression with the term builders,
qualified as `Markets.supply` and so on outside the macro.

# Fields
$(TYPEDFIELDS)

Every ceiling of the model is one: the capacity of a plant, an allowable cut,
what a port can ship. A limit that only caps one variable — one supply, one
activity, one flow — is set as a bound on it, which costs the solver nothing;
any other is a constraint of the model.

A limit is **hard** by default: no price can buy past it. With a finite
`penalty` it is **soft**: it can be crossed, at that cost per unit beyond the
bound, the way a sector model holds trade near its recent levels without
forbidding a larger move. Below the bound nothing changes; past it the rent
is exactly the penalty, which is also the most it can ever be. The penalty is
a device of the model, not a cost anyone pays, so the results report it apart,
in `Results.penalties`, and leave it out of `welfare`; a cost that is really
paid past a threshold — overtime, an emergency import — belongs in the
economy instead, as a second process or route at the higher cost.

# Example
```julia
## an allowable cut, and the capacity of the mills, region by region
@limit "cut SEF" supply(region = "SEF", product = "logs") <= 3.5
[@limit("mill \$r", activity(region = r, process = "sawmill") <= cap[r]) for r in keys(cap)]

## a port that ships at most 6 thousand containers, counted by product
@limit "port CN" sum(w * exports(region = "CN", product = p) for (p, w) in containers) <= 6.0

## two products harvested from the same stands cannot take more than they hold
Limit(name = "SEF stands",
      expression = Markets.supply(region = "SEF", product = ["swr", "smallwood"]),
      upper = 12.4)

## no more than a third of the sawnwood used in OF may come from abroad
@limit "OF imports" imports(region = "OF", product = "sawn_sw") <=
                    consumption(region = "OF", product = "sawn_sw") / 3

## imports held near last year's, unless it is worth 50 a unit to go further
@limit "pulp imports DE" imports(region = "DE", product = "pulp") <= 1.1 * last_year penalty = 50
```
"""
Base.@kwdef struct Limit
    "Name of the limit, used to report it"
    name::String
    """
    What is bounded: a [`LimitExpression`](@ref), without its constant, which
    has been moved into the bounds
    """
    expression::LimitExpression
    "Lower bound, `-Inf` (the default) for none"
    lower::Float64 = -Inf
    "Upper bound, `Inf` (the default) for none"
    upper::Float64 = Inf
    """
    Cost per unit of crossing a bound: `Inf` (the default) for a hard limit,
    which cannot be crossed, or a positive number for a soft one
    """
    penalty::Float64 = Inf

    function Limit(name, expression, lower, upper, penalty)
        e = as_expression(expression)
        terms = [t for t in e.terms if !iszero(last(t))]
        isempty(terms) && throw(ArgumentError(
            "limit $name bounds no quantity of the model"))
        isfinite(lower) || isfinite(upper) || throw(ArgumentError(
            "limit $name has neither a lower nor an upper bound"))
        lower <= upper || throw(ArgumentError(
            "limit $name has a lower bound ($lower) above its upper bound ($upper)"))
        penalty > 0 || throw(ArgumentError(
            "the penalty of limit $name must be positive (Inf for a hard limit), got $penalty"))
        new(nm(name), LimitExpression(terms, 0.0), lower - e.constant, upper - e.constant,
            penalty)
    end
end

# the options a limit takes after its inequality, as `penalty = 50`
function limit_options(options)
    map(options) do o
        o isa Expr && o.head in (:(=), :kw) && length(o.args) == 2 && o.args[1] === :penalty ||
            throw(ArgumentError("after its inequality a limit takes only `penalty = ...`, " *
                                "got $(repr(o))"))
        Expr(:kw, :penalty, esc(o.args[2]))
    end
end

"""
    @limit name lhs <= rhs
    @limit name lhs >= rhs
    @limit name lhs == rhs
    @limit name lower <= expression <= upper
    @limit name lhs <= rhs penalty = p

A [`Limit`](@ref) written as the inequality it stands for. Inside the macro the
term builders — `supply`, `consumption`, `production`, `activity`, `inputs`,
`flow`, `exports` and `imports` — are those of Markets.jl, whatever the names
mean outside it, so they need no qualification; quantities may appear on both
sides, and any constant is moved into the bound. `name` is any expression
giving a string. A trailing `penalty = p` makes the limit soft: it can be
crossed at a cost of `p` per unit (see [`Limit`](@ref)).

# Example
```julia
@limit "SEF stands" supply(region = "SEF", product = "swr") +
                    supply(region = "SEF", product = "smallwood") <= 12.4

## the same, for every subset of the products a region's stands can supply
[@limit "\$r: \$(join(C, '+'))" sum(supply(region = r, product = p) for p in C) <= stock(r, C)
 for r in regions for C in subsets]

@limit "milled at home" activity(region = "SWF", process = "sawmill") >=
                        0.5 * supply(region = "SWF", product = "logs")
@limit "export band" 1.0 <= exports(region = "RW", product = "logs") <= 3.0
```
"""
macro limit(name, ex, options...)
    extra = limit_options(options)
    qualify(x) = !(x isa Expr) ? x :
        x.head == :call && x.args[1] in QUANTITY_BUILDERS ?
            Expr(:call, GlobalRef(@__MODULE__, x.args[1]), map(qualify, x.args[2:end])...) :
            Expr(x.head, map(qualify, x.args)...)
    build(expression; bounds...) = Expr(:call, GlobalRef(@__MODULE__, :Limit),
        Expr(:parameters, Expr(:kw, :name, esc(name)),
             Expr(:kw, :expression, expression),
             (Expr(:kw, k, v) for (k, v) in bounds)..., extra...))
    lower_ops, upper_ops = (:>=, :≥), (:<=, :≤)
    if ex isa Expr && ex.head == :call && length(ex.args) == 3 &&
       ex.args[1] in (lower_ops..., upper_ops..., :(==))
        op, a, b = ex.args
        difference = :($(GlobalRef(@__MODULE__, :as_expression))($(esc(qualify(a))) - $(esc(qualify(b)))))
        op in upper_ops && return build(difference; upper = 0.0)
        op in lower_ops && return build(difference; lower = 0.0)
        return build(difference; lower = 0.0, upper = 0.0)
    elseif ex isa Expr && ex.head == :comparison && length(ex.args) == 5
        lo, op1, e, op2, hi = ex.args
        expression = :($(GlobalRef(@__MODULE__, :as_expression))($(esc(qualify(e)))))
        op1 in upper_ops && op2 in upper_ops &&
            return build(expression; lower = esc(lo), upper = esc(hi))
        op1 in lower_ops && op2 in lower_ops &&
            return build(expression; lower = esc(hi), upper = esc(lo))
    end
    throw(ArgumentError("@limit takes a name and an inequality or an equality, such as " *
                        "`@limit \"cap\" supply(product = \"logs\") <= 10`, got $(repr(ex))"))
end

"""
    @limits name [i in I, j in J, ...] lhs <= rhs
    @limits name [i in I, j in J; condition] lhs <= rhs
    @limits name [i in I, j in J] lhs <= rhs penalty = p

A family of [`Limit`](@ref)s, one for every combination of the indices, as a
constraint family of JuMP or an indexed equation of GAMS: the same inequality,
written as for [`@limit`](@ref), with the indices free in it. Each member is
named after the family and its indices — `cut[SEF]`, `stands[SEF, swr]` — and
the result is a plain vector of limits. An index may be written `i in I`,
`i ∈ I` or `i = I`, and a condition after a semicolon keeps only the
combinations that meet it. A trailing `penalty = p` makes every member soft.

To read a family back from the results, select its members by name:
`res.limits[startswith.(res.limits.limit, "cut["), :]`.

# Example
```julia
@limits "cut" [r in france] supply(region = r, product = "swr") <= cap[r]

## one limit per region and product, leaving out the products a region lacks
@limits "stands" [r in regions, p in primary; haskey(stock, (r, p))] supply(region = r, product = p) <= stock[(r, p)]
```
"""
macro limits(name, indices, ex, options...)
    limit_options(options)          # checked here, so the error names the family
    fail() = throw(ArgumentError("@limits takes a name, indices such as `[r in regions]` " *
                                 "or `[r in regions, p in products; condition]`, and an " *
                                 "inequality, got $(repr(indices))"))
    indices isa Expr || fail()
    condition = nothing
    if indices.head == :vcat && length(indices.args) == 2      # [i in I; condition]
        ranges, condition = Any[indices.args[1]], indices.args[2]
    elseif indices.head == :vect
        ranges = Any[a for a in indices.args if !(a isa Expr && a.head == :parameters)]
        params = [a for a in indices.args if a isa Expr && a.head == :parameters]
        isempty(params) || (condition = only(only(params).args))
    else
        fail()
    end
    loops = map(ranges) do r
        r isa Expr || fail()
        if r.head == :call && length(r.args) == 3 && r.args[1] in (:in, :∈)
            r.args[2], r.args[3]
        elseif r.head == :(=) && length(r.args) == 2
            r.args[1], r.args[2]
        else
            fail()
        end
    end
    vars = first.(loops)
    label = :(string($name, "[", join(string.(($(vars...),)), ", "), "]"))
    member = Expr(:macrocall, GlobalRef(@__MODULE__, Symbol("@limit")), __source__, label, ex,
                  options...)
    # the last index varies fastest, as in nested loops written in that order
    iteration = Expr(:block, (Expr(:(=), v, collection) for (v, collection) in reverse(loops))...)
    generator = condition === nothing ?
        Expr(:generator, member, iteration.args...) :
        Expr(:generator, member, Expr(:filter, condition, iteration.args...))
    # several indices make a product, collected as an array of their shape and
    # handed back flat
    return esc(Expr(:call, GlobalRef(Base, :vec),
                    Expr(:call, GlobalRef(Base, :collect), GlobalRef(@__MODULE__, :Limit),
                         generator)))
end

# --- display ------------------------------------------------------------------

function Base.show(io::IO, q::QuantityRef)
    shown = [string(k, " = ", v isa AbstractString ? v : join(v, "|"))
             for (k, v) in pairs(q.selection) if v != (:)]
    print(io, q.kind, "(", join(shown, ", "), ")")
end

function Base.show(io::IO, e::LimitExpression)
    isempty(e.terms) && return print(io, e.constant)
    for (i, (q, w)) in enumerate(e.terms)
        sign = w < 0 ? (i == 1 ? "−" : " − ") : (i == 1 ? "" : " + ")
        a = abs(w)
        print(io, sign, isone(a) ? "" : string(round(a; sigdigits = 4), " "), q)
    end
    iszero(e.constant) || print(io, e.constant < 0 ? " − " : " + ", abs(e.constant))
end

function Base.show(io::IO, l::Limit)
    print(io, "Limit(", repr(l.name), ": ")
    l.lower == l.upper ? print(io, l.expression, " = ", l.upper) :
        print(io, isfinite(l.lower) ? string(l.lower, " ≤ ") : "", l.expression,
              isfinite(l.upper) ? string(" ≤ ", l.upper) : "")
    isfinite(l.penalty) && print(io, ", penalty ", l.penalty)
    print(io, ")")
end
