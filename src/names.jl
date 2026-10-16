# ----------------------------------------------------------------------------
#  Names
#
#  Every name in the API — a region, a product, a process, a tax — is a string.
#  Symbols are accepted wherever a name is wanted, and converted on the spot,
#  so `:EU` and `"EU"` may be mixed freely; what gets stored, compared and
#  reported back is always the string.
#
#  Symbols would make the dictionary lookups in the model loop cheaper, by
#  roughly 7 ns each.  A model of the size this package is for spends two
#  seconds in the solver, so even a hundred thousand lookups would save under
#  a millisecond of it — and the price would be that `spec.region == "EU"`
#  quietly answers `false` for anyone who reads a field back.  Nothing in a
#  model that warns about silent misreadings should ship that trade.
# ----------------------------------------------------------------------------

"""
What the API accepts wherever it wants a name: an `AbstractString` or a
`Symbol`.  `:EU` and `"EU"` mean the same thing everywhere, and both are
stored, compared and reported as the string `"EU"`.
"""
const Name = Union{Symbol,AbstractString}

"""
$(TYPEDSIGNATURES)

The name `x` as the string the package stores.  Throws an `ArgumentError` on
anything that is not a name, which is where a misplaced number or a forgotten
`=>` is caught.
"""
nm(x::AbstractString) = String(x)
nm(x::Symbol)         = String(x)
nm(x) = throw(ArgumentError(
    "a name must be a string or a symbol, got $(repr(x))::$(typeof(x))"))

"The names `xs` as strings."
nms(xs) = String[nm(x) for x in xs]

"""
$(TYPEDSIGNATURES)

A selector: one name, several names in a vector, or `"all"` — the default of
every selector field — for every one.

`"all"` is a sentinel, so a region or product genuinely called `all` cannot be
singled out by naming it; select it as the one-element vector `["all"]`, or
rename it.
"""
selector(x::Name)            = nm(x)
selector(xs::AbstractVector) = nms(xs)
selector(xs::Tuple)          = nms(collect(xs))
selector(x) = throw(ArgumentError(
    "a selector must be a name, a vector of names, or \"all\", got " *
    "$(repr(x))::$(typeof(x))"))

"Whether a selector field matches a single name."
matches(field, x) = field == "all" || (field isa AbstractString ? field == x : x in field)

"A `name => number` dictionary, keyed by strings."
namekeys(d::AbstractDict) = Dict{String,Float64}(nm(k) => float(v) for (k, v) in d)

"A dictionary keyed by a pair of names, e.g. `(origin, destination)`."
namekeys2(d::AbstractDict) = Dict{Tuple{String,String},Float64}(
    (nm(k[1]), nm(k[2])) => float(v) for (k, v) in d)

"A dictionary keyed by a triple of names, e.g. `(product, from, to)`."
namekeys3(d::AbstractDict) = Dict{Tuple{String,String,String},Float64}(
    (nm(k[1]), nm(k[2]), nm(k[3])) => float(v) for (k, v) in d)

"A vector of `name => number` pairs, keyed by strings."
namepairs(ps) = Pair{String,Float64}[nm(first(p)) => float(last(p)) for p in ps]

"""
A bound that is either one number for every region or a `region => bound`
dictionary.
"""
namebound(x::Real)         = float(x)
namebound(d::AbstractDict) = namekeys(d)
namebound(x) = throw(ArgumentError(
    "a bound must be a number or a `name => number` dictionary, got " *
    "$(repr(x))::$(typeof(x))"))
