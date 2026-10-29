# Part of Markets.jl package [https://github.com/sylvaticus/Markets.jl] - MIT Licence
# File verified 20261029 Antonello LOBIANCO

# ----------------------------------------------------------------------------
#  Data schema  (this is all a user touches to describe an economy)
#
#  All the types and the builders below take keyword arguments only.
# ----------------------------------------------------------------------------

"""
$(TYPEDEF)

A group of products that substitute for each other in final demand.

Demand in this model is integrable by construction: it is the gradient of the
consumer-benefit term of the objective, which is what makes the
surplus-maximisation equivalence hold. An arbitrary matrix of cross-price
elasticities is **not** integrable — it is the same obstruction that rules out
an arbitrary matrix of Armington elasticities — so substitution between
products is expressed the same way, as a tree. Products that substitute closely
go in a nest with a high `sigma`; the nest as a whole then substitutes against
the rest at the lower `sigma` of the level above, and the whole bundle faces one
[`DemandCurve`](@ref) curve.

The implied cross-price effects are the usual CES ones: within a nest the
compensated cross-price elasticity of product `i` with respect to the price of
`j` is ``\\sigma s_j``, with `s_j` the value share of `j`, so a large supplier
transmits its price more strongly than a small one at the same `sigma`.

Substitution against goods **outside** the model — public transport against
cars, say — needs none of this. An exogenous price is a constant, so its effect
is a shift of the whole curve: use `shift` with [`exogenous_shift`](@ref) and
the cross-price elasticity as the driver elasticity.

# Fields
$(TYPEDFIELDS)

# Example
```julia
## petrol and diesel cars are close substitutes; an electric car is a more
## distant one for both
DemandNest(sigma = 1.6, products = ["electric_car"], shares = [0.25],
           nests = [DemandNest(sigma = 5.0, share = 0.75,
                               products = ["petrol_car", "diesel_car"],
                               shares   = [0.65, 0.35])])
```
"""
Base.@kwdef struct DemandNest
    """
    Elasticity of substitution ``\\sigma`` between the members of the group. It
    must be `≥ 0`, and at least the `sigma` of the nest enclosing it. Below 1
    the members are poor substitutes, or complements; `0` makes them be
    consumed in fixed proportions (Leontief), `1` is Cobb–Douglas, and `Inf`
    makes them perfect substitutes — one good to the consumer, bought for
    whichever of them is cheapest
    """
    sigma::Float64
    "Products belonging directly to the group"
    products::Vector{String} = String[]
    """
    Value shares of those products at equal prices, in the same order. They are
    normalised, together with the `share` of each sub-group, over the members
    of this node. Equal shares when omitted
    """
    shares::Vector{Float64} = fill(1.0, length(products))
    "Sub-groups of the group, for a deeper tree"
    nests::Vector{DemandNest} = DemandNest[]
    "Weight of this group among the members of the group enclosing it"
    share::Float64 = 1.0

    function DemandNest(sigma, products, shares, nests, share)
        check_nesting(sigma, nests, "a demand group")
        isempty(products) && isempty(nests) && throw(ArgumentError(
            "a demand group must contain at least one product or sub-group"))
        length(products) == length(shares) || throw(ArgumentError(
            "a demand group needs one share per product, got $(length(products)) " *
            "products and $(length(shares)) shares"))
        all(>(0), shares) && share > 0 || throw(ArgumentError(
            "the shares of a demand group must be positive"))
        allunique([nms(products); leaves(nests)]) || throw(ArgumentError(
            "a product appears more than once in a demand group"))
        new(sigma, nms(products), shares, nests, share)
    end
end

# The two kinds of group a user writes — products that substitute in final
# demand, origins that substitute in a destination — share their shape: an
# elasticity, direct members, sub-groups. What differs is only what the leaves
# are called, so the traversals and the nesting rule are written once.
"The direct members of a group: its products, or its origins."
members(n::DemandNest) = n.products

"""
$(TYPEDSIGNATURES)

Every leaf — product or origin — appearing anywhere in the group `n`, in a
vector of groups, or a bare name, which is its own leaf.
"""
leaves(n, acc = String[]) = (append!(acc, members(n)); foreach(sub -> leaves(sub, acc), n.nests); acc)
leaves(ns::AbstractVector, acc = String[]) = (foreach(n -> leaves(n, acc), ns); acc)
leaves(p::AbstractString, acc = String[]) = push!(acc, String(p))

"""
$(TYPEDSIGNATURES)

Check that every sub-group of a group of elasticity `sigma` substitutes at
least as easily inside as with the outside: the substitution between two
leaves is the elasticity of the smallest group containing both, which only
means something if the elasticities fall from the leaves to the root.
"""
function check_nesting(sigma, nests, what)
    sigma >= 0 || throw(ArgumentError("$what's elasticity must be ≥ 0, got $sigma"))
    for n in nests
        n.sigma >= sigma || throw(ArgumentError(
            "a sub-group of $what substitutes less inside ($(n.sigma)) than with the " *
            "outside ($sigma): members of a group must substitute for each other at " *
            "least as easily as for non-members"))
    end
end

"""
$(TYPEDEF)

A demand or a supply curve for one product in one region: `side` is `:demand`
or `:supply`, and the two aliases [`DemandCurve`](@ref) and
[`SupplyCurve`](@ref) are how a curve is written and read. The curve is given
in one of two ways:

* by a **reference point** `(q0, p0)` and the own-price `elasticity` there.
  The curve has constant elasticity by default, `D = q0 · (p/p0)^(-η)` or
  `S = q0 · (p/p0)^ε`, or is the straight line through the point with that
  slope, with `form = :linear`;
* by the **coefficients** of a straight line, `Q = alpha + beta · p`, with
  `beta < 0` for demand and `beta > 0` for supply. A demand `alpha` is the
  quantity demanded at a zero price and must be positive; a supply `alpha` may
  have either sign: negative, nothing is supplied below the price
  `-alpha/beta`; positive, something is supplied even at a zero price.

Either way, `shift` multiplies the quantity at every price.

# Fields
$(TYPEDFIELDS)
"""
struct Curve{side}
    """
    The product, or, for a demand curve, a [`DemandNest`](@ref) of products
    that substitute for each other, in which case the curve describes the
    **composite** of the group
    """
    product::Union{String,DemandNest}
    "Region where the product is consumed or supplied"
    region::String
    """
    Shape of the curve: `:isoelastic`, the same elasticity at every price, or
    `:linear`, the same slope. Giving `alpha` and `beta` makes it `:linear`
    """
    form::Symbol
    "Reference price, i.e. the price at which the quantity is `q0`; `NaN` for a line given by its coefficients"
    p0::Float64
    "Reference quantity, i.e. the quantity at price `p0`; `NaN` for a line given by its coefficients"
    q0::Float64
    """
    Own-price elasticity at the reference point — ``\\eta`` of demand,
    ``\\varepsilon`` of supply — which must be positive. Any value is allowed:
    a demand with ``\\eta \\le 1`` makes the benefit term of the objective an
    antiderivative of the inverse demand rather than the integral from zero,
    which leaves the equilibrium and the prices unchanged but makes the
    reported objective value a welfare *difference* rather than a level.
    `Inf` makes the curve flat, in either form: any quantity is bought or sold
    at `p0`, as an export market or imports at a given world price are (`q0`,
    which then plays no role, must still be positive, and `shift` has no
    effect). A flat supply that can reach a flat demand paying more than it
    costs, transport and taxes included, leaves the market without a finite
    answer. `NaN` for a line given by its coefficients
    """
    elasticity::Float64
    """
    Quantity at a zero price, ``\\alpha`` in ``Q = \\alpha + \\beta p``, for a
    linear curve; given, or derived from the reference point. `NaN` for a
    constant-elasticity curve
    """
    alpha::Float64
    """
    Slope ``\\beta = dQ/dp`` of a linear curve, negative for demand and
    positive for supply; given, or derived from the reference point. `NaN` for
    a constant-elasticity curve
    """
    beta::Float64
    """
    Multiplicative shifter of the whole curve: the quantity at any given price
    is `shift` times what the curve otherwise implies. On a line, it scales
    `alpha` and `beta` together, so the price at which the quantity is zero
    does not move; a flat curve has no quantity to scale and refuses any
    shifter but 1. It is how everything other than the price enters — income,
    population and the price of a substitute outside the model for demand; the
    size of the resource, the roads and ports, the labour and machinery
    available for supply — leaving the reference point or the coefficients to
    hold the calibration. [`exogenous_shift`](@ref) builds one from drivers and
    their elasticities. The default, 1, is the calibrated curve itself. A hard
    ceiling on a quantity — an allowable cut, a quota, an exhausted resource —
    is a [`Limit`](@ref) instead
    """
    shift::Float64

    function Curve{side}(product, region, form, p0, q0, elasticity, alpha, beta,
                         shift) where side
        side in (:demand, :supply) || throw(ArgumentError(
            "a curve is a Curve{:demand} or a Curve{:supply}, got Curve{$(repr(side))}"))
        side === :supply && product isa DemandNest && throw(ArgumentError(
            "a supply curve describes one product, not a group of substitutes"))
        what = "$(product isa DemandNest ? "a group" : product) in $region"
        shift > 0 || throw(ArgumentError(
            "the $side shifter of $what must be positive, got $shift"))
        isinf(elasticity) && shift != 1 && throw(ArgumentError(
            "the $side of $what is flat (elasticity Inf): it takes any quantity at " *
            "p0, so a shifter has nothing to scale. Change p0 to move it"))
        new{side}(product isa DemandNest ? product : nm(product), nm(region),
                  form, p0, q0, elasticity, alpha, beta, shift)
    end
end

"""
    DemandCurve(; product, region, p0, q0, elasticity, form = :isoelastic, shift = 1)
    DemandCurve(; product, region, alpha, beta, shift = 1)
    DemandCurve(c::DemandCurve; shift)

*Final* demand for one product — or for a [`DemandNest`](@ref) of substitutes
— in one region: a [`Curve`](@ref) with `side = :demand`, given by a reference
point and its elasticity or by the coefficients of a straight line. The third
form is the curve `c` with the shifter `shift` in place of its own: the same
calibration, whatever its form, moved by a new multiplier, which is how a
scenario moves a curve.

# Example
```julia
DemandCurve(product = "gasoline", region = "EU", p0 = 900, q0 = 80, elasticity = 0.3)
DemandCurve(product = "gasoline", region = "EU", p0 = 900, q0 = 80, elasticity = 0.3,
            form = :linear)                       # the tangent line at the same point
DemandCurve(product = "wheat", region = "A", alpha = 200, beta = -9)
thirsty = DemandCurve(milk; shift = 1.2)          # a fifth more milk at every price
```
"""
const DemandCurve = Curve{:demand}

"""
    SupplyCurve(; product, region, p0, q0, elasticity, form = :isoelastic, shift = 1)
    SupplyCurve(; product, region, alpha, beta, shift = 1)
    SupplyCurve(c::SupplyCurve; shift)

*Primary* supply (e.g. crude oil from wells, ore from a mine, fish from the
sea, crops from land) of one product in one region: a [`Curve`](@ref) with
`side = :supply`, given by a reference point and its elasticity or by the
coefficients of a straight line. The third form is the curve `c` with the
shifter `shift` in place of its own: the same calibration, whatever its form,
moved by a new multiplier, which is how a scenario moves a curve.

# Example
```julia
SupplyCurve(product = "crude", region = "US", p0 = 550, q0 = 600, elasticity = 0.4)
SupplyCurve(product = "crude", region = "US", p0 = 550, q0 = 600, elasticity = 0.4,
            form = :linear)                       # the tangent line at the same point
SupplyCurve(product = "wheat", region = "A", alpha = -50, beta = 10.5)
salvage = SupplyCurve(logs; shift = 3.0)          # three times the logs at every price
```
"""
const SupplyCurve = Curve{:supply}

Curve{side}(; product, region, p0 = nothing, q0 = nothing, elasticity = nothing,
            alpha = nothing, beta = nothing, form = nothing, shift = 1.0) where side =
    Curve{side}(product, region,
                curve_parameters(side, "$(product isa DemandNest ? "a group" : product) in $region",
                                 p0, q0, elasticity, alpha, beta, form)...,
                shift)

Curve{side}(c::Curve{side}; shift::Real) where side =
    Curve{side}(c.product, c.region, c.form, c.p0, c.q0, c.elasticity, c.alpha, c.beta,
                float(shift))

# Julia picks a method by the types of the positional arguments, never by the
# names of the keywords, so the two ways of giving a curve share one keyword
# constructor, and this sorts out which one was used: it returns the fields
# `(form, p0, q0, elasticity, alpha, beta)`.
function curve_parameters(side::Symbol, what::String, p0, q0, elasticity, alpha, beta, form)
    given(x) = x !== nothing
    point, line = (p0, q0, elasticity), (alpha, beta)
    shape = form === nothing ? nothing : Symbol(form)
    shape in (nothing, :isoelastic, :linear) || throw(ArgumentError(
        "the $side of $what has form $(repr(form)): it must be :isoelastic or :linear"))
    if any(given, line)
        any(given, point) && throw(ArgumentError(
            "the $side of $what is given both by a reference point (p0, q0, elasticity) " *
            "and by the coefficients of a line (alpha, beta): give one or the other"))
        all(given, line) || throw(ArgumentError(
            "the linear $side of $what needs both alpha and beta"))
        shape === :isoelastic && throw(ArgumentError(
            "alpha and beta describe a straight line, so the $side of $what cannot " *
            "be :isoelastic"))
        α, β = float(alpha), float(beta)
        isfinite(α) && isfinite(β) || throw(ArgumentError(
            "the coefficients of the $side of $what must be finite, got alpha = $α, beta = $β"))
        if side === :demand
            α > 0 && β < 0 || throw(ArgumentError(
                "a linear demand D = alpha + beta·p needs alpha > 0 and beta < 0; " *
                "the demand for $what has alpha = $α, beta = $β"))
        else
            β > 0 || throw(ArgumentError(
                "a linear supply S = alpha + beta·p needs beta > 0; " *
                "the supply of $what has beta = $β"))
        end
        return (:linear, NaN, NaN, NaN, α, β)
    end
    all(given, point) || throw(ArgumentError(
        "the $side of $what needs either a reference point and its elasticity " *
        "(p0, q0, elasticity) or the coefficients of a line (alpha, beta)"))
    p, q, e = float(p0), float(q0), float(elasticity)
    p > 0 && q > 0 || throw(ArgumentError(
        "the reference point of the $side of $what must be positive, got p0 = $p, q0 = $q"))
    e > 0 || throw(ArgumentError(
        "the $side elasticity of $what must be positive, got $e"))
    # a flat curve is the same in both forms, and has no finite slope
    (shape === :linear && isfinite(e)) || return (:isoelastic, p, q, e, NaN, NaN)
    β = (side === :demand ? -e : e) * q / p
    return (:linear, p, q, e, q - β * p, β)
end

# The sign with which the elasticity enters the exponent: the inverse demand
# falls with the quantity, the inverse supply rises.
exponent_sign(::Curve{:demand}) = -1
exponent_sign(::Curve{:supply}) = +1

# Where a curve sits, as a quantity: its reference quantity, or, for a line
# given by its coefficients, a quantity of its order of magnitude (half what is
# demanded at a zero price; for supply what is supplied at a zero price, or
# offered at twice the price where supply starts). It starts the solver and
# scales the numerical tolerances; it does not calibrate anything.
reference_quantity(c::DemandCurve) = isnan(c.q0) ? c.alpha / 2 : c.q0
reference_quantity(c::SupplyCurve) = isnan(c.q0) ? (iszero(c.alpha) ? c.beta : abs(c.alpha)) : c.q0

"""
$(TYPEDSIGNATURES)

The price along the curve `c` at quantity `q`, before any tax: the willingness
to pay for the `q`-th unit on a demand curve, the marginal cost of it on a
supply curve.
"""
curve_price(c::Curve, q::Real) =
    c.form === :linear ? (q / c.shift - c.alpha) / c.beta :
                         c.p0 * (q / (c.shift * c.q0))^(exponent_sign(c) / c.elasticity)

"""
$(TYPEDSIGNATURES)

The integral of the inverse curve `c` up to the quantity `q`, as a function of
`q` that may be a model variable: the gross consumer benefit of a demand curve,
the resource cost of a supply curve. For a constant-elasticity curve it is an
antiderivative, not the integral from zero, which diverges for a demand with
``\\eta \\le 1``: the constant changes neither the optimum nor the prices.
"""
function curve_integral(c::Curve, q)
    c.form === :linear && return (q^2 / (2 * c.shift) - c.alpha * q) / c.beta
    e = 1 + exponent_sign(c) / c.elasticity
    a = c.p0 * (c.shift * c.q0)^(-exponent_sign(c) / c.elasticity)
    # at η = 1 the exponent is zero and the power form is undefined: the
    # antiderivative is then a·log(q)
    return abs(e) < 1e-8 ? a * log(q) : (a / e) * q^e
end

function Base.show(io::IO, c::Curve{side}) where side
    print(io, side === :demand ? "DemandCurve" : "SupplyCurve", "(product = ", repr(c.product),
          ", region = ", repr(c.region))
    if isnan(c.p0)
        print(io, ", alpha = ", c.alpha, ", beta = ", c.beta)
    else
        print(io, ", p0 = ", c.p0, ", q0 = ", c.q0, ", elasticity = ", c.elasticity)
        c.form === :linear && print(io, ", form = :linear")
    end
    c.shift == 1 || print(io, ", shift = ", c.shift)
    print(io, ")")
end

"""
$(TYPEDSIGNATURES)

A log-linear shifter of a supply or demand curve, built from the exogenous
drivers that move it and the elasticity of the curve with respect to each:

```math
\\text{shift} \\;=\\; \\prod_j \\left(\\frac{Z_j}{Z_j^0}\\right)^{\\gamma_j}
```

Each driver is given as a `NamedTuple` with its current `level`, its
`reference` level — the one the curve was calibrated at, which gives a shift of
1 — and the `elasticity` of the curve with respect to it. Pass them as keyword
arguments, named for readability only.

This is the form that comes out of the usual estimation. Regressing
``\\log S`` on ``\\log P`` and on the ``\\log Z_j`` gives the price
elasticity and the ``\\gamma_j`` in one go, so the shifters cost nothing extra to obtain
once the curve itself has been estimated.

The shifter is the same on a straight line: it multiplies the quantity at
every price, so ``\\gamma_j`` is the elasticity of the quantity to driver
``j`` along the whole line, as along a constant-elasticity curve. A line
estimated in levels, ``Q = \\alpha + \\beta P + \\gamma Z``, has its
drivers in the intercept instead: there the driver moves `alpha`, by
``\\gamma (Z - Z^0)``, and is given as a new `alpha` rather than as a shifter.

# Example
```julia
# a cod stock 18% larger than at calibration, and a slightly larger fleet
shift = exogenous_shift(
    stock = (level = 118.0, reference = 100.0, elasticity = 0.6),
    fleet = (level = 1.05,  reference = 1.0,   elasticity = 0.25))

SupplyCurve(product = "cod", region = "EU", p0 = 3000, q0 = 300,
           elasticity = 0.40, shift = shift)
```
"""
function exogenous_shift(; drivers...) # verified 20261016 Antonello
    isempty(drivers) && throw(ArgumentError("give at least one driver"))
    shift = 1.0
    for (name, d) in pairs(drivers)
        haskey(d, :level) && haskey(d, :reference) && haskey(d, :elasticity) ||
            throw(ArgumentError(
                "driver $name needs `level`, `reference` and `elasticity`, got $(keys(d))"))
        d.level > 0 && d.reference > 0 || throw(ArgumentError(
            "the level and reference of driver $name must be positive, " *
            "got $(d.level) and $(d.reference)"))
        shift *= (d.level / d.reference)^d.elasticity
    end
    return shift
end

"""
$(TYPEDEF)

One *input requirement* of a [`Process`](@ref). Inputs in the same nest
substitute for each other; different nests of the same process are needed in
fixed proportion to the process activity.

* A **Leontief** (fixed-coefficient) requirement is a nest with a single
  product: `composite` units of it are needed per unit of process activity.
  Build it with [`leontief`](@ref).
* A **CES** (smoothly substitutable) requirement bundles several `products`
  into a composite via a CES aggregator with value `shares` (δ, Σδ = 1),
  elasticity of substitution `sigma` (σ) and scale `phi` (φ):

      composite · z  ≤  φ · ( Σᵢ δᵢ^(1/σ) · qᵢ^ρ )^(1/ρ),     ρ = (σ-1)/σ

  As an input's price rises the optimiser smoothly shifts the mix toward
  cheaper inputs. Build it with [`ces`](@ref). `sigma` may be any value
  `≥ 0`: the formula is replaced by its limit at `σ = 0` (fixed proportions,
  `minᵢ qᵢ/δᵢ`), at `σ = 1` (Cobb–Douglas, `Πᵢ (qᵢ/δᵢ)^δᵢ`) and at `σ = ∞`
  (perfect substitutes, `Σᵢ qᵢ`).

Use the [`leontief`](@ref) and [`ces`](@ref) builders rather than this
constructor: they normalise the shares and check the arguments.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct Nest
    "Units of the (composite) input needed per unit of process activity"
    composite::Float64
    "Products of the nest. A single product makes the nest a Leontief one"
    products::Vector{String}
    """
    Value shares ``\\delta`` of each product, summing to 1 (equal shares by
    default). They set the input mix when all input prices are equal. Ignored
    for a single-product (Leontief) nest
    """
    shares::Vector{Float64} = fill(1 / length(products), length(products))
    """
    Elasticity of substitution ``\\sigma`` between the products of the nest.
    It must be `≥ 0`: `0` for inputs needed in fixed proportions, `1` for
    Cobb–Douglas, `Inf` for perfect substitutes. It is ignored (and
    conventionally set to 0) for a single-product (Leontief) nest
    """
    sigma::Float64
    "Scale parameter ``\\phi`` of the CES aggregator"
    phi::Float64 = 1.0

    function Nest(composite, products, shares, sigma, phi)
        composite > 0 || throw(ArgumentError(
            "a nest must need a positive quantity of its input, got $composite"))
        isempty(products) && throw(ArgumentError("a nest must have at least one product"))
        length(products) == length(shares) || throw(ArgumentError(
            "a nest needs one share per product, got $(length(products)) products " *
            "and $(length(shares)) shares"))
        allunique(products) || throw(ArgumentError(
            "the products of a nest must differ, got $products"))
        all(>(0), shares) || throw(ArgumentError("the shares of a nest must be positive"))
        phi > 0 || throw(ArgumentError("the scale of a nest must be positive, got $phi"))
        # a single-product nest is a fixed coefficient and never consults sigma
        length(products) == 1 || sigma >= 0 || throw(ArgumentError(
            "a nest of several products needs an elasticity of substitution ≥ 0, " *
            "got $sigma"))
        # ... nor phi, which it would otherwise swallow in silence
        length(products) == 1 && phi != 1 && throw(ArgumentError(
            "a nest of the single product $(only(products)) is a fixed coefficient, " *
            "and its scale phi = $phi would be ignored: fold it into the coefficient"))
        new(composite, nms(products), shares, sigma, phi)
    end
end

"""
$(TYPEDSIGNATURES)

Build a fixed-coefficient (Leontief) input requirement: `coeff` units of
`product` per unit of process activity.

# Example
```julia
leontief(product = "flour", coeff = 0.75)    # 0.75 t of flour per t of bread
```
"""
leontief(; product::Name, coeff::Real) =
    Nest(composite = coeff, products = [product], shares = [1.0], sigma = 0.0, phi = 1.0)

"""
$(TYPEDSIGNATURES)

Build a smoothly substitutable (CES) input requirement of `composite` units of
a bundle of `products` per unit of process activity.

`shares` are the value shares of each input when all input prices are equal
(equal shares by default); they need not sum to 1, as they are normalised
internally. `sigma` is the elasticity of substitution and must be `≥ 0`:
below 1 the inputs are poor substitutes, `0` needs them in fixed proportions,
`1` is Cobb–Douglas and `Inf` makes them perfect substitutes.

# Example
```julia
## 1.1 t of a mix of scrap and pig iron per tonne of steel
ces(composite = 1.10, products = ["scrap", "pig_iron"],
    shares = [0.4, 0.6], sigma = 2.0)
```
"""
function ces(; composite::Real, products::AbstractVector,
               shares::Vector{<:Real} = fill(1.0, length(products)),
               sigma::Real, phi::Real = 1.0)
    length(products) == length(shares) || throw(ArgumentError(
        "ces needs one share per product, got $(length(products)) products and " *
        "$(length(shares)) shares"))
    sigma >= 0 || throw(ArgumentError(
        "ces needs an elasticity of substitution ≥ 0, got $sigma"))
    s = collect(float.(shares)); s ./= sum(s)
    Nest(composite = composite, products = products, shares = s, sigma = sigma, phi = phi)
end

"""
$(TYPEDEF)

A transformation activity, operated at a level of activity chosen by the model.

# Fields
$(TYPEDFIELDS)

# Example
```julia
## a refinery splitting each tonne of crude into three fuels
Process(name    = "refinery",
        inputs  = [leontief(product = "crude", coeff = 1.0)],
        outputs = ["gasoline" => 0.45, "diesel" => 0.35, "fuel_oil" => 0.15],
        vacost  = 40)
```
"""
Base.@kwdef struct Process
    "Name of the process"
    name::String
    """
    Input requirements, as a vector of [`Nest`](@ref)s built with
    [`leontief`](@ref) or [`ces`](@ref). Each nest is needed in fixed
    proportion to the process activity
    """
    inputs::Vector{Nest}
    """
    Outputs, as `product => yield` pairs per unit of activity. Several pairs
    make a joint-production process (e.g. petrol *and* diesel from one refinery)
    """
    outputs::Vector{Pair{String,Float64}}
    "Regions where the process may operate, or `:` (the default) for every region"
    regions::Selector = (:)
    """
    Value-added cost per unit of activity: the cost of everything not modelled
    as an explicit input (labour, energy, capital, ...). A single number
    applies to every region the process operates in; a `region => cost`
    dictionary sets them one by one, and must then cover every such region.
    Wages, energy and capital do not cost the same everywhere, and a cost
    calibrated on base-year margins differs from region to region. The
    installed capacity of the plants, which bounds the activity, is a
    [`Limit`](@ref) on it
    """
    vacost::Union{Float64,Dict{String,Float64}} = 0.0

    function Process(name, inputs, outputs, regions, vacost)
        new(nm(name), inputs, namepairs(outputs), selector(regions), namebound(vacost))
    end
end

"Value-added cost per unit of activity of process `p` in region `r`."
vacost_of(p::Process, r::String) = p.vacost isa Dict ? p.vacost[r] : p.vacost

"""
$(TYPEDEF)

A group of origins that substitute for each other more closely than they do for
the rest, inside an [`Armington`](@ref) composite.

Regions whose products are technically interchangeable — because they share
grades, standards or certification, or because the same non-tariff measures
let them through — belong in one group. Within it buyers substitute with the
group's own `sigma`; against anything outside it they substitute with the
`sigma` of the enclosing nest, which is lower.

Groups may be nested further through their own `nests` field, so the structure
is a tree whose leaves are origins. This is the general shape an elasticity
specification may take: the substitution between any two origins is the
elasticity of the smallest group containing both, and a group must be at least
as substitutable inside as it is with the outside.

# Fields
$(TYPEDFIELDS)

# Example
```julia
# in this market, EU and US products are near-interchangeable
OriginNest(sigma = 12, origins = ["EU", "US"])
```
"""
Base.@kwdef struct OriginNest
    """
    Elasticity of substitution ``\\sigma`` between the members of the group. It
    must be `≥ 0`, and at least the `sigma` of the nest that encloses it.
    `Inf` makes the members perfect substitutes: the group is pooled into a
    single good, as in the homogeneous case, while still being a distinct
    variety to the outside
    """
    sigma::Float64
    "Origins belonging to the group"
    origins::Vector{String} = String[]
    "Sub-groups of the group, for a deeper tree"
    nests::Vector{OriginNest} = OriginNest[]

    function OriginNest(sigma, origins, nests)
        check_nesting(sigma, nests, "an origin group")
        isempty(origins) && isempty(nests) && throw(ArgumentError(
            "an origin group must contain at least one origin or sub-group"))
        new(sigma, nms(origins), nests)
    end
end
members(n::OriginNest) = n.origins

"""
$(TYPEDSIGNATURES)

A process that consumes a product and yields nothing: the way to give a
by-product somewhere to go.

Material balances here are equalities, so every unit produced must be used,
traded or consumed. A by-product nobody wants is then a *bad*, and its
equilibrium price is negative — the amount somebody must be paid to take it.
That is the right answer to the question as posed, but it is rarely the
question meant, and under a large supply shock it can drag a local market
somewhere surprising.

Adding a disposal route restores the usual assumption. With `cost = 0` the
price of the product settles at zero rather than going negative, which is free
disposal; with a positive `cost` it settles at minus that cost, which is
disposal that has to be paid for — a landfill charge, an incinerator, a
tailings dam.

# Example
```julia
processes = [refinery, steel_mill,
             disposal(product = "fuel_oil"),           # free disposal
             disposal(product = "slag", cost = 12)]    # paid disposal
```
"""
disposal(; product::Name, cost::Real = 0.0, regions = (:),
           name::Name = string("dispose_", product)) =
    Process(name = name, inputs = [leontief(product = product, coeff = 1.0)],
            outputs = Pair{String,Float64}[], regions = regions, vacost = cost)

"""
$(TYPEDEF)

Imperfect substitution between the regional varieties of one product
(Armington, 1969).

By default a product is **homogeneous**: it is the same good wherever it comes
from, regional prices differ by at most the transport cost, and two regions
never ship it to each other both ways (a region on the cheapest path between
two others may carry it through). Listing a product here instead makes the varieties
from the different origins imperfect substitutes: each region uses a CES
composite of them, so a region can import and export the same product at once
and its price responds to supply and demand everywhere, not only to the cost of
the cheapest source.

Each destination `r` draws on the varieties of the origins it can buy from: its
own, plus every origin `o` with a transport route `(product, o, r)`. The
composite is

    A ≤ ( Σₒ δₒ^(1/σ) · Xₒ^ρ )^(1/ρ),     ρ = (σ-1)/σ

with `Xₒ` the quantity bought from origin `o` and `δₒ` its value share. The
price a region pays for the composite is then the CES price index of the
delivered prices,

    Pᶜ = ( Σₒ δₒ · (Pₒ + τₒ)^(1-σ) )^(1/(1-σ))

A `scale` φ multiplies the aggregate, and divides this price, by the same
factor.

which tends to the cheapest delivered price as `σ → ∞`: `sigma = Inf` is
exactly the homogeneous (Samuelson spatial price equilibrium) case, and a large
finite `sigma` approaches it. At the other end, `sigma = 1` is the
Cobb–Douglas composite `Πₒ (Xₒ/δₒ)^δₒ`, priced at `Πₒ (Pₒ + τₒ)^δₒ`, and
`sigma = 0` the Leontief one `minₒ Xₒ/δₒ`, priced at `Σₒ δₒ (Pₒ + τₒ)`.

A single `sigma` makes every origin substitute equally well for every other,
which is the defining restriction of a flat CES. Two ways out, both available
here:

* **Groups of origins** ([`OriginNest`](@ref)): origins whose products are
  mutually interchangeable are nested together with a higher `sigma`, and
  substitute with the rest at the lower `sigma` of this specification.
* **One specification per destination** (`destination`): the structure may
  differ from one importing market to the next, which is what technical
  standards and non-tariff measures actually do, since the importer sets them.

# Fields
$(TYPEDFIELDS)

# Example
```julia
# moderate substitution, equal shares, same in every market
Armington(product = "wine", sigma = 4)

# calibrated shares: 70% domestic, the rest split between the two other regions
Armington(product = "wine", sigma = 4,
          shares = merge(Dict((r, r) => 0.7 for r in regions),
                         Dict((o, r) => 0.15 for o in regions, r in regions if o != r)))

# in the EU market, EU and US cars are near-interchangeable (shared safety
# standards) while CN cars are a distant substitute; elsewhere everything
# substitutes freely
Armington(product = "cars", destination = "EU", sigma = 2.5,
          nests = [OriginNest(sigma = 12, origins = ["EU", "US"])])
Armington(product = "cars", destination = ["US", "CN"], sigma = 8)
```
"""
Base.@kwdef struct Armington
    "Product whose regional varieties are imperfect substitutes"
    product::String
    """
    Elasticity of substitution ``\\sigma`` between the varieties of the
    different origins. It must be `≥ 0`. Below 1 the varieties are poor
    substitutes, or complements; `0` makes each destination buy them in fixed
    proportions (Leontief) and `1` is Cobb–Douglas. `Inf` (the default for
    every product not listed) means perfect substitutes, i.e. a homogeneous
    product
    """
    sigma::Float64 = Inf
    """
    Value shares ``\\delta`` of each origin in each destination's composite,
    keyed `(origin, destination)` — the shares that would be observed if the
    delivered prices of all origins were equal. They are normalised per
    destination. An origin with no share given (or a share of 0) is left out of
    that destination's composite. When this is empty, every origin available to
    a destination gets an equal share
    """
    shares::Dict{Tuple{String,String},Float64} = Dict{Tuple{String,String},Float64}()
    """
    Importing regions this specification applies to, or `:` (the default)
    for every region without one of their own. A product may have one
    specification per destination plus a `:` fallback, which is how the
    structure, the elasticity and the shares can differ between importing
    markets
    """
    destination::Selector = (:)
    """
    Groups of origins that substitute for each other more closely than for the
    rest, as a vector of [`OriginNest`](@ref)s. Every available origin not
    placed in a group stays a direct member of the composite. Empty (the
    default) gives the flat CES in which all origins substitute equally
    """
    nests::Vector{OriginNest} = OriginNest[]
    """
    Scale ``\\phi`` of the composite, one number for every destination or a
    `destination => scale` dictionary (a destination absent from it gets 1).
    With shares that sum to 1 the CES aggregate equals the plain sum of the
    varieties bought only where their delivered prices are equal, and falls
    short of it where they differ. A scale of the base-year sum over the
    base-year aggregate makes the composite count physical units again at the
    base year, at a price equal to the average delivered unit value, which is
    what lets base-year data in m³ or tonnes balance exactly;
    [`armington_calibration`](@ref) sets it. The default, 1, leaves the
    aggregate as it is
    """
    scale::Union{Float64,Dict{String,Float64}} = 1.0

    function Armington(product, sigma, shares, destination, nests, scale)
        sigma >= 0 || throw(ArgumentError(
            "the Armington elasticity of $product must be ≥ 0 (or Inf for a " *
            "homogeneous product), got $sigma"))
        all(≥(0), values(shares)) || throw(ArgumentError(
            "the Armington shares of $product must not be negative"))
        check_nesting(sigma, nests, "the origins of $product")
        allunique(leaves(nests)) || throw(ArgumentError(
            "an origin appears in more than one group of $product"))
        all(>(0), scale isa AbstractDict ? values(scale) : (scale,)) || throw(ArgumentError(
            "the Armington scale of $product must be positive, got $scale"))
        new(nm(product), sigma, namekeys2(shares), selector(destination), nests,
            namebound(scale))
    end
end

"Scale of the composite of Armington specification `a` in destination `r`."
scale_of(a::Armington, r::String) = a.scale isa Dict ? get(a.scale, r, 1.0) : a.scale

"""
$(TYPEDSIGNATURES)

Calibrate a CES group from the base-year quantities and prices of its members,
so that at those prices the model buys exactly those quantities. The group may
be any CES of the package: substitutes in final demand ([`DemandNest`](@ref)),
inputs of a process ([`ces`](@ref)) or the origins of a product
([`Armington`](@ref), for which [`armington_shares`](@ref) does this for a
whole trade matrix at once).

`quantities` and `prices` give the base-year quantity and price of each member,
in the same order, and `sigma` is the elasticity of substitution the group is
calibrated for. The result is a named tuple of three things the model needs and
nobody observes:

* `shares`: the value shares ``\\delta`` of the members, i.e. the shares they
  would have **if all their prices were equal**. Inverting the CES demand
  condition ``q_i = A\\,\\delta_i (P/p_i)^\\sigma`` gives
  ``\\delta_i \\propto q_i\\, p_i^{\\sigma}``, normalised to sum to 1. They are
  not the observed shares: a member bought little *because* it is dear gets a
  large share, since it would be bought much more at an equal price.
* `price`: the CES price index ``P`` of the group, the price of one unit of its
  composite at the base-year prices.
* `quantity`: the quantity ``A`` of the composite, the base-year expenditure
  ``\\sum_i p_i q_i`` divided by `price`. It is an index number, not a physical
  quantity: m³ of firewood and tonnes of pellets do not add up, while their
  values at a common price do.
* `scale`: ``\\sum_i q_i`` divided by `quantity`. Where the members are
  measured in one unit — m³ of wood of different kinds in a mill's input
  bundle — this is the scale ``\\phi`` that makes the composite count that
  unit again at the base year, at the average unit price
  ``\\sum_i p_i q_i / \\sum_i q_i``: pass it as the `phi` of [`ces`](@ref),
  with the base-year ``\\sum_i q_i`` per unit of activity as `composite`.
  Where the members do not add up it means nothing, and is ignored.

`price` and `quantity` are the reference point of whatever the group feeds:
the `p0` and `q0` of the [`DemandCurve`](@ref) over a [`DemandNest`](@ref), or,
for a group nested inside another, the price and quantity of the group as one
member of the level above.

Where the CES takes its limiting form, so does the calibration:

| `sigma` | `shares` | `price` | `quantity` |
|:--------|:---------|:--------|:-----------|
| `0`   | quantity shares ``q_i / \\sum q`` | ``\\sum_i \\delta_i p_i`` | ``\\sum_i q_i`` |
| `1`   | value shares ``p_i q_i / \\sum p q`` | ``\\prod_i p_i^{\\delta_i}`` | expenditure / `price` |
| `Inf` | quantity shares, which do not enter the model | the common price | ``\\sum_i q_i`` |

Perfect substitutes are bought only at the lowest price, so with `sigma = Inf`
every member with a positive quantity must have that price, to within the
relative tolerance `rtol`. A member with a zero quantity gets a zero share,
which [`DemandNest`](@ref) and [`ces`](@ref) do not accept: leave it out of the
group.

# Example
```julia
## the base year in one region: 7 units of firewood at 55, 1.4 of pellets at 250
fuel = ces_calibration(quantities = [7.0, 1.4], prices = [55.0, 250.0], sigma = 3)
## (shares = [0.0505, 0.9495], price = 177.05, quantity = 4.151, scale = 2.023)
DemandCurve(region = "OF", p0 = fuel.price, q0 = fuel.quantity, elasticity = 0.35,
           product = DemandNest(sigma = 3, products = ["firewood", "pellets"],
                                shares = fuel.shares))
```
"""
function ces_calibration(; quantities::AbstractVector{<:Real},
                           prices::AbstractVector{<:Real}, sigma::Real,
                           rtol::Real = 1e-8)
    q, p = collect(float.(quantities)), collect(float.(prices))
    length(q) == length(p) || throw(ArgumentError(
        "ces_calibration needs one price per quantity, got $(length(q)) quantities " *
        "and $(length(p)) prices"))
    all(>=(0), q) || throw(ArgumentError("quantities must not be negative, got $q"))
    any(>(0), q) || throw(ArgumentError("no positive quantity to calibrate on"))
    all(>(0), p) || throw(ArgumentError("prices must be positive, got $p"))
    sigma >= 0 || throw(ArgumentError("sigma must be ≥ 0, got $sigma"))
    if isinf(sigma)
        # perfect substitutes: the composite is the plain sum, bought at the
        # lowest price, and the shares never enter the model
        price = minimum(p)
        all(isapprox(p[i], price; rtol) for i in eachindex(q) if q[i] > 0) ||
            throw(ArgumentError(
                "with sigma = Inf the members are perfect substitutes, bought only at " *
                "the lowest price, but members bought in the base year have the " *
                "prices $(p[q .> 0])"))
        return (shares = q ./ sum(q), price = price, quantity = sum(q), scale = 1.0)
    end
    weights = q .* p .^ sigma
    shares  = weights ./ sum(weights)
    # the CES price index, and its Cobb–Douglas limit where 1/(1-σ) is undefined
    price = isone(sigma) ? prod(p .^ shares) :
                           sum(shares .* p .^ (1 - sigma))^(1 / (1 - sigma))
    quantity = sum(q .* p) / price
    return (shares = shares, price = price, quantity = quantity, scale = sum(q) / quantity)
end

"""
$(TYPEDSIGNATURES)

Calibrate the [`Armington`](@ref) shares of one product from a base-year trade
matrix, so that the model reproduces the observed sourcing at the observed
prices.

`flows` maps `(origin, destination)` to the quantity the destination bought
from that origin in the base year, its own output included as `(r, r)`.
`prices` maps the same pairs to the **delivered** price of that variety there,
i.e. the producer price of the origin plus the cost of carrying it to the
destination. `sigma` is the elasticity the shares are being calibrated for,
any finite value `≥ 0`.

The shares are *not* the observed quantity or value shares. Inverting the CES
demand condition ``X_o = A\\,\\delta_o (P^c/q_o)^\\sigma`` gives

```math
\\delta_{o,r} \\;\\propto\\; X_{o,r}\\; q_{o,r}^{\\,\\sigma}
```

which is what this returns, normalised per destination: it is
[`ces_calibration`](@ref) applied to the origins of each destination. Feeding observed shares
in directly would only be right if all delivered prices were equal, since
``\\delta`` is by definition the share at equal delivered prices. From observed
*value* shares `s` rather than quantities the equivalent form is
``\\delta \\propto s\\, q^{\\sigma-1}``.

Pairs with a non-positive flow are left out, and an origin left out of a
destination's shares is left out of its composite: a route unused in the base
year stays unused, which is the well-known small-shares limitation of Armington
calibration.

# Example
```julia
shares = armington_shares(
    flows  = Dict(("EU", "EU") => 80.0, ("US", "EU") => 15.0, ("CN", "EU") => 5.0),
    prices = Dict(("EU", "EU") => 210.0, ("US", "EU") => 235.0, ("CN", "EU") => 250.0),
    sigma  = 4)
Armington(product = "wine", sigma = 4, shares = shares)
```
"""
function armington_shares(; flows::AbstractDict, prices::AbstractDict,
                            sigma::Real)
    0 <= sigma < Inf || throw(ArgumentError("sigma must be finite and ≥ 0, got $sigma"))
    flows, prices = namekeys2(flows), namekeys2(prices)
    used = [key for (key, x) in flows if x > 0]
    isempty(used) && throw(ArgumentError("no positive flow to calibrate on"))
    for key in used
        haskey(prices, key) || throw(ArgumentError(
            "no delivered price given for the flow $(key[1]) → $(key[2])"))
        prices[key] > 0 || throw(ArgumentError(
            "the delivered price of $(key[1]) → $(key[2]) must be positive, " *
            "got $(prices[key])"))
    end
    δ = Dict{Tuple{String,String},Float64}()
    for dest in unique(k[2] for k in used)        # one CES group per destination
        keys_d = [k for k in used if k[2] == dest]
        c = ces_calibration(quantities = [flows[k] for k in keys_d],
                            prices = [prices[k] for k in keys_d], sigma = sigma)
        merge!(δ, Dict(zip(keys_d, c.shares)))
    end
    return δ
end

"""
$(TYPEDSIGNATURES)

The [`Armington`](@ref) specification of `product` that reproduces a base-year
trade matrix exactly: its shares, through every level of the origin groups
`nests`, and the `scale` of each destination's composite.

`flows` maps `(origin, destination)` pairs to the base-year quantity shipped,
the diagonal included — what each region keeps of its own production, which
trade statistics do not report and production minus exports gives. `prices`
maps the same pairs to the **delivered** price: the producer price at the
origin plus the cost of carrying the product to the destination, plus any
duty charged on the way. `sigma` and `nests` are the structure the
specification will have.

Each destination's tree is calibrated from the bottom up, by
[`ces_calibration`](@ref) at every node: a group's price index and composite
quantity make it one member of the group above, and an origin's share is the
product of the shares on its path to the root. At the root, the scale is the
base-year sum of the flows over the composite quantity, so that the composite
counts the physical units the flows are measured in and its price is the
average delivered unit value,

```math
P^c_r = \\frac{\\sum_o X_{o,r}\\, q_{o,r}}{\\sum_o X_{o,r}},
```

which is the consumer price the base year implies. Calibrated on that price,
the destination's demand curves and the input bundles of its processes buy
exactly the base-year quantities.

A group of perfect substitutes (`sigma = Inf`) is a pool, bought only at its
lowest delivered price: every member shipping to a destination must have that
price there, to within `rtol`, and the members' shares do not matter, only
their sum. An origin of a pool that ships nothing to a destination is kept in
its composite if `prices` gives its delivered price there, so that it can
start shipping when it becomes the cheapest, and the calibration checks that
it is not the cheapest already. Any other origin that ships nothing is left
out of that destination for good, as [`armington_shares`](@ref) does.

Every origin given a share must be able to deliver to its destination in the
economy the specification goes into — produce the product, and have a route
there — or the engine drops it and the shares no longer add up to what was
calibrated.

# Example
```julia
## France pooled, then with the rest of the EU, then against the rest
flows  = Dict(("FR", "EU") => 4.0,  ("EU", "EU") => 70.0, ("US", "EU") => 6.0)
prices = Dict(("FR", "EU") => 214.0, ("EU", "EU") => 205.0, ("US", "EU") => 238.0)
armington_calibration(product = "sawnwood", flows = flows, prices = prices,
                      sigma = 3, nests = [OriginNest(sigma = 9, origins = ["FR", "EU"])])
```
"""
function armington_calibration(; product::Name, flows::AbstractDict,
                                 prices::AbstractDict, sigma::Real,
                                 nests::Vector{OriginNest} = OriginNest[],
                                 rtol::Real = 1e-8)
    0 <= sigma < Inf || throw(ArgumentError("sigma must be finite and ≥ 0, got $sigma"))
    flows, prices = namekeys2(flows), namekeys2(prices)
    flow(o, r) = get(flows, (o, r), 0.0)
    grouped = leaves(nests)

    shares = Dict{Tuple{String,String},Float64}()
    scale  = Dict{String,Float64}()
    for r in unique(k[2] for (k, x) in flows if x > 0)
        delivered(o) = haskey(prices, (o, r)) ? prices[(o, r)] : throw(ArgumentError(
            "no delivered price of $product given for the flow $o → $r"))

        # A node of the tree is calibrated into its price, its composite
        # quantity and the share of each origin beneath it within the node, or
        # `nothing` when nothing beneath it is bought. Inside a pool an origin
        # that is not bought stays as long as it has a price.
        function calibrate(o::String, pooled::Bool)
            flow(o, r) > 0 && return (price = delivered(o), quantity = flow(o, r),
                                      within = Dict(o => 1.0))
            pooled && haskey(prices, (o, r)) &&
                return (price = prices[(o, r)], quantity = 0.0, within = Dict(o => 1.0))
            return nothing
        end
        function calibrate(n::OriginNest, pooled::Bool)
            inside = pooled || isinf(n.sigma)
            members = filter(!isnothing, Any[[calibrate(o, inside) for o in n.origins];
                                             [calibrate(sub, inside) for sub in n.nests]])
            return combine(members, n.sigma)
        end
        function combine(members, σ)
            # a group nothing is bought from is dropped, as the engine drops it,
            # and a group of one is that one member
            any(m -> m.quantity > 0, members) || return nothing
            length(members) == 1 && return members[1]
            c = ces_calibration(quantities = [m.quantity for m in members],
                                prices = [m.price for m in members], sigma = σ, rtol = rtol)
            # in a pool the shares are never consulted, only their sum, and a
            # member not bought must keep a positive one to stay in it
            δ = isinf(σ) ? fill(1 / length(members), length(members)) : c.shares
            within = Dict{String,Float64}()
            for (m, s) in zip(members, δ), (o, w) in m.within
                within[o] = s * w
            end
            return (price = c.price, quantity = c.quantity, within = within)
        end

        top = filter(!isnothing, Any[[calibrate(n, false) for n in nests];
                                     [calibrate(o, false) for o in unique(k[1] for k in keys(flows))
                                      if !(o in grouped)]])
        root = combine(top, sigma)
        for (o, w) in root.within
            shares[(o, r)] = w
        end
        scale[r] = sum(x for (k, x) in flows if k[2] == r && x > 0) / root.quantity
    end
    return Armington(product = product, sigma = sigma, shares = shares, nests = nests,
                     scale = scale)
end

"""
$(TYPEDSIGNATURES)

A trade matrix for one product where the statistics do not give one: the
flows that follow a gravity pattern as closely as possible while adding up to
everything that *is* known.

`production` and `use` give, by region, how much was made and how much was
used there — consumed plus taken in by its mills. Every region ships what it
makes somewhere, itself included, and receives what it uses from somewhere,
so these are the row and column totals of the matrix. They seldom agree
exactly; the uses are scaled to the total production so that they do. The
known totals below are met as given, so one that fixes the whole use of a
region has to agree with that use once scaled.

The pattern is the gravity one: a flow between two regions falls off as
``e^{-\\text{decay} \\cdot d}`` with their `distance` ``d``, a function
`(origin, destination) -> d` or a dictionary keyed `(origin, destination)` in
either order, and a region's deliveries to itself are weighted by
`home_bias`.

`totals` adds whatever else is known, as `(origins, destinations) => quantity`,
each side one region or a vector of them: a single observed flow,
`("NA", "EU") => 9.6`, or a total over a group, such as national exports that
customs report without saying which region of the country they left from,
`(["SEF", "SWF", "GEF", "OF"], "EU") => 12.5`.

The flows are found by iterative proportional fitting — scaling them in turn
to meet each total, until all are met to within `rtol` — which gives the
matrix closest to the gravity pattern, in the entropy sense, that meets every
total. A set of totals that contradict each other, such as exports larger
than production, cannot all be met, and is an error after `maxiter` rounds.

The result, diagonal included, is what [`armington_calibration`](@ref) takes as
`flows`, alone or merged with the flows that were observed.

# Example
```julia
## logs among the regions of France, from what each harvested and sawed and
## from national exports to the rest of the EU
france = ["SEF", "SWF", "GEF", "OF"]
impute_flows(production = Dict("SEF" => 3.6, "SWF" => 8.9, "GEF" => 3.0, "OF" => 3.7),
             use        = Dict("SEF" => 2.4, "SWF" => 5.1, "GEF" => 1.6, "OF" => 1.9,
                               "EU" => 8.2),
             distance   = (o, r) -> o == r ? 0.0 : 0.1,
             decay = 1.4, home_bias = 3.0,
             totals = Dict((france, "EU") => 8.2))
```
"""
function impute_flows(; production::AbstractDict, use::AbstractDict, distance,
                        decay::Real, home_bias::Real = 1.0,
                        totals::AbstractDict = Dict{Tuple{String,String},Float64}(),
                        rtol::Real = 1e-10, maxiter::Integer = 10_000)
    decay >= 0 || throw(ArgumentError("decay must be ≥ 0, got $decay"))
    home_bias > 0 || throw(ArgumentError("home_bias must be positive, got $home_bias"))
    made, used = namekeys(production), namekeys(use)
    all(>=(0), values(made)) && all(>=(0), values(used)) || throw(ArgumentError(
        "production and use must not be negative"))
    sum(values(made); init = 0.0) > 0 && sum(values(used); init = 0.0) > 0 ||
        throw(ArgumentError("nothing is produced or nothing is used"))
    d(o, r) = o == r ? 0.0 :
        distance isa AbstractDict ?
            (haskey(distance, (o, r)) ? distance[(o, r)] : haskey(distance, (r, o)) ?
             distance[(r, o)] : throw(ArgumentError("no distance given between $o and $r"))) :
        distance(o, r)

    # the cells of the matrix, and the gravity pattern they start from
    origins      = sort([o for (o, q) in made if q > 0])
    destinations = sort([r for (r, q) in used if q > 0])
    cells = [(o, r) for o in origins for r in destinations]
    x = [exp(-decay * d(o, r)) * (o == r ? home_bias : 1.0) for (o, r) in cells]

    # every total, as the cells it sums over and the quantity it must reach
    scaled = sum(values(made)) / sum(values(used))
    constraints = Tuple{String,Vector{Int},Float64}[]
    for o in origins
        push!(constraints, ("the production of $o", findall(c -> c[1] == o, cells), made[o]))
    end
    for r in destinations
        push!(constraints, ("the use of $r", findall(c -> c[2] == r, cells), scaled * used[r]))
    end
    side(s) = s isa Name ? [nm(s)] : nms(s)
    for ((from, to), q) in totals
        q >= 0 || throw(ArgumentError("a known total must not be negative, got $q"))
        os, rs = side(from), side(to)
        covered = findall(c -> c[1] in os && c[2] in rs, cells)
        isempty(covered) && q > 0 && throw(ArgumentError(
            "the total from $(join(os, ", ")) to $(join(rs, ", ")) covers no region " *
            "that produces and no region that uses the product"))
        push!(constraints, ("the flows from $(join(os, ", ")) to $(join(rs, ", "))",
                            covered, float(q)))
    end

    miss(idx, q) = abs(sum(x[idx]; init = 0.0) - q) / max(q, eps())
    for _ in 1:maxiter
        for (_, idx, q) in constraints
            s = sum(x[idx]; init = 0.0)
            s > 0 && (x[idx] .*= q / s)
        end
        maximum(miss(idx, q) for (_, idx, q) in constraints) <= rtol &&
            return Dict(zip(cells, x))
    end
    worst = argmax(c -> miss(c[2], c[3]), constraints)
    throw(ArgumentError("the known totals contradict each other: after $maxiter rounds " *
                        "$(worst[1]) is still $(round(100 * miss(worst[2], worst[3]); sigdigits = 3))% " *
                        "away from its value"))
end

"""
$(TYPEDEF)

A tax or a subsidy — a subsidy being a negative rate — driving a wedge between
what one side of a transaction pays and what the other receives.

`on` says what the wedge falls on:

| `on` | the wedge sits between | selected by |
|:---|:---|:---|
| `"supply"` | the market price and what primary producers receive | `product`, `region` |
| `"demand"` | what consumers pay and the market price | `product`, `region` |
| `"process"` | a transformation and its value-added cost | `process`, `region` |
| `"trade"` | the price at the origin and the cost of landing a good | `product`, `from`, `to` |

`rate` is per unit of the base when `relative` is `false` — per tonne produced, per
unit of activity, per tonne shipped — and a fraction of it when `relative` is
`true`: of the price received for `"supply"`, of the market price for `"demand"`
— so a rate of 0.25 puts what consumers pay a quarter above the market price —
of the value-added cost for `"process"`, of the price at the origin for `"trade"`.
Selectors left at `:` match everything, and several taxes may fall on the
same base, in which case the absolute rates add and so do the relative ones,
whose sum must stay above -1 as each rate does.

!!! warning "A taxed equilibrium is not a welfare maximum"
    The quantities and prices a taxed model returns are the right ones: the
    wedges enter the program exactly as a cost would, and its first-order
    conditions are the distorted market's. The **objective value** is not
    welfare any more, because it treats tax revenue as a cost consumed rather
    than as a transfer received. To get back to welfare, add the revenue
    reported in `Results.taxes` to the objective.

!!! note "Ad valorem rates on trade are solved by iteration"
    An ad valorem tariff is a fraction of the price at the origin, which is a
    dual of the problem and so cannot be written into its objective. For those,
    and only those, [`solve_market`](@ref) solves repeatedly, each time holding
    the tariff at the money equivalent implied by the previous solution's
    prices, until the two agree. Every other combination is exact in one solve.

# Fields
$(TYPEDFIELDS)

# Example
```julia
Tax(name = "oil_royalty",   on = "supply",  product = "crude", region = "US", rate = 4.0)
Tax(name = "fuel_vat",      on = "demand",  product = "gasoline", rate = 0.2, relative = true)
Tax(name = "green_steel",   on = "process", process = "steel_mill", region = "EU", rate = -15.0)
Tax(name = "steel_tariff",  on = "trade",   product = "steel", to = "US", rate = 0.25,
    relative = true)
```
"""
Base.@kwdef struct Tax
    "Name of the tax, used to report its revenue"
    name::String = "tax"
    "What the wedge falls on: `\"supply\"`, `\"demand\"`, `\"process\"` or `\"trade\"`"
    on::String
    "Size of the wedge, negative for a subsidy"
    rate::Float64
    "Whether `rate` is a fraction of the base rather than an amount per unit"
    relative::Bool = false
    "Products it falls on, or `:` for every one"
    product::Selector = (:)
    "Regions it falls on, or `:` for every one — the producing or consuming region"
    region::Selector = (:)
    "Processes it falls on, or `:` for every one"
    process::Selector = (:)
    "Origins, for a tax on trade"
    from::Selector = (:)
    "Destinations, for a tax on trade"
    to::Selector = (:)

    function Tax(name, on, rate, relative, product, region, process, from, to)
        base = nm(on)
        base in ("supply", "demand", "process", "trade") || throw(ArgumentError(
            "a tax falls on \"supply\", \"demand\", \"process\" or \"trade\", got $on"))
        !relative || rate > -1 || throw(ArgumentError(
            "a relative rate must be greater than -1, got $rate"))
        new(nm(name), base, rate, relative, selector(product), selector(region),
            selector(process), selector(from), selector(to))
    end
end

"""
$(TYPEDSIGNATURES)

The absolute and relative rates of the taxes of `taxes` that fall on a base,
summed. `kind` is the `on` field to match and the keywords are the selectors to
check, each against the field of the same name.

Each relative rate is above -1 on its own, but several can add up past it, and
the model multiplies the base by `1 + relative`: at or below zero a cost turns
into a gain, or the benefit of consumers into a loss, and the solver may still
report success. That is an `ArgumentError` here, where every summed rate the
model uses is computed, and [`validate`](@ref) asks for the rates of every base
so that it is reported before anything is built.
"""
function tax_rates(taxes, kind::AbstractString; product = nothing, region = nothing,
                   process = nothing, from = nothing, to = nothing)
    absolute = 0.0
    relative = 0.0
    relative_taxes = String[]
    for t in taxes
        t.on == kind || continue
        product === nothing || matches(t.product, product) || continue
        region  === nothing || matches(t.region, region)   || continue
        process === nothing || matches(t.process, process) || continue
        from    === nothing || matches(t.from, from)       || continue
        to      === nothing || matches(t.to, to)           || continue
        t.relative ? (relative += t.rate; push!(relative_taxes, "$(t.name) ($(t.rate))")) :
                     (absolute += t.rate)
    end
    relative > -1 || let base = join(("$k $v" for (k, v) in pairs((; product, region,
                                     process, from, to)) if v !== nothing), ", ")
        throw(ArgumentError("the relative rates of the taxes $(join(relative_taxes, ", ")) " *
                            "on $kind ($base) add up to $relative; together they must " *
                            "stay above -1, as each does on its own"))
    end
    return (absolute = absolute, relative = relative)
end

"""
$(TYPEDEF)

The complete description of an economy, to be passed to [`solve_market`](@ref).

# Fields
$(TYPEDFIELDS)

# Example
```julia
MarketData(regions   = ["EU", "US"],
           products  = ["crude", "gasoline"],
           demand    = demand,      # vector of DemandCurve
           supply    = supply,      # vector of SupplyCurve
           processes = processes,   # vector of Process
           transport_costs = Dict(("crude", "EU", "US") => 25.0,
                                  ("crude", "US", "EU") => 25.0))
```
"""
Base.@kwdef struct MarketData
    "Regions of the model"
    regions::Vector{String}
    "Products of the model"
    products::Vector{String}
    "Final demand curves, one [`DemandCurve`](@ref) per consumed product and region"
    demand::Vector{DemandCurve} = DemandCurve[]
    "Primary supply curves, one [`SupplyCurve`](@ref) per supplied product and region"
    supply::Vector{SupplyCurve} = SupplyCurve[]
    "Transformation activities, a vector of [`Process`](@ref)es"
    processes::Vector{Process} = Process[]
    """
    Products whose regional varieties are imperfect substitutes, as a vector of
    [`Armington`](@ref) specifications. A product that is not listed here (the
    default) is homogeneous: the same good wherever it comes from
    """
    armington::Vector{Armington} = Armington[]
    """
    Linear constraints on any of the quantities of the model, as a vector of
    [`Limit`](@ref)s: the capacity of a plant, an allowable cut, a ceiling on
    what a port can ship, a resource several products share, a quota, a
    minimum share, a fixed sum
    """
    limits::Vector{Limit} = Limit[]
    """
    Taxes and subsidies, as a vector of [`Tax`](@ref)es, on production, on
    consumption, on transformation or on trade
    """
    taxes::Vector{Tax} = Tax[]
    """
    Unit transport costs, mapping `(product, from_region, to_region)` to the
    cost per unit of quantity shipped. This is also what declares trade: a
    route that is missing from this dictionary cannot be used, and a product
    with no route at all is not traded. Routes are directional, so a product
    that moves both ways needs an entry each way
    """
    transport_costs::Dict{Tuple{String,String,String},Float64} =
        Dict{Tuple{String,String,String},Float64}()

    function MarketData(regions, products, demand, supply, processes, armington,
                        limits, taxes, transport_costs)
        new(nms(regions), nms(products), demand, supply, processes, armington,
            limits, taxes, namekeys3(transport_costs))
    end
end

"The regions process `p` operates in, in the economy `d`."
regions_of(p::Process, d::MarketData) = each(p.regions, d.regions)
