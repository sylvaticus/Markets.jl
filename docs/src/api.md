# API reference

```@docs
Markets
```

## Index

```@index
Pages = ["api.md"]
```

## Public API

The types and functions exported by the module: everything needed to describe
an economy, solve it and read the results.

```@autodocs
Modules = [Markets]
Private = false
Order   = [:constant, :type, :function, :macro]
```

## Re-exported from JuMP

`latex_formulation(model)` typesets a JuMP model, such as the one
[`jump_model`](@ref) returns, as LaTeX. It is JuMP's own function, re-exported
unchanged, and is documented in
[JuMP's API reference](https://jump.dev/JuMP.jl/stable/api/JuMP/#latex_formulation).
JuMP 1.32 drops the minus sign of some negative coefficients when it does
this ([JuMP.jl#4268](https://github.com/jump-dev/JuMP.jl/issues/4268)), so
read the signs from `print(model)`.

## Internals

Not exported, and subject to change without notice.

```@autodocs
Modules = [Markets]
Public  = false
Order   = [:constant, :type, :function, :macro]
```
