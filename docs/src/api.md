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

## Internals

Not exported, and subject to change without notice.

```@autodocs
Modules = [Markets]
Public  = false
Order   = [:constant, :type, :function, :macro]
```
