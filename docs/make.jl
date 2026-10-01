using Documenter
using Markets

makedocs(
    sitename = "Markets.jl",
    authors  = "Antonello Lobianco",
    modules  = [Markets],
    pages = [
        "Home"                => "index.md",
        "Using the module"    => "usage.md",
        "Modelling choices"   => "modelling.md",
        "Code implementation" => [
            "implementation.md",
            "API reference" => "api.md",
        ],
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
)

# Deploy to gh-pages (only acts when run on CI)
deploydocs(
    repo = "github.com/sylvaticus/Markets.jl.git",
    devbranch = "main"
)
