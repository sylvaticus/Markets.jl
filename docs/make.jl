using Documenter
using Literate
using Markets

# The forest example lives in examples/ as a runnable script with its
# explanation in comments; Literate turns it into a page whose code Documenter
# then executes, so the tables on the page are the solver's own output.
const EXAMPLE   = joinpath(@__DIR__, "..", "examples", "forest", "forest_market.jl")
const GENERATED = joinpath(@__DIR__, "src", "generated")

rm(GENERATED; recursive = true, force = true)
Literate.markdown(EXAMPLE, GENERATED; documenter = true,
                  repo_root_url = "https://github.com/sylvaticus/Markets.jl/blob/main")

makedocs(
    sitename = "Markets.jl",
    authors  = "Antonello Lobianco",
    modules  = [Markets],
    pages = [
        "Home"                => "index.md",
        "Using the module"    => "usage.md",
        "Forest example"      => "generated/forest_market.md",
        "Modelling choices"   => "modelling.md",
        "Code implementation" => [
            "implementation.md",
            "API reference" => "api.md",
        ],
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
)

# Deploy to gh-pages (only acts when run on CI).
#
# `versions` has to name patch releases explicitly: Documenter's default
# `"v#.#"` scheme drops every 0.0.x version, which would leave v0.0.1 and
# v0.0.2 — the two states of the old model — out of the version selector.
deploydocs(
    repo = "github.com/sylvaticus/Markets.jl.git",
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
)
