using Documenter
using Literate
using Markets

# The examples live in examples/ as runnable scripts with their explanation in
# comments; Literate turns each into a page whose code Documenter then
# executes, so the tables on those pages are the solver's own output.
# The three forest pages run in this order: the last reads what the first two
# saved in examples/forest/results.
const EXAMPLES = [joinpath(@__DIR__, "..", "examples", "minimal", "minimal_market.jl"),
                  joinpath(@__DIR__, "..", "examples", "forest",  "forest_market.jl"),
                  joinpath(@__DIR__, "..", "examples", "forest",  "forest_calibration.jl"),
                  joinpath(@__DIR__, "..", "examples", "forest",  "forest_full_data.jl"),
                  joinpath(@__DIR__, "..", "examples", "spatequ", "spatequ.jl"),
                  joinpath(@__DIR__, "..", "examples", "egypt",   "egypt.jl")]
const GENERATED = joinpath(@__DIR__, "src", "generated")

rm(GENERATED; recursive = true, force = true)

# The logos and the favicon live in assets/imgs/logos, beside the README that
# shows them too; Documenter finds them under the names it looks for in
# docs/src/assets (logo.svg, logo-dark.svg for the dark themes).
const LOGOS = joinpath(@__DIR__, "..", "assets", "imgs", "logos")
for (from, to) in ("markets-jl-logo.svg" => "logo.svg", "markets-jl-logo-dark.svg" => "logo-dark.svg",
                   "favicon.ico" => "favicon.ico")
    cp(joinpath(LOGOS, from), joinpath(@__DIR__, "src", "assets", to); force = true)
end
for example in EXAMPLES
    Literate.markdown(example, GENERATED; documenter = true,
                      repo_root_url = "https://github.com/sylvaticus/Markets.jl/blob/main")
end

makedocs(
    sitename = "Markets.jl",
    authors  = "Antonello Lobianco",
    modules  = [Markets],
    pages = [
        "Home"                => "index.md",
        "Using the module"    => "usage.md",
        "Examples" => [
            "Two regions and a dairy (minimal)" => "generated/minimal_market.md",
            "Forest sector economics (advanced)" => [
                "Writing the model"        => "generated/forest_market.md",
                "Calibrating with partial data" => "generated/forest_calibration.md",
                "Complete data, and what calibration changes" => "generated/forest_full_data.md",
            ],
            "Ports of published models" => [
                "Takayama and Judge (GAMS spatequ)" => "generated/spatequ.md",
                "Egyptian agriculture (GAMS egypt)" => "generated/egypt.md",
            ],
        ],
        "Modelling choices"   => [
            "modelling/index.md",
            "Equilibrium, curves and welfare" => "modelling/equilibrium.md",
            "Processes and substitution"      => "modelling/processes.md",
            "Taxes and trade"                 => "modelling/trade.md",
            "Accounts, calibration and scope" => "modelling/scope.md",
        ],
        "Code implementation" => [
            "implementation.md",
            "API reference" => "api.md",
        ],
        "References"          => "references.md",
    ],
    # numbered section headings, page by page (see the stylesheet); the logo
    # already carries the name, so the sidebar does not repeat it
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                             assets = ["assets/numbering.css", "assets/favicon.ico"],
                             sidebar_sitename = false),
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
