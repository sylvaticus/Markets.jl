using Documenter
using Markets

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "Markets.jl documentation",
    author = "Antonello Lobianco",
    pages = [
        "Index" => "index.md",
        "Model description" => "model_description.md",
        "Reference API" => "markets.md"
    ],
    format = Documenter.HTML(prettyurls = false),
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/sylvaticus/Markets.jl.git",
    devbranch = "main"
)
