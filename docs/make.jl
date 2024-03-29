using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

push!(LOAD_PATH, "$(@__DIR__)/..")

using Documenter
using Quac

makedocs(
    sitename = "Quac",
    pages = [
        "Home" => "index.md",
        "API Reference" =>
            ["Gates" => "api/gates.md", "Circuit" => "api/circuit.md", "Algorithms" => "api/algorithms.md"],
        "Ecosystem" => "ecosystem.md",
    ],
)

deploydocs(repo = "github.com/bsc-quantic/Quac.jl.git", devbranch = "master", push_preview = true)
