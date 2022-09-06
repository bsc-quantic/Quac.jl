if pwd() == @__DIR__
    # executing inside 'docs/' folder => load parent
    using Pkg
    Pkg.activate(".")

    push!(LOAD_PATH, dirname(pwd()))

elseif pwd() == dirname(@__DIR__)
    # executing in project root folder => activate 'docs' env
    using Pkg
    Pkg.activate("docs")

    push!(LOAD_PATH, pwd())

end

using Documenter
using Quac

makedocs(
    sitename="Quac",
    pages=[
        "index.md",
        "API Reference" => [
            "Gates" => "api/gates.md",
            "Circuit" => "api/circuit.md",
            "Algorithms" => "api/algorithms.md",
        ],
        "Ecosystem" => "ecosystem.md",
    ]
)