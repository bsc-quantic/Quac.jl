using Test
using Quac

@testset "Unit tests" verbose = true begin
    include("Circuit_test.jl")
end

using Aqua
Aqua.test_all(Quac, piracy = false, stale_deps = false)
Aqua.test_ambiguities(Quac)
