using Test
using Quac

@testset "Unit tests" verbose = true begin
    include("Operator_test.jl")
    include("Gate_test.jl")
    include("Circuit_test.jl")
end

using Aqua
Aqua.test_all(Quac, piracy = false, stale_deps = false)
Aqua.test_ambiguities(Quac)
