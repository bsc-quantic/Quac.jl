using Test

include("Circuit_test.jl")

import Quac
using Aqua
Aqua.test_all(Quac, piracy = false, stale_deps = false)
Aqua.test_ambiguities(Quac)
