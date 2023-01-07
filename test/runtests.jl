using Test

include("Circuit_test.jl")

import Quac
using Aqua
Aqua.test_all(Quac, ambiguities = false, piracy = false)
Aqua.test_ambiguities(Quac)
