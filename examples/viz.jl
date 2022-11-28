push!(LOAD_PATH, "$(@__DIR__)/..")

using Quac

circ = Quac.Algorithms.QFT(4)

draw(circ)