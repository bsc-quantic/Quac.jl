module Algorithms

using Quac
using Graphs
using Random: randperm

export QFT
export QuantumVolume

"""
    QFT(n)

Generate a Quantum Fourier Transform circuit of `n` qubits.
"""
function QFT(n::Int)
    circuit = Circuit(n)
    for target in 1:n
        push!(circuit, H(target))

        for control in target+1:n
            m = control - target
            push!(circuit, CRz(control, target, θ = 2π / 2^m))
        end
    end

    circuit
end

"""
    QuantumVolume(n, depth)

Generate a Quantum Volume circuit of `n` qubits and `depth` layers. See [1] for more details.

# References

[1] Cross, Andrew W., et al. "Validating quantum computers using randomized model circuits." Physical Review A 100.3 (2019): 032328.
"""
function QuantumVolume(n, depth)
    circuit = Circuit(n)

    for _ in 1:depth
        # Generate a random permutation for this layer
        permutation = randperm(n)

        for i in 1:2:n-1
            q1 = permutation[i]
            q2 = permutation[i+1]

            push!(circuit, rand(SU{2}, q1, q2))
        end
    end

    return circuit
end

function RQC(grid::Quac.Sycamore, depth::Int; sequence = "ABCDCDAB")
    sequence ∈ ["ABCDCDAB", "ABCDABCD"] || throw(ArgumentError("Invalid sequence: $sequence"))

    circuit = Circuit(nv(grid))

    # initial layer of Hadamards
    for qid in vertices(grid)
        push!(circuit, H(grid[qid]))
    end

    randomops1 = [√X(), √Y(), U3(; θ = π, ϕ = π / 4, λ = 3π / 4)] # TODO is U3 parameters correct?

    for layer in 1:depth
        for qid in vertices(grid)
            push!(circuit, rand(randomops1)(grid[qid]))
        end

        mod = layer % length(sequence)
        target_edges = if sequence[mod] == 'A'
            Iterators.filter(Quac.hedges(grid)) do edge
                i, _ = Tuple(edge[1])
                isodd(i)
            end |> collect
        elseif sequence[mod] == 'B'
            Iterators.filter(Quac.hedges(grid)) do edge
                i, _ = Tuple(edge[1])
                iseven(i)
            end |> collect
        elseif sequence[mod] == 'C'
            Iterators.filter(Quac.vedges(grid)) do edge
                i, _ = Tuple(edge[1])
                isodd(i)
            end |> collect
        elseif sequence[mod] == 'D'
            Iterators.filter(Quac.vedges(grid)) do edge
                i, _ = Tuple(edge[1])
                iseven(i)
            end |> collect
        end

        for edge in target_edges
            i, j = Tuple(edge)
            push!(circuit, rand(FSim)(grid[i], grid[j]))
        end
    end

    # last layer of Hadamards
    for qid in vertices(grid)
        push!(circuit, H(grid[qid]))
    end

    circuit
end

end
