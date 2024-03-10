module Algorithms

using Quac
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
            q2 = permutation[i + 1]

            push!(circuit, rand(SU{2}, q1, q2))
        end
    end

    return circuit
end

end