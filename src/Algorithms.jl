module Algorithms

using Quac

export QFT

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
            push!(circuit, CRz(control, target, (θ = 2π / 2^m,)))
        end
    end

    circuit
end

end