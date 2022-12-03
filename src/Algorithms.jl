module Algorithms

using Quac

export QFT

"""
    QFT(n)

Generate a Quantum Fourier Transform circuit of `n` qubits.
"""
QFT(n::Int) = begin
    circuit = Circuit(n)
    for target in 1:n
        push!(circuit, H(target))

        for control in target+1:n
            m = control - target
            θ = 2π / 2^m
            push!(circuit, CRz(control, target, (θ = θ,)))
        end
    end

    circuit
end

end