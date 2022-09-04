module Algorithms

using Quac

export QFT

"""
	QFT(n)

Generate a Quantum Fourier Transform circuit of `n` qubits.
"""
QFT(n::Int) = begin
    circ = Circuit(n)
    for target in 1:n
        push!(circ, H(target))

        for control in target+1:n
            m = control - target
            θ = 2π / 2^m
            push!(circ, CRz(control, target, θ))
        end
    end

    circ
end

end