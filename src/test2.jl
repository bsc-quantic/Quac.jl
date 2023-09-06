using LinearAlgebra

# Function to check if a matrix is one-qubit
function is_one_qubit(U)
    return size(U) == (2, 2)
end

# Function to decompose a one-qubit gate into ZYZ rotation
function zyz_decomposition(U)
    # Decompose U into ZYZ rotation angles
    α = angle(U[1, 1] + U[2, 2])
    β = angle(-U[1, 2] + U[2, 1])
    γ = angle(U[1, 2] + U[2, 1])
    δ = angle(U[1, 1] - U[2, 2])

    θ = acos(real(0.5 * (trace(U) - 1)))
    φ = 0.5 * (β + γ)
    λ = 0.5 * (γ - β)

    return (θ, φ, λ)
end

# Recursive function to decompose a unitary matrix into one-qubit gates and CNOT gates
using LinearAlgebra

# Function to check if a matrix is one-qubit
function is_one_qubit(U)
    return size(U) == (2, 2)
end

# Function to decompose a one-qubit gate into ZYZ rotation
function zyz_decomposition(U)
    # Decompose U into ZYZ rotation angles
    α = angle(U[1, 1] + U[2, 2])
    β = angle(-U[1, 2] + U[2, 1])
    γ = angle(U[1, 2] + U[2, 1])
    δ = angle(U[1, 1] - U[2, 2])

    θ = acos(real(0.5 * (trace(U) - 1)))
    φ = 0.5 * (β + γ)
    λ = 0.5 * (γ - β)

    return (θ, φ, λ)
end

# Function to find the Gray code for the control and target qubits
function gray_code(control, target)
    return xor(control, target)
end

# Recursive function to decompose a unitary matrix into one-qubit gates and CNOT gates
function decompose(U)
    if is_one_qubit(U)
        # If U is a one-qubit gate, decompose it using the ZYZ method
        angles = zyz_decomposition(U)
        return [("U", angles)]
    else
        n = Int(log2(size(U)[1]))

        # Divide U into 4 smaller unitary matrices
        k = Int(n/2)
        U00, U01, U10, U11 = [U[1:2^k, 1:2^k], U[1:2^k, 2^k+1:end],
                              U[2^k+1:end, 1:2^k], U[2^k+1:end, 2^k+1:end]]

        # Recursive QSD
        V0 = decompose(U00 * U11 + U01 * U10)
        V1 = decompose(U00 * U10 - U01 * U11)

        # Reconstruct the CNOT gates sequence
        CNOT_sequence = []
        for i in 1:k
            control = gray_code(i - 1, i) - 1
            target = i - 1
            push!(CNOT_sequence, ("CNOT", control, target))
        end

        # Combine the decompositions of V0 and V1 with the CNOT sequence
        combined_decomposition = []
        for i in 1:length(V0)
            push!(combined_decomposition, V0[i])
            push!(combined_decomposition, CNOT_sequence...)
            push!(combined_decomposition, V1[i])
            push!(combined_decomposition, CNOT_sequence...)
        end

        return combined_decomposition
    end
end
