# using LinearAlgebra
# using Random

# function to_su2(u)
#     return sqrt(1 / det(u)) * u
# end

# function trace_dist(u, v)
#     return real(0.5 * tr(sqrt((u - v)' * (u - v))))
# end

# function create_unitaries(base, limit)
#     gate_list = []
#     for width in 1:limit
#         for bits in bitprod(width)
#             u = Matrix{ComplexF32}(LinearAlgebra.I, 2, 2)
#             for bit in bits
#                 u = u * base[bit]
#             end
#             push!(gate_list, u)
#         end
#     end
#     return gate_list
# end

# function find_closest_u(gate_list, u)
#     min_dist, min_u = 10, Matrix{ComplexF32}(LinearAlgebra.I, 2, 2)
#     for gate in gate_list
#         tr_dist = trace_dist(gate, u)
#         if tr_dist < min_dist
#             min_dist, min_u = tr_dist, gate
#         end
#     end
#     return min_u
# end

# function u_to_bloch(u)
#     angle = real(acos((u[1, 1] + u[2, 2]) / 2))
#     sin_val = sin(angle)
#     if sin_val < 1e-10
#         axis = [0, 0, 1]
#     else
#         nx = (u[1, 2] + u[2, 1]) / (2im * sin_val)
#         ny = (u[1, 2] - u[2, 1]) / (2 * sin_val)
#         nz = (u[1, 1] - u[2, 2]) / (2im * sin_val)
#         axis = [nx, ny, nz]
#     end
#     return real.(axis), 2 * angle
# end

# function gc_decomp(u)
#     function diagonalize(u)
#         _, v = eigen(u)
#         return Matrix(v)
#     end

#     axis, theta = u_to_bloch(u)

#     phi = 2.0 * asin(sqrt(sqrt((0.5 - 0.5 * cos(theta / 2)))))

#     v = exp(-1im * phi / 2 * [0 -1; 1 0])
#     if real(axis[3]) > 0
#         w = exp(-1im * (2 * pi - phi) / 2 * [0 -1; 1 0])
#     else
#         w = exp(-1im * phi / 2 * [0 -1; 1 0])
#     end

#     ud = diagonalize(u)
#     vwvdwd = diagonalize(v * w * v' * w')
#     s = ud * vwvdwd'

#     v_hat = s * v * s'
#     w_hat = s * w * s'
#     return v_hat, w_hat
# end

# function sk_algo(u, gates, n)
#     if n == 0
#         return find_closest_u(gates, u)
#     else
#         u_next = sk_algo(u, gates, n - 1)
#         v, w = gc_decomp(u * u_next')
#         v_next = sk_algo(v, gates, n - 1)
#         w_next = sk_algo(w, gates, n - 1)
#         return v_next * w_next * v_next' * w_next' * u_next
#     end
# end


using LinearAlgebra

function gray_code(n)
    if n == 1
        return ["0", "1"]
    else
        prev_gray = gray_code(n - 1)
        return vcat([g * "0" for g in prev_gray], [g * "1" for g in reverse(prev_gray)])
    end
end

function controlled_gray_code(n)
    gc = gray_code(n)
    result = Matrix{Bool}(undef, length(gc), n)
    for (i, g) in enumerate(gc)
        result[i, :] = [c == '1' for c in g]
    end
    return result
end

function isometry_decomposition(U)
    m, n = size(U)
    println("size U: $m x $n")

    input_qubits = Int(log2(m))
    output_qubits = Int(log2(n))

    if m == 1 || output_qubits == 1
        return [([], U)]
    end

    gc = controlled_gray_code(output_qubits - 1)
    gates = []
    V = U
    for k in 1:length(gc)
        control = findall(x -> x == 1, gc[k])
        if isempty(control)
            gate = Matrix{Complex{Float64}}(LinearAlgebra.I, n, n)
            gate[1:m, 1] = V[:, 1]
            push!(gates, ([], gate))
            V = V[:, 2:end]
        else
            c = control[1]
            W_rows = 2^(input_qubits - 1)
            W_cols = 2^c
            W = V[1:W_rows, 1:W_cols]
            new_gates = isometry_decomposition(W)
            for g in new_gates
                push!(gates, (control, g))
            end
            g = new_gates[end] # The last gate in the new_gates list

            I_block = Matrix{Complex{Float64}}(LinearAlgebra.I, 2^(c - 1), 2^(c - 1))
            V1 = V[:, 1:W_cols]
            V2 = V[:, W_cols+1:end]

            V = hcat(V2, V1 * g' * kron(I_block, Matrix{Complex{Float64}}(LinearAlgebra.I, W_rows, W_rows)))
        end
    end
    return gates
end

_X = [0 1; 1 0]
_Y = [0 -im; im 0]
_Z = [1 0; 0 -1]
_H = 1/sqrt(2) * [1 1; 1 -1]
_I2 = Matrix{ComplexF64}(LinearAlgebra.I, 2, 2)
_CNOT = Matrix{ComplexF64}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0])
function unitary_distance(U1::Matrix{ComplexF64}, U2::Matrix{ComplexF64})
    return norm(U1 - U2)
end
function generate_gate_sequences(length::Int)
    gates = [ kron(_I2, _X), kron(_I2, _Y), kron(_I2, _Z), kron(_I2, _H),kron(_X,_I2), kron( _Y, _I2), kron( _Z, _I2), kron( _H, _I2), _CNOT]
    seqs = [[]]

    for _ in 1:length
        new_seqs = []
        for seq in seqs
            for gate in gates
                push!(new_seqs, [seq..., gate])
            end
        end
        seqs = new_seqs
    end

    return seqs
end
function solovay_kitaev(U::Matrix{ComplexF64}, depth::Int)
    if depth == 0
        best_sequence = []
        best_distance = Inf
        for seq in generate_gate_sequences(1)
            V = reduce(*, seq)
            dist = unitary_distance(U, V)
            if dist < best_distance
                best_distance = dist
                best_sequence = seq
            end
        end
        return best_sequence
    else
        # Apply the Solovay-Kitaev algorithm recursively
        approx_seq = solovay_kitaev(U, depth - 1)
        approx_U = reduce(*, approx_seq)
        U_diff = inv(approx_U) * U

        # Find the best correction sequence
        best_corr_seq = []
        best_corr_distance = Inf
        for seq in generate_gate_sequences(1)
            V = reduce(*, seq)
            dist = unitary_distance(U_diff, V)
            if dist < best_corr_distance
                best_corr_distance = dist
                best_corr_seq = seq
            end
        end

        # Return the combined sequence
        combined_sequence = [approx_seq..., best_corr_seq...]
        return combined_sequence
    end
end


using LinearAlgebra
using TensorOperations

function qsd(U::Matrix{ComplexF64})
    n = Int(log2(size(U, 1)))
    m = 2^n

    # Initialize the output matrices
    V = Array{Matrix{ComplexF64}, 1}(undef, n)
    W = Array{Matrix{ComplexF64}, 1}(undef, n)

    for k = 1:n
        V[k] = Matrix{ComplexF64}(LinearAlgebra.I, 2, 2)
        W[k] = Matrix{ComplexF64}(LinearAlgebra.I, 2, 2)
    end

    V, W, _ = qsd_recursive(U, V, W, n, m)

    return V, W
end

function qsd_recursive(U::Matrix{ComplexF64}, V::Array{Matrix{ComplexF64}, 1}, W::Array{Matrix{ComplexF64}, 1}, n::Int, m::Int)
    if n == 1
        V[1] = U[1:2, 1:2]
        W[1] = U[1:2, 3:4]
        return V, W, U
    else
        k = 2^(n - 1)
        U11 = U[1:k, 1:k]
        U12 = U[1:k, (k+1):m]
        U21 = U[(k+1):m, 1:k]
        U22 = U[(k+1):m, (k+1):m]

        # The recursive step
        V, W, R11 = qsd_recursive(U11, V, W, n-1, k)
        V, W, R12 = qsd_recursive(U12, V, W, n-1, k)
        V, W, R21 = qsd_recursive(U21, V, W, n-1, k)
        V, W, R22 = qsd_recursive(U22, V, W, n-1, k)

        # Combine the results
        R = zeros(ComplexF64, m, m)
        R[1:k, 1:k] = R11
        R[1:k, (k+1):m] = R12
        R[(k+1):m, 1:k] = R21
        R[(k+1):m, (k+1):m] = R22

        # Tensor product of 2x2 identity matrix
        I2 = Matrix{ComplexF64}(LinearAlgebra.I, 2, 2)
        @tensoropt V[n][a,b] := I2[a,c] * V[n-1][b,d] * conj(I2[c,d])
        @tensoropt W[n][a,b] := I2[a,c] * W[n-1][b,d] * conj(I2[c,d])

        return V, W, R
    end
end
