import LinearAlgebra: Matrix, Diagonal, eigvals, eigvecs, eigen
using LinearAlgebra: Eigen, LinearAlgebra

# preferred representation
function arraytype end
export arraytype

arraytype(::T) where {T<:AbstractGate} = arraytype(T)
arraytype(::Type{G}) where {G<:AbstractGate} = Array{T,2 * lanes(G)} where {T}

for G in [I, Z, S, Sd, T, Td, Rz]
    @eval arraytype(::Type{$G}) = Diagonal
end

# TODO Array{T,N} where {T} instead of Matrix
# TODO N-dim Diagonal type
arraytype(::Type{T}) where {T<:Control} = arraytype(op(T)) == Diagonal ? Diagonal : Matrix

Matrix(x::AbstractGate) = Matrix{ComplexF32}(x)
Matrix(::Type{T}) where {T<:AbstractGate} = Matrix{ComplexF32}(T)

for G in [I, X, Y, Z, H, S, Sd, T, Td, Swap]
    @eval Matrix{T}(_::$G) where {T} = Matrix{T}($G)
end

Matrix{T}(::Type{I}) where {T} = Matrix{T}(LinearAlgebra.I, 2, 2)
Matrix{T}(::Type{X}) where {T} = Matrix{T}([0 1; 1 0])
Matrix{T}(::Type{Y}) where {T} = Matrix{T}([0 -1im; 1im 0])
Matrix{T}(::Type{Z}) where {T} = Matrix{T}([1 0; 0 -1])
Matrix{T}(::Type{H}) where {T} = Matrix{T}([1 1; 1 -1] ./ sqrt(2))
Matrix{T}(::Type{S}) where {T} = Matrix{T}([1 0; 0 1im])
Matrix{T}(::Type{Sd}) where {T} = Matrix{T}([1 0; 0 -1im])
Matrix{F}(::Type{T}) where {F} = Matrix{F}([1 0; 0 cispi(1 // 4)])
Matrix{F}(::Type{Td}) where {F} = Matrix{F}([1 0; 0 cispi(-1 // 4)])

Matrix{T}(::Type{Swap}) where {T} = Matrix{T}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])

Matrix{T}(g::Rx) where {T} = Matrix{T}([
    cos(g[:θ] / 2) -1im*sin(g[:θ] / 2)
    -1im*sin(g[:θ] / 2) cos(g[:θ] / 2)
])
Matrix{T}(g::Ry) where {T} = Matrix{T}([
    cos(g[:θ] / 2) -sin(g[:θ] / 2)
    sin(g[:θ] / 2) cos(g[:θ] / 2)
])
Matrix{T}(g::Rz) where {T} = Matrix{T}([1 0; 0 cis(g[:θ])])

Matrix{T}(g::U2) where {T} = 1 / sqrt(2) * Matrix{T}([
    1 -cis(g[:λ])
    cis(g[:ϕ]) cis(g[:ϕ] + g[:λ])
])
Matrix{T}(g::U3) where {T} = Matrix{T}([
    cos(g[:θ] / 2) -cis(g[:λ])*sin(g[:θ] / 2)
    cis(g[:ϕ])*sin(g[:θ] / 2) cis(g[:ϕ] + g[:λ])*cos(g[:θ] / 2)
])

function Matrix{T}(::Type{Control{G}}) where {T,G}
    n = lanes(Control{G})
    t = lanes(G)

    M = Matrix{T}(LinearAlgebra.I, 2^n, 2^n)

    M[(2^n-2^t+1):end, (2^n-2^t+1):end] = Matrix{T}(G)

    return M
end

function Matrix{T}(g::Control) where {T}
    n = (length ∘ lanes)(g)

    M = Matrix{T}(LinearAlgebra.I, 2^n, 2^n)

    M[(2^n-2+1):end, (2^n-2+1):end] = Matrix{T}(op(g))

    return M
end

# diagonal matrices
# NOTE efficient multiplication due to no memory swap needed: plain element-wise multiplication
Diagonal(x::AbstractGate) = Diagonal{ComplexF32}(x)
Diagonal(::Type{T}) where {T<:AbstractGate} = Diagonal{ComplexF32}(T)

for G in [I, Z, S, Sd, T, Td]
    @eval Diagonal{T}(_::$G) where {T} = Diagonal{T}($G)
end

Diagonal{T}(::Type{I}) where {T} = Diagonal{T}(LinearAlgebra.I, 2)
Diagonal{T}(::Type{Z}) where {T} = Diagonal{T}([1, -1])
Diagonal{T}(::Type{S}) where {T} = Diagonal{T}([1, 1im])
Diagonal{T}(::Type{Sd}) where {T} = Diagonal{T}([1, -1im])
Diagonal{F}(::Type{T}) where {F} = Diagonal{F}([1, cispi(1 // 4)])
Diagonal{T}(::Type{Td}) where {T} = Diagonal{T}([1, cispi(-1 // 4)])

Diagonal{T}(g::Rz) where {T} = Diagonal{T}([1 0; 0 cis(g[:θ])])

# permutational matrices (diagonal + permutation)
# Permutation(x::AbstractGate) = Permutation{ComplexF32}(x)
# Permutation{T}(_::AbstractGate) where {T} = error("Implementation not found")

# Linear Algebra factorizations
eigen(::Type{T}) where {T<:AbstractGate} = Eigen(eigvals(T), eigvecs(T))
eigen(::T) where {T<:AbstractGate} = eigen(T)
eigen(g::T) where {T<:AbstractParametricGate} = Eigen(eigvals(T), eigvecs(T))

eigvals(::Type{T}) where {T<:AbstractGate} = eigvals(Matrix(T))
eigvals(::T) where {T<:AbstractGate} = eigvals(T)
eigvals(g::T) where {T<:AbstractParametricGate} = eigvals(Matrix(g))

eigvecs(::Type{T}) where {T<:AbstractGate} = eigvecs(Matrix(T))
eigvecs(::T) where {T<:AbstractGate} = eigvecs(T)
eigvecs(g::T) where {T<:AbstractParametricGate} = eigvecs(Matrix(g))

eigvals(::Type{I}) = [1, 1]
eigvecs(::Type{I}) = [1 0; 0 1]

for G in [X, Y, Z, H]
    @eval eigvals(::$G) = eigvals($G)
    @eval eigvals(::Type{$G}) = [-1, 1]
end

eigvecs(::Type{X}) = sqrt(2) / 2 .* [-1 1; 1 1]
eigvecs(::Type{Y}) = sqrt(2) / 2 .* [-1im -1im; -1 1]
eigvecs(::Type{Z}) = [0 1; 1 0]
eigvecs(::Type{H}) = (m = [1-sqrt(2) 1+sqrt(2); 1 1]; m ./ sqrt.(sum(m .^ 2, dims = 1)))

eigvals(::Type{S}) = [1im, 1]
eigvals(::Type{Sd}) = [-1im, 1]
eigvals(::Type{T}) = [sqrt(2) / 2 + 1im * sqrt(2) / 2, 1]
eigvals(::Type{Td}) = [sqrt(2) / 2 - 1im * sqrt(2) / 2, 1]

for G in [S, Sd, T, Td]
    @eval eigvecs(::Type{$G}) = eigvecs(Z)
end

for G in [Rx, Ry, Rz]
    @eval eigen(g::$G) = Eigen(eigvals(g), eigvecs($G))
end

eigvals(g::Rx) = [cis(-g[:θ] / 2), cis(g[:θ] / 2)]
eigvals(g::Ry) = [cis(-g[:θ] / 2), cis(g[:θ] / 2)]
eigvals(g::Rz) = [1, cis(g[:θ])]

eigvecs(g::Rx) = eigvecs(Rx)
eigvecs(::Type{Rx}) = (α = sqrt(2) / 2; α * [1 1; 1 -1])
eigvecs(g::Ry) = eigvecs(Ry)
eigvecs(::Type{Ry}) = (α = sqrt(2) / 2; α * [-1im 1; 1 -1im])
eigvecs(g::Rz) = eigvecs(Rz)
eigvecs(::Type{Rz}) = [1 0; 0 1]

eigvals(::Type{Swap}) = [-1, 1, 1, 1]
eigvecs(::Type{Swap}) = (α = sqrt(2) / 2;
[0 1 0 0; α 0 α 0; -α 0 α 0; 0 0 0 1])

# TODO eigen of `Control`?
