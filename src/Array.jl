import Base: Matrix, Array
import LinearAlgebra: Diagonal, diag, eigvals, eigvecs, eigen
using LinearAlgebra: Eigen, LinearAlgebra, qr

# preferred representation
function arraytype end
export arraytype

arraytype(::Gate{Op}) where {Op<:Operator} = arraytype(Op)
# arraytype(::Type{<:Gate}) = Array{T} where {T}

for Op in [I, Z, S, Sd, T, Td, Rz]
    @eval arraytype(::$Op) = Diagonal
end

for Op in [X, Y, H, Rx, Ry]
    @eval arraytype(::$Op) = Matrix
end

# TODO arraytype(::Type{T}) where {T<:Control} = arraytype(op(T)) == Diagonal ? Diagonal : Matrix

Matrix(x::Operator) = Matrix{ComplexF64}(x)

Matrix{T}(::I) where {T} = Matrix{T}(LinearAlgebra.I, 2, 2)
Matrix{T}(::X) where {T} = Matrix{T}([0 1; 1 0])
Matrix{T}(::Y) where {T} = Matrix{T}([0 -1im; 1im 0])
Matrix{T}(::Z) where {T} = Matrix{T}([1 0; 0 -1])
Matrix{T}(::H) where {T} = Matrix{T}([1 1; 1 -1] ./ sqrt(2))
Matrix{T}(::S) where {T} = Matrix{T}([1 0; 0 1im])
Matrix{T}(::Sd) where {T} = Matrix{T}([1 0; 0 -1im])
Matrix{F}(::T) where {F} = Matrix{F}([1 0; 0 cispi(1 // 4)])
Matrix{F}(::Td) where {F} = Matrix{F}([1 0; 0 cispi(-1 // 4)])

Matrix{T}(::Swap) where {T} = Matrix{T}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])

Matrix{T}(op::Rx) where {T} = Matrix{T}([
    cos(op.θ / 2) -1im*sin(op.θ / 2)
    -1im*sin(op.θ / 2) cos(op.θ / 2)
])
Matrix{T}(op::Ry) where {T} = Matrix{T}([
    cos(op.θ / 2) -sin(op.θ / 2)
    sin(op.θ / 2) cos(op.θ / 2)
])
Matrix{T}(op::Rz) where {T} = Diagonal{T}(op)

Matrix{T}(op::Rxx) where {T} = Matrix{T}(
    [
        cos(op.θ / 2) 0 0 -1im*sin(op.θ / 2)
        0 cos(op.θ / 2) -1im*sin(op.θ / 2) 0
        0 -1im*sin(op.θ / 2) cos(op.θ / 2) 0
        -1im*sin(op.θ / 2) 0 0 cos(op.θ / 2)
    ],
)

Matrix{T}(op::Ryy) where {T} = Matrix{T}(
    [
        cos(op.θ / 2) 0 0 1im*sin(op.θ / 2)
        0 cos(op.θ / 2) -1im*sin(op.θ / 2) 0
        0 -1im*sin(op.θ / 2) cos(op.θ / 2) 0
        1im*sin(op.θ / 2) 0 0 cos(op.θ / 2)
    ],
)

Matrix{T}(op::Rzz) where {T} = Matrix{T}([
    cis(-op.θ / 2) 0 0 0
    0 cis(op.θ / 2) 0 0
    0 0 cis(op.θ / 2) 0
    0 0 0 cis(-op.θ / 2)
])

Matrix{T}(op::U2) where {T} = 1 / sqrt(2) * Matrix{T}([
    1 -cis(op.λ)
    cis(op.ϕ) cis(op.ϕ + op.λ)
])

Matrix{T}(op::U3) where {T} = Matrix{T}([
    cos(op.θ / 2) -cis(op.λ)*sin(op.θ / 2)
    cis(op.ϕ)*sin(op.θ / 2) cis(op.ϕ + op.λ)*cos(op.θ / 2)
])

Matrix{T}(op::Hz) where {T} = Matrix{T}(
    [
        cis(op.θ / 2)*cos(op.θ / 2) -1im*cis(-op.θ / 2 + op.ϕ)*sin(op.θ / 2)
        -1im*cis(-op.θ / 2 - op.ϕ)*sin(op.θ / 2) cis(op.θ / 2)*cos(op.θ / 2)
    ],
)

Matrix{T}(op::FSim) where {T} = Matrix{T}([
    1 0 0 0
    0 cos(op.θ) -1im*sin(op.ϕ) 0
    0 -1im*sin(op.ϕ) cos(op.θ) 0
    0 0 0 1
])

function Matrix{T}(op::Op) where {T,Op<:Control}
    N = length(Op)
    t = length(targettype(Op))

    M = Matrix{T}(LinearAlgebra.I, 2^N, 2^N)

    M[(2^N-2^t+1):end, (2^N-2^t+1):end] = Matrix{T}(target(op))

    return M
end

function Matrix{T}(op::SU) where {T}
    return Matrix{T}(op.matrix)
end

Array(x::Operator) = Array{ComplexF64}(x)

Array{T}(op::Op) where {T,Op<:Operator} = Array{T,2 * length(Op)}(op)

# NOTE multidimensional `Array` literal concatenation was introduced in 1.7
# TODO clean code when we stop supporting Julia 1.6
Array{T,4}(::Swap) where {T} = Array{T}([1; 0;; 0; 0;;; 0; 0;; 1; 0;;;; 0; 1;; 0; 0;;; 0; 0;; 0; 1])

Array{T,4}(op::FSim) where {T} =
    Array{T}([1; 0;; 0; 0;;; 0; cos(op.θ);; -1im*sin(op.ϕ); 0;;;; 0; -1im*sin(op.ϕ);; cos(op.θ); 0;;; 0; 0;; 0; 1])

Array{T,4}(op::Rxx) where {T} = Array{T,4}(
    [
        cos(op.θ / 2); 0;; 0; -1im*sin(op.θ / 2);;; 0; cos(op.θ / 2);; -1im*sin(op.θ / 2); 0;;;;
        0; -1im*sin(op.θ / 2);; cos(op.θ / 2); 0;;; -1im*sin(op.θ / 2); 0;; 0; cos(op.θ / 2)
    ],
)

Array{T,4}(op::Ryy) where {T} = Array{T,4}(
    [
        cos(op.θ / 2); 0;; 0; 1im*sin(op.θ / 2);;; 0; cos(op.θ / 2);; -1im*sin(op.θ / 2); 0;;;;
        0; -1im*sin(op.θ / 2);; cos(op.θ / 2); 0;;; 1im*sin(op.θ / 2); 0;; 0; cos(op.θ / 2)
    ],
)

Array{T,4}(op::Rzz) where {T} = Array{T,4}(
    [
        cis(-op.θ / 2); 0;; 0; 0;;; 0; cis(op.θ / 2);; 0; 0;;;;
        0; 0;; cis(op.θ / 2); 0;;; 0; 0;; 0; cis(-op.θ / 2)
    ],
)

Array{T}(op::Op) where {T,Op<:Control} = Array{T,2 * length(Op)}(reshape(Matrix{T}(op), fill(2, 2 * length(Op))...))

Array{T}(op::SU{N}) where {T,N} = Array{T,2N}(reshape(Matrix{T}(op), fill(2, 2N)...))

# diagonal matrices
# NOTE efficient multiplication due to no memory swap needed: plain element-wise multiplication
Diagonal(x::Operator) = Diagonal{ComplexF64}(x)

Diagonal{T}(::I) where {T} = Diagonal{T}(LinearAlgebra.I, 2)
Diagonal{T}(::Z) where {T} = Diagonal{T}([1, -1])
Diagonal{T}(::S) where {T} = Diagonal{T}([1, 1im])
Diagonal{T}(::Sd) where {T} = Diagonal{T}([1, -1im])
Diagonal{F}(::T) where {F} = Diagonal{F}([1, cispi(1 // 4)])
Diagonal{T}(::Td) where {T} = Diagonal{T}([1, cispi(-1 // 4)])

Diagonal{T}(op::Rz) where {T} = Diagonal{T}([1, cis(op.θ)])
Diagonal{T}(op::Rzz) where {T} = Diagonal{T}([1, cis(op.θ), cis(op.θ), 1])

Diagonal{T}(op::Op) where {T,Op<:Control} =
    Diagonal{T}([fill(one(T), 2^length(Op) - 2^length(targettype(Op)))..., diag(Diagonal{T}(target(op)))...])

# permutational matrices (diagonal + permutation)
# Permutation(x::Gate) = Permutation{ComplexF64}(x)
# Permutation{T}(_::Gate) where {T} = error("Implementation not found")

# Linear Algebra factorizations
eigen(::Type{T}) where {T<:Gate} = Eigen(eigvals(T), eigvecs(T))
eigen(g::T) where {T<:Gate} = isparametric(T) ? Eigen(eigvals(g), eigvecs(g)) : eigen(T)

eigvals(::Type{T}) where {T<:Gate} = eigvals(Matrix(T))
eigvals(g::T) where {T<:Gate} = isparametric(T) ? eigvals(Matrix(g)) : eigvals(T)

eigvecs(::Type{T}) where {T<:Gate} = eigvecs(Matrix(T))
eigvecs(g::T) where {T<:Gate} = isparametric(T) ? eigvecs(Matrix(g)) : eigvecs(T)

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

eigvals(::Type{T}) where {T<:Control} = [eigvals(op(T))..., fill(1, 2^(lanes(T) - 1))...]
eigvals(g::T) where {T<:Control} = [eigvals(op(g))..., fill(1, 2^(lanes(T) - 1))...]

# TODO eigenvecs of `Control`?
