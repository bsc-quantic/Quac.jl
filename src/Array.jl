import Base: Matrix, Array
import LinearAlgebra: Diagonal, diag
using LinearAlgebra: Eigen, LinearAlgebra, qr

# preferred representation
function arraytype end
export arraytype

arraytype(gate::Gate) = arraytype(operator(gate))

arraytype(::I) = Diagonal{Bool}
arraytype(::X) = Matrix{Bool}
arraytype(::Y) = Matrix{Int}
arraytype(::Z) = Diagonal{Int}

for Op in [S, Sd, T, Td, Rz]
    @eval arraytype(::$Op) = Diagonal{ComplexF64}
end

arraytype(::H) = Matrix{Float64}

for Op in [Rx, Hz, U2, U3]
    @eval arraytype(::$Op) = Matrix{ComplexF64}
end

arraytype(::Ry) = Matrix{Float64}

for Op in [Rxx, Ryy, Rzz, FSim]
    @eval arraytype(::$Op) = Array{ComplexF64,4}
end

arraytype(::Swap) = Array{Bool,4}
arraytype(::SU{N}) where {N} = Array{ComplexF64,2N}

arraytype(op::Control) = Array{ComplexF64,2 * length(op)}

(::Type{A})(x::Gate) where {A<:Array} = A(operator(x))
Diagonal(x::Quac.Gate) = Diagonal(operator(x))
Diagonal{T}(x::Quac.Gate) where {T} = Diagonal{T}(operator(x))
Diagonal{T,V}(x::Quac.Gate) where {T,V<:AbstractVector{T}} = Diagonal{T,V}(operator(x))

# TODO arraytype(::Type{T}) where {T<:Control} = arraytype(op(T)) == Diagonal ? Diagonal : Matrix

Matrix(x::Operator) = Matrix{ComplexF64}(x)

Matrix{T}(::I) where {T} = Matrix{T}(LinearAlgebra.I, 2, 2)
Matrix{T}(::X) where {T} = Matrix{T}([0 1; 1 0])
Matrix{T}(::Y) where {T} = Matrix{T}([0 -1; 1 0])
Matrix{T}(::Z) where {T} = Matrix{T}([1 0; 0 -1])
Matrix{T}(::H) where {T} = Matrix{T}([1 1; 1 -1] ./ sqrt(2))
Matrix{T}(::S) where {T} = Matrix{T}([1 0; 0 1im])
Matrix{T}(::Sd) where {T} = Matrix{T}([1 0; 0 -1im])
Matrix{F}(::T) where {F} = Matrix{F}([1 0; 0 cispi(1 // 4)])
Matrix{F}(::Td) where {F} = Matrix{F}([1 0; 0 cispi(-1 // 4)])

Matrix{T}(op::Rx) where {T} = Matrix{T}([
    cos(op.θ / 2) -1im*sin(op.θ / 2)
    -1im*sin(op.θ / 2) cos(op.θ / 2)
])
Matrix{T}(op::Ry) where {T} = Matrix{T}([
    cos(op.θ / 2) -sin(op.θ / 2)
    sin(op.θ / 2) cos(op.θ / 2)
])
Matrix{T}(op::Rz) where {T} = Matrix{T}(Diagonal{T}(op))

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

Matrix{T}(::Swap) where {T} = Matrix{T}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])

Matrix{T}(::ISwap) where {T} = Matrix{T}([1 0 0 0; 0 0 1im 0; 0 1im 0 0; 0 0 0 1])

Matrix{T}(op::FractionalISwap) where {T} = Matrix{T}(
    [
        1 0 0 0
        0 cospi(op.exponent / 2) 1im*sinpi(op.exponent / 2)*cispi(2 * op.phase_exponent) 0
        0 1im*sinpi(op.exponent / 2)*cispi(-2 * op.phase_exponent) cospi(op.exponent / 2) 0
        0 0 0 1
    ],
)

Matrix{T}(::FSwap) where {T} = Matrix{T}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1])

Matrix{T}(op::FSim) where {T} = Matrix{T}([
    1 0 0 0
    0 cos(op.θ) -1im*sin(op.θ) 0
    0 -1im*sin(op.θ) cos(op.θ) 0
    0 0 0 cis(-op.ϕ)
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

Array{T,4}(::Swap) where {T} = Array{T}([1; 0;; 0; 0;;; 0; 0;; 1; 0;;;; 0; 1;; 0; 0;;; 0; 0;; 0; 1])

Array{T,4}(::ISwap) where {T} = Array{T}([1; 0;; 0; 0;;; 0; 0;; 1im; 0;;;; 0; 1im;; 0; 0;;; 0; 0;; 0; 1])

Array{T,4}(op::FractionalISwap) where {T} = Array{T}(
    [
        1; 0;; 0; 0;;; 0; cospi(op.exponent / 2);; 1im*sinpi(op.exponent / 2)*cispi(-2 * op.phase_exponent); 0;;;;
        0; 1im*sinpi(op.exponent / 2)*cispi(-2 * op.phase_exponent);; cospi(op.exponent / 2); 0;;; 0; 0;; 0; 1
    ],
)

Array{T,4}(::Swap) where {T} = Array{T}([1; 0;; 0; 0;;; 0; 0;; 1; 0;;;; 0; 1;; 0; 0;;; 0; 0;; 0; -1])

Array{T,4}(op::FSim) where {T} = Array{T}(
    [1; 0;; 0; 0;;; 0; cos(op.θ);; -1im*sin(op.θ); 0;;;; 0; -1im*sin(op.θ);; cos(op.θ); 0;;; 0; 0;; 0; cis(-op.θ)],
)

Array{T}(op::Op) where {T,Op<:Control} = Array{T,2 * length(Op)}(op)
Array{T,N}(op::Control) where {T,N} = Array{T,N}(reshape(Matrix{T}(op), fill(2, N)...))

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
LinearAlgebra.eigen(op::Operator) = Eigen(eigvals(op), eigvecs(op))
LinearAlgebra.eigvals(op::Operator) = eigvals(Matrix(op))
LinearAlgebra.eigvecs(op::Operator) = eigvecs(Matrix(op))

LinearAlgebra.eigen(g::Gate) = eigen(operator(g))
LinearAlgebra.eigvals(g::Gate) = eigvals(operator(g))
LinearAlgebra.eigvecs(g::Gate) = eigvecs(operator(g))

LinearAlgebra.eigvals(::I) = [1, 1]
LinearAlgebra.eigvecs(::I) = [1 0; 0 1]

for Op in [X, Y, Z, H]
    @eval LinearAlgebra.eigvals(::$Op) = [-1, 1]
end

LinearAlgebra.eigvecs(::X) = sqrt(2) / 2 .* [-1 1; 1 1]
LinearAlgebra.eigvecs(::Y) = sqrt(2) / 2 .* [-1im -1im; -1 1]
LinearAlgebra.eigvecs(::Z) = [0 1; 1 0]
LinearAlgebra.eigvecs(::H) = (m = [1-sqrt(2) 1+sqrt(2); 1 1]; m ./ sqrt.(sum(m .^ 2, dims = 1)))

LinearAlgebra.eigvals(::S) = [1im, 1]
LinearAlgebra.eigvals(::Sd) = [-1im, 1]
LinearAlgebra.eigvals(::T) = [sqrt(2) / 2 + 1im * sqrt(2) / 2, 1]
LinearAlgebra.eigvals(::Td) = [sqrt(2) / 2 - 1im * sqrt(2) / 2, 1]

for G in [S, Sd, T, Td]
    @eval LinearAlgebra.eigvecs(::$G) = eigvecs(Z())
end

LinearAlgebra.eigvals(op::Rx) = [cis(-op.θ / 2), cis(op.θ / 2)]
LinearAlgebra.eigvals(op::Ry) = [cis(-op.θ / 2), cis(op.θ / 2)]
LinearAlgebra.eigvals(op::Rz) = [1, cis(op.θ)]

LinearAlgebra.eigvecs(::Rx) = (α = sqrt(2) / 2; α * [1 1; 1 -1])
LinearAlgebra.eigvecs(::Ry) = (α = sqrt(2) / 2; α * [-1im 1; 1 -1im])
LinearAlgebra.eigvecs(::Rz) = [1 0; 0 1]

LinearAlgebra.eigvals(::Swap) = [-1, 1, 1, 1]
LinearAlgebra.eigvecs(::Swap) = (α = sqrt(2) / 2;
[0 1 0 0; α 0 α 0; -α 0 α 0; 0 0 0 1])

LinearAlgebra.eigvals(op::T) where {T<:Control} = [eigvals(op(g))..., fill(1, 2^(lanes(T) - 1))...]

# TODO eigenvecs of `Control`?
