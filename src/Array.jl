import LinearAlgebra: Matrix, Diagonal
import LinearAlgebra

Matrix(x::AbstractGate) = Matrix{ComplexF32}(x)

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

Matrix{T}(::Type{Swap}) = Matrix{T}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])

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

# TODO Matrix{T}(g::Control) = ...

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
Diagonal{T}(::Type{T}) where {T} = Diagonal{T}([1, cispi(1 // 4)])
Diagonal{T}(::Type{Td}) where {T} = Diagonal{T}([1, cispi(-1 // 4)])

Diagonal{T}(g::Rz) where {T} = Diagonal{T}([1 0; 0 cis(g[:θ])])

# permutational matrices (diagonal + permutation)
# Permutation(x::AbstractGate) = Permutation{ComplexF32}(x)
# Permutation{T}(_::AbstractGate) where {T} = error("Implementation not found")