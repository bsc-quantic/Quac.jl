import LinearAlgebra: Matrix, Diagonal

Matrix(x::AbstractGate) = Matrix{ComplexF32}(x)
Matrix{T}(_::AbstractGate) where {T} = error("Implementation not found")

Matrix{T}(_::I) where {T} = Matrix{T}(LinearAlgebra.I, 2, 2)
Matrix{T}(_::X) where {T} = Matrix{T}([0 1; 1 0])
Matrix{T}(_::Y) where {T} = Matrix{T}([0 -1im; 1im 0])
Matrix{T}(_::Z) where {T} = Matrix{T}([1 0; 0 -1])
Matrix{T}(_::H) where {T} = Matrix{T}([1 1; 1 -1] ./ sqrt(2))
Matrix{T}(_::S) where {T} = Matrix{T}([1 0; 0 1im])
Matrix{T}(_::Sd) where {T} = Matrix{T}([1 0; 0 -1im])
Matrix{F}(_::T) where {F} = Matrix{F}([1 0; 0 cispi(1 // 4)])
Matrix{F}(_::Td) where {F} = Matrix{F}([1 0; 0 cispi(-1 // 4)])

Matrix{T}(g::Rx) where {T} =
    Matrix{T}([cos(g[:θ] / 2) -1im*sin(g[:θ] / 2); -1im*sin(g[:θ] / 2) cos(g[:θ] / 2)])
Matrix{T}(g::Ry) where {T} =
    Matrix{T}([cos(g[:θ] / 2) -sin(g[:θ] / 2); sin(g[:θ] / 2) cos(g[:θ] / 2)])
Matrix{T}(g::Rz) where {T} = Matrix{T}([1 0; 0 cis(g[:θ])])

Matrix{T}(g::U2) where {T} = 1 / sqrt(2) * Matrix{T}([1 -cis(g[:λ]); cis(g[:ϕ]) cis(g[:ϕ] + g[:λ])])
Matrix{T}(g::U3) where {T} = Matrix{T}(
    [
        cos(g[:θ] / 2) -cis(g[:λ])*sin(g[:θ] / 2)
        cis(g[:ϕ])*sin(g[:θ] / 2) cis(g[:ϕ] + g[:λ])*cos(g[:θ] / 2)
    ],
)

# TODO Matrix{T}(g::Control) = ...

Matrix{T}(g::Swap) = Matrix{T}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])

# diagonal matrices
# NOTE efficient multiplication due to no memory swap needed: plain element-wise multiplication
Diagonal(x::AbstractGate) = Diagonal{ComplexF32}(x)
Diagonal{T}(_::AbstractGate) where {T} = error("Implementation not found")

Diagonal{T}(_::I) where {T} = Diagonal{T}(LinearAlgebra.I, 2)
Diagonal{T}(_::Z) where {T} = Diagonal{T}([1, -1])
Diagonal{T}(_::S) where {T} = Diagonal{T}([1, 1im])
Diagonal{T}(_::Sd) where {T} = Diagonal{T}([1, -1im])
Diagonal{T}(_::T) where {T} = Diagonal{T}([1, cispi(1 // 4)])
Diagonal{T}(_::Td) where {T} = Diagonal{T}([1, cispi(-1 // 4)])
Diagonal{T}(_::Rz) where {T} = Diagonal{T}([1 0; 0 cis(g[:θ])])

# permutational matrices (diagonal + permutation)
# Permutation(x::AbstractGate) = Permutation{ComplexF32}(x)
# Permutation{T}(_::AbstractGate) where {T} = error("Implementation not found")