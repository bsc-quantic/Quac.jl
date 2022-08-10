using LinearAlgebra

Matrix(x::AbstractGate) = Matrix{Float32}(x)
Matrix{T}(::AbstractGate) where {T} = error("Implementation not found")

Matrix{T}(::I) where {T} = Matrix{T}(LinearAlgebra.I, 2, 2)
Matrix{T}(::X) where {T} = Matrix{T}([0 1; 1 0])
Matrix{T}(::Y) where {T} = Matrix{T}([0 -1im; 1im 0])
Matrix{T}(::Z) where {T} = Matrix{T}([1 0; 0 -1])
Matrix{T}(::H) where {T} = Matrix{T}([1 1; 1 -1] ./ sqrt(2))
Matrix{T}(::S) where {T} = Matrix{T}([1 0; 0 1im])
Matrix{T}(::Sd) where {T} = Matrix{T}([1 0; 0 -1im])
Matrix{F}(::T) where {F} = Matrix{F}([1 0; 0 cispi(1 // 4)])
Matrix{F}(::Td) where {F} = Matrix{F}([1 0; 0 cispi(-1 // 4)])

Matrix{T}(g::Rx) where {T} =
    Matrix{T}([cos(g.θ / 2) -1im*sin(g.θ / 2); -1im*sin(g.θ / 2) cos(g.θ / 2)])
Matrix{T}(g::Rx) where {T} =
    Matrix{T}([cos(g.θ / 2) -1im*sin(g.θ / 2); -1im*sin(g.θ / 2) cos(g.θ / 2)])
Matrix{T}(g::Rz) where {T} = Matrix{T}([1 0; 0 cis(g.θ)])

Matrix{T}(g::U2) where {T} = 1 / sqrt(2) * Matrix{T}([1 -cis(g.λ); cis(g.ϕ) cis(g.ϕ + g.λ)])
Matrix{T}(g::U3) where {T} = Matrix{T}(
    [
        cos(g.θ / 2) -cis(g.λ)*sin(g.θ / 2)
        cis(g.ϕ)*sin(g.θ / 2) cis(g.ϕ + g.λ)*cos(g.θ / 2)
    ],
)

# TODO Matrix{T}(g::Control) = ...

Matrix{T}(g::Swap) = Matrix{T}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])