import Base: show, adjoint, rand

export AbstractGate
export lane
export I, X, Y, Z, H, S, Sd, T, Td
export AbstractParametricGate
export isparametric, parameters
export Rx, Ry, Rz, U1, U2, U3
export Control, Swap
export CX, CY, CZ, CRx, CRy, CRz

"""
Expected requirements:
- A field called `lane`
"""
abstract type AbstractGate end

lanes(x::AbstractGate) = (x.lane,)
Base.adjoint(x::T) where {T<:AbstractGate} = Base.adjoint(T)(lanes(x))

"""
    I(lane)

The ``σ₀`` Pauli matrix gate.
"""
struct I <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{I}) = I

"""
    X(lane)

The ``σ₁`` Pauli matrix gate.
"""
struct X <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{X}) = X

"""
    Y(lane)

The ``σ₂`` Pauli matrix gate.
"""
struct Y <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{Y}) = Y

"""
    Z(lane)

The ``σ₂`` Pauli matrix gate.
"""
struct Z <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{Z}) = Z

"""
    H(lane)

The Hadamard gate.
"""
struct H <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{H}) = H

"""
    S(lane)

The S gate or ``π/2`` rotation around Z-axis.
"""
struct S <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{S}) = Sd

"""
    Sd(lane)

The S† gate or ``-π/2`` rotation around Z-axis.
"""
struct Sd <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{Sd}) = S

"""
    T(lane)

The T gate or ``π/4`` rotation around Z-axis.
"""
struct T <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{T}) = Td

"""
    Td(lane)

The T† gate or ``-π/4`` rotation around Z-axis.
"""
struct Td <: AbstractGate
    lane::Int
end

Base.adjoint(::Type{Td}) = T

abstract type AbstractParametricGate <: AbstractGate end

isparametric(::Type{<:AbstractGate}) = false
isparametric(::Type{<:AbstractParametricGate}) = true

parameters(::Type{T}) where {T<:AbstractParametricGate} = keys(fieldtype(T, :param))
parameters(x::T) where {T<:AbstractParametricGate} = x.param

Base.getindex(x::T, key::Symbol) where {T<:AbstractParametricGate} = x.param[key]

Base.adjoint(::Type{T}) where {T<:AbstractParametricGate} = T
Base.adjoint(x::T) where {T<:AbstractParametricGate} = Base.adjoint(T)(lanes())

Base.rand(::NamedTuple{N,T}) where {N,T} = NamedTuple{N}(rand(type) for type in T.parameters)
Base.rand(::Type{T}, lane::Int) where {T<:AbstractParametricGate} = T(lane, rand(fieldtype(T, :param)))

"""
    Rx(lane, (θ,))

The ``\\theta`` rotation around the X-axis gate.
"""
struct Rx <: AbstractParametricGate
    lane::Int
    param::NamedTuple{(:θ,),Tuple{Float32}}
end

"""
    Ry(lane, (θ,))

The ``\\theta`` rotation around the Y-axis gate.
"""
struct Ry <: AbstractParametricGate
    lane::Int
    param::NamedTuple{(:θ,),Tuple{Float32}}
end

"""
    Rz(lane, (θ,))

The ``\\theta`` rotation around the Z-axis gate.

# Notes
- The `U1` gate is an alias of `Rz`.
"""
struct Rz <: AbstractParametricGate
    lane::Int
    param::NamedTuple{(:θ,),Tuple{Float32}}
end

U1 = Rz

struct U2 <: AbstractParametricGate
    lane::Int
    param::NamedTuple{(:ϕ, :λ),Tuple{Float32,Float32}}
end

struct U3 <: AbstractParametricGate
    lane::Int
    param::NamedTuple{(:θ, :ϕ, :λ),Tuple{Float32,Float32,Float32}}
end

struct Control{T} <: AbstractGate
    lane::Int
    op::T
end

Control{T}(control::Integer, target::Integer) where {T<:AbstractGate} = Control(control, T(target))
Control{T}(lanes::Integer...) where {T<:AbstractGate} = Control(first(lanes), Control{T}(Iterators.drop(lanes, 1)...))

CX(control, target) = Control(control, X(target))
CY(control, target) = Control(control, Y(target))
CZ(control, target) = Control(control, Z(target))
CRx(control, target, θ) = Control(control, Rx(target, θ))
CRy(control, target, θ) = Control(control, Ry(target, θ))
CRz(control, target, θ) = Control(control, Rz(target, θ))

control(g::Control{T}) where {T} = g.lane
control(g::Control{T}) where {T<:Control} = (g.lane, control(g.op)...)
target(g::Control{T}) where {T} = lanes(g.op)
target(g::Control{T}) where {T<:Control} = target(g.op)
lanes(g::Control{T}) where {T} = (control(g)..., target(g)...)

Base.adjoint(::Type{Control{T}}) where {T<:AbstractGate} = Control{adjoint(T)}
Base.adjoint(g::Control{T}) where {T<:AbstractGate} = Control(lanes(g), adjoint(g.op))

const Toffoli{T} = Control{Control{T}}

# special case for Control{T} where {T<:AbstractParametricGate}, as it is parametric
isparametric(::Type{Control{<:AbstractParametricGate}}) = true
isparametric(::Type{Control{T}}) where {T<:Control} = isparametric(T)
parameters(g::Control{<:AbstractParametricGate}) = parameters(g.op)
parameters(g::Control{<:Control}) = parameters(g.op)

"""
    Swap(lanes)

The SWAP gate.
"""
struct Swap <: AbstractGate
    lane::NTuple{2,Int}
end

Base.adjoint(::Type{Swap}) = Swap