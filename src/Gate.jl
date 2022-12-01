import Base: show, adjoint, rand

export AbstractGate
export lanes
export I, X, Y, Z, H, S, Sd, T, Td
export AbstractParametricGate
export isparametric, parameters
export Rx, Ry, Rz, U1, U2, U3
export Control, Swap
export CX, CY, CZ, CRx, CRy, CRz
export control, target, op

"""
    AbstractGate

Type of gates. Any gate type must fulfill the following requirements:

  - A `lane` field (or property) of type `T <: Union{Int, NTuple{N,Int}} where {N}`
    
      + This requirement can be bypassed by specializing `lanes`.

  - A specialized method for `Base.adjoint(::Type{T})` where `T` is the gate type.
"""
abstract type AbstractGate end

function lanes end

lanes(x::AbstractGate) = (x.lane,)
Base.adjoint(x::T) where {T<:AbstractGate} = Base.adjoint(T)(lanes(x)...)

"""
    I(lane)

The ``σ₀`` Pauli matrix gate.
"""
struct I <: AbstractGate
    lane::Int
end

"""
    X(lane)

The ``σ₁`` Pauli matrix gate.
"""
struct X <: AbstractGate
    lane::Int
end

"""
    Y(lane)

The ``σ₂`` Pauli matrix gate.
"""
struct Y <: AbstractGate
    lane::Int
end

"""
    Z(lane)

The ``σ₂`` Pauli matrix gate.
"""
struct Z <: AbstractGate
    lane::Int
end

"""
    H(lane)

The Hadamard gate.
"""
struct H <: AbstractGate
    lane::Int
end

for G in [I, X, Y, Z, H]
    @eval Base.adjoint(::Type{$G}) = $G
end

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

"""
    AbstractParametricGate

The type of parametric gates.
"""
abstract type AbstractParametricGate <: AbstractGate end

isparametric(::T) where {T<:AbstractGate} = isparametric(T)

isparametric(::Type{<:AbstractGate}) = false
isparametric(::Type{<:AbstractParametricGate}) = true

parameters(::Type{T}) where {T<:AbstractParametricGate} = fieldtype(T, :param).parameters[1]
parameters(x::T) where {T<:AbstractParametricGate} = x.param

Base.propertynames(::Type{T}) where {T<:AbstractParametricGate} = parameters(T)
Base.propertynames(x::T) where {T<:AbstractParametricGate} = parameters(T)

Base.getindex(x::T, key::Symbol) where {T<:AbstractParametricGate} = parameters(x)[key]

Base.adjoint(::Type{T}) where {T<:AbstractParametricGate} = T
Base.adjoint(x::T) where {T<:AbstractParametricGate} =
    Base.adjoint(T)(lanes(x)..., NamedTuple{parameters(T)}(.-(values(parameters(x)))))

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

"""
    U2

The ``U2`` gate.
"""
struct U2 <: AbstractParametricGate
    lane::Int
    param::NamedTuple{(:ϕ, :λ),Tuple{Float32,Float32}}
end

"""
    U3

The ``U3`` gate.
"""
struct U3 <: AbstractParametricGate
    lane::Int
    param::NamedTuple{(:θ, :ϕ, :λ),Tuple{Float32,Float32,Float32}}
end

for G in [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, U2, U3]
    @eval lanes(::Type{$G}) = 1
end

"""
    Control

A controlled gate.
"""
struct Control{T<:AbstractGate} <: AbstractGate
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

control(g::Control{T}) where {T} = (g.lane,)
control(g::Control{T}) where {T<:Control} = (g.lane, control(g.op)...)
target(g::Control{T}) where {T} = lanes(g.op)
target(g::Control{T}) where {T<:Control} = target(g.op)
lanes(g::Control{T}) where {T} = (control(g)..., target(g)...)
lanes(::Type{Control{T}}) where {T} = 1 + lanes(T)

op(g::Control{T}) where {T} = g.op
op(g::Control{T}) where {T<:Control} = op(g.op)
op(::Type{Control{T}}) where {T} = T
op(::Type{Control{T}}) where {T<:Control} = op(T)

Base.adjoint(::Type{Control{T}}) where {T<:AbstractGate} = Control{adjoint(T)}
Base.adjoint(g::Control{T}) where {T<:AbstractGate} = Control(g.lane, op(g)')

const Toffoli{T} = Control{Control{T}}

isparametric(::Type{Control{T}}) where {T} = isparametric(T)
parameters(g::Control) = parameters(op(g))
parameters(::Type{T}) where {T<:Control} = parameters(op(T))

"""
    Swap(lanes)

The SWAP gate.
"""
struct Swap <: AbstractGate
    lane::NTuple{2,Int}
end

Base.adjoint(::Type{Swap}) = Swap
lanes(::Type{Swap}) = 2