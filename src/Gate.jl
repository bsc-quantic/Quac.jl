import Base: adjoint, rand

export Gate
export lanes
export X, Y, Z, H, S, Sd, T, Td
export ParametricGate
export isparametric, parameters
export Rx, Ry, Rz, U1, U2, U3
export Control, Swap
export CX, CY, CZ, CRx, CRy, CRz
export control, target, op

"""
    Gate

Type of gates. Any gate type must fulfill the following requirements:

  - A `lane` field (or property) of type `T <: Union{Int, NTuple{N,Int}} where {N}`
    
      + This requirement can be bypassed by specializing `lanes`.

  - A specialized method for `Base.adjoint(::Type{T})` where `T` is the gate type.
"""
abstract type Gate end

function lanes end

lanes(x::Gate) = (x.lane...,)
Base.adjoint(x::T) where {T<:Gate} = Base.adjoint(T)(lanes(x)...)

"""
    I(lane)

The ``σ_0`` Pauli matrix gate.

# Note

Due to name clashes with `LinearAlgebra.I`, `Quac.I` is not exported by default.
"""
struct I <: Gate
    lane::Int
end

"""
    X(lane)

The ``σ_1`` Pauli matrix gate.
"""
struct X <: Gate
    lane::Int
end

"""
    Y(lane)

The ``σ_2`` Pauli matrix gate.
"""
struct Y <: Gate
    lane::Int
end

"""
    Z(lane)

The ``σ_3`` Pauli matrix gate.
"""
struct Z <: Gate
    lane::Int
end

const Pauli = Union{I,X,Y,Z}

"""
    H(lane)

The Hadamard gate.
"""
struct H <: Gate
    lane::Int
end

for G in [I, X, Y, Z, H]
    @eval Base.adjoint(::Type{$G}) = $G
end

"""
    S(lane)

The ``S`` gate or ``\\frac{π}{2}`` rotation around Z-axis.
"""
struct S <: Gate
    lane::Int
end

Base.adjoint(::Type{S}) = Sd

"""
    Sd(lane)

The ``S^\\dagger`` gate or ``-\\frac{π}{2}`` rotation around Z-axis.
"""
struct Sd <: Gate
    lane::Int
end

Base.adjoint(::Type{Sd}) = S

"""
    T(lane)

The ``T`` gate or ``\\frac{π}{4}`` rotation around Z-axis.
"""
struct T <: Gate
    lane::Int
end

Base.adjoint(::Type{T}) = Td

"""
    Td(lane)

The ``T^\\dagger`` gate or ``-\\frac{π}{4}`` rotation around Z-axis.
"""
struct Td <: Gate
    lane::Int
end

Base.adjoint(::Type{Td}) = T

const Phase = Union{Z,S,Sd,T,Td,Rz}

"""
    ParametricGate

The type of parametric gates.
"""
abstract type ParametricGate <: Gate end

isparametric(::T) where {T<:Gate} = isparametric(T)

isparametric(::Type{<:Gate}) = false
isparametric(::Type{<:ParametricGate}) = true

parameters(::Type{T}) where {T<:ParametricGate} = fieldtype(T, :param).parameters[1]
parameters(x::T) where {T<:ParametricGate} = x.param

Base.propertynames(::Type{T}) where {T<:ParametricGate} = parameters(T)
Base.propertynames(x::T) where {T<:ParametricGate} = parameters(T)

Base.getindex(x::T, key::Symbol) where {T<:ParametricGate} = parameters(x)[key]

Base.adjoint(::Type{T}) where {T<:ParametricGate} = T
Base.adjoint(x::T) where {T<:ParametricGate} =
    Base.adjoint(T)(lanes(x)..., NamedTuple{parameters(T)}(.-(values(parameters(x)))))

Base.rand(::NamedTuple{N,T}) where {N,T} = NamedTuple{N}(rand(type) for type in T.parameters)
Base.rand(::Type{T}, lane::Int) where {T<:ParametricGate} = T(lane, rand(fieldtype(T, :param)))

"""
    Rx(lane, (θ,))

The ``\\theta`` rotation around the X-axis gate.
"""
struct Rx <: ParametricGate
    lane::Int
    param::NamedTuple{(:θ,),Tuple{Float32}}
end

"""
    Ry(lane, (θ,))

The ``\\theta`` rotation around the Y-axis gate.
"""
struct Ry <: ParametricGate
    lane::Int
    param::NamedTuple{(:θ,),Tuple{Float32}}
end

"""
    Rz(lane, (θ,))

The ``\\theta`` rotation around the Z-axis gate.

# Notes

  - The `U1` gate is an alias of `Rz`.
"""
struct Rz <: ParametricGate
    lane::Int
    param::NamedTuple{(:θ,),Tuple{Float32}}
end

const U1 = Rz

"""
    U2(lane, (ϕ, λ))

The ``U2`` gate.
"""
struct U2 <: ParametricGate
    lane::Int
    param::NamedTuple{(:ϕ, :λ),Tuple{Float32,Float32}}
end

"""
    U3(lane, (θ, ϕ, λ))

The ``U3`` gate.
"""
struct U3 <: ParametricGate
    lane::Int
    param::NamedTuple{(:θ, :ϕ, :λ),Tuple{Float32,Float32,Float32}}
end

for G in [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, U2, U3]
    @eval lanes(::Type{$G}) = 1
end

"""
    Control(lane, op::Gate)

A controlled gate.
"""
struct Control{T<:Gate} <: Gate
    lane::Int
    op::T
end
Control{T}(control::Integer, target::Integer) where {T<:Gate} = Control(control, T(target))
Control{T}(args...) where {T} = Control{T}(args[1], T(args[2:end]...))

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

Base.adjoint(::Type{Control{T}}) where {T<:Gate} = Control{adjoint(T)}
Base.adjoint(g::Control{T}) where {T<:Gate} = Control(g.lane, op(g)')

const Toffoli{T} = Control{Control{T}}

isparametric(::Type{Control{T}}) where {T} = isparametric(T)
parameters(g::Control) = parameters(op(g))
parameters(::Type{T}) where {T<:Control} = parameters(op(T))

"""
    Swap(lane1, lane2)

The SWAP gate.
"""
struct Swap <: Gate
    lane::NTuple{2,Int}

    function Swap(a, b)
        new((a, b))
    end
end

Base.adjoint(::Type{Swap}) = Swap
lanes(::Type{Swap}) = 2