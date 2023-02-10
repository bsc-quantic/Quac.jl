import Base: adjoint, rand
using Base: front, tail

export Gate
export lanes
export X, Y, Z, H, S, Sd, T, Td
export isparametric, parameters
export Rx, Ry, Rz, U1, U2, U3
export Control, Swap
export CX, CY, CZ, CRx, CRy, CRz
export control, target, operator

abstract type Operator{Params<:NamedTuple} end

# `Operator` with no parameters
const StaticOperator = Operator{NamedTuple{(),Tuple{}}}

parameters(::Type{<:Operator{Params}}) where {Params} = Params
isparametric(::Type{T}) where {T<:Operator} = parameters(T) !== NamedTuple{(),Tuple{}}

for Op in [:I, :X, :Y, :Z, :H, :S, :Sd, :T, :Td]
    @eval abstract type $Op <: StaticOperator end
    @eval Base.length(::Type{$Op}) = 1
end

abstract type Rx <: Operator{NamedTuple{(:θ,),Tuple{Float64}}} end
abstract type Ry <: Operator{NamedTuple{(:θ,),Tuple{Float64}}} end
abstract type Rz <: Operator{NamedTuple{(:θ,),Tuple{Float64}}} end

for Op in [:Rx, :Ry, :Rz]
    @eval Base.length(::Type{$Op}) = 1
end

const U1 = Rz

abstract type U2 <: Operator{NamedTuple{(:ϕ, :λ),Tuple{Float64,Float64}}} end
Base.length(::Type{U2}) = 1

abstract type U3 <: Operator{NamedTuple{(:θ, :ϕ, :λ),Tuple{Float64,Float64,Float64}}} end
Base.length(::Type{U3}) = 1

abstract type Swap <: StaticOperator end
Base.length(::Type{Swap}) = 2

abstract type Control{Op<:Operator} <: Operator{NamedTuple{(:target,),Tuple{Operator}}} end
Base.length(::Type{Control{T}}) where {T<:Operator} = 1 + length(T)
parameters(::Type{Control{Op}}) where {Op} = parameters(Op)

for Op in [:X, :Y, :Z, :Rx, :Ry, :Rz]
    @eval const $(Symbol("C" * String(Op))) = Control{$Op}
end

# adjoints
for Op in [:I, :X, :Y, :Z, :Rx, :Ry, :Rz, :H, :Swap]
    @eval Base.adjoint(::Type{$Op}) = $Op
end

Base.adjoint(::Type{S}) = Sd
Base.adjoint(::Type{Sd}) = S
Base.adjoint(::Type{T}) = Td
Base.adjoint(::Type{Td}) = T

Base.adjoint(::Type{Control{Op}}) where {Op} = Control{adjoint(Op)}

# operator sets
const Pauli = Union{I,X,Y,Z}
const Phase = Union{Z,S,Sd,T,Td,Rz}

"""
    Gate{Operator}(lanes...; parameters...)

An `Operator` located at some `lanes` and configured with some `parameters`.
"""
struct Gate{Op,N,Params}
    lanes::NTuple{N,Int}
    parameters::Params

    function Gate{Op}(lanes...; params...) where {Op<:Operator}
        N = length(Op)
        P = parameters(Op)
        params = NamedTuple{tuple(keys(params)...),Tuple{typeof.(collect(values(params)))...}}((values(params)...,))
        new{Op,N,P}(lanes, params)
    end
end

# constructor aliases
for Op in filter(x -> x isa DataType, subtypes(Operator))
    @eval $(Symbol(Op))(lanes...; params...) = Gate{$Op}(lanes...; params...)
end

# TODO Gate{Control} constructor

lanes(g::Gate) = g.lanes
Base.length(::Type{Gate{Op}}) where {Op} = length(Op)
operator(::Gate{Op}) where {Op} = Op

parameters(g::Gate) = g.parameters
parameters(::Type{Gate{Op}}) where {Op} = parameters(Op)
Base.propertynames(::Type{Gate{Op}}) where {Op} = (first(parameters(Op).parameters)...,)
Base.propertynames(::G) where {G<:Gate{Op}} where {Op} = propertynames(G)
Base.getproperty(g::Gate{Op}, i::Symbol) where {Op} = i ∈ propertynames(g) ? parameters(g)[i] : getfield(g, i)

Base.adjoint(g::Gate{Op}) where {Op} = Gate{Op'}(lanes(g)...; [key => -val for (key, val) in pairs(parameters(g))]...)

# NOTE useful type piracy
Base.rand(::Type{NamedTuple{N,T}}) where {N,T} = NamedTuple{N}(rand(type) for type in T.parameters)

Base.rand(::Type{Op}) where {Op<:Operator} = rand(parameters(Op))
Base.rand(::Type{Gate{Op}}, lanes::Integer...) where {Op} = Gate{Op}(lanes...; rand(Op)...)

# Gate{Control}
op(::Type{Op}) where {Op<:Operator} = Op
op(::Type{Control{Op}}) where {Op} = Op
op(::Type{Control{Op}}) where {Op<:Control} = op(Op)
op(::Type{<:Gate{Op}}) where {Op} = op(Op)

control(g::G) where {G<:Gate{<:Control}} = lanes(g)[1:end-length(op(G))]
target(g::G) where {G<:Gate{<:Control}} = lanes(g)[end-length(op(G))+1:end]
