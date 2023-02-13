export Gate
export lanes
export X, Y, Z, H, S, Sd, T, Td
export isparametric, parameters
export Rx, Ry, Rz, U1, U2, U3
export Control, Swap
export CX, CY, CZ, CRx, CRy, CRz
export control, target, operator
export Pauli, Phase

abstract type Operator{Params<:NamedTuple} end

# NOTE useful type piracy
Base.keys(::Type{<:NamedTuple{K}}) where {K} = K

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
const Phase = Union{I,Z,S,Sd,T,Td,Rz}

"""
    Gate{Operator}(lanes...; parameters...)

An `Operator` located at some `lanes` and configured with some `parameters`.
"""
struct Gate{Op<:Operator,N,Params}
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
for Op in [:I, :X, :Y, :Z, :H, :S, :Sd, :T, :Td, :U2, :U3, :Rx, :Ry, :Rz, :Swap]
    @eval $Op(lanes...; params...) = Gate{$Op}(lanes...; params...)
end

Control{Op}(lanes...; params...) where {Op} = Gate{Control{Op}}(lanes...; params...)

lanes(g::Gate) = g.lanes
Base.length(::Type{Gate{Op}}) where {Op} = length(Op)
operator(::Type{<:Gate{Op}}) where {Op} = Op
operator(::Gate{Op}) where {Op} = operator(Gate{Op})

function Base.summary(io::IO, gate::Gate)
    flatten = Iterators.flatten
    map = Iterators.map

    Op = operator(gate)
    Args = join(flatten((lanes(gate), map(x -> "$(x[1])=$(x[2])", pairs(parameters(gate))))), ",")
    print(io, "$Op($Args)")
end

Base.show(io::IO, ::MIME"text/plain", gate::Gate) = summary(io, gate)

parameters(g::Gate) = g.parameters
parameters(::Type{<:Gate{Op}}) where {Op} = parameters(Op)
Base.propertynames(::Type{<:Gate{Op}}) where {Op} = (keys(parameters(Op))...,)
Base.propertynames(::G) where {G<:Gate{Op}} where {Op} = propertynames(G)
Base.getproperty(g::Gate{Op}, i::Symbol) where {Op} = i ∈ propertynames(g) ? parameters(g)[i] : getfield(g, i)

Base.adjoint(::Type{<:Gate{Op}}) where {Op} = Gate{adjoint(Op)}
Base.adjoint(::Type{Gate{Op,N,P}}) where {Op,N,P} = Gate{adjoint(Op),N,P}
Base.adjoint(g::Gate{Op}) where {Op} = Gate{Op'}(lanes(g)...; [key => -val for (key, val) in pairs(parameters(g))]...)

# NOTE useful type piracy
Base.rand(::Type{NamedTuple{N,T}}) where {N,T} = NamedTuple{N}(rand(type) for type in T.parameters)

Base.rand(::Type{Op}) where {Op<:Operator} = rand(parameters(Op))
Base.rand(::Type{Gate{Op}}, lanes::Integer...) where {Op} = Gate{Op}(lanes...; rand(Op)...)

# Gate{Control}
targettype(::Type{Op}) where {Op<:Operator} = Op
targettype(::Type{Control{Op}}) where {Op} = Op
targettype(::Type{Control{Op}}) where {Op<:Control} = targettype(Op)
targettype(::Type{<:Gate{Op}}) where {Op} = targettype(Op)

control(g::G) where {G<:Gate{<:Control}} = lanes(g)[1:end-length(targettype(G))]
target(g::G) where {G<:Gate{<:Control}} = lanes(g)[end-length(targettype(G))+1:end]
