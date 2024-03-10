export Gate
export lanes
export X, Y, Z, H, S, Sd, T, Td
export isparametric, parameters
export Rx, Ry, Rz, U1, U2, U3, Hz
export Rxx, Ryy, Rzz
export Control, Swap, FSim, SU
export CX, CY, CZ, CRx, CRy, CRz
export control, target, operator
export Pauli, Phase

abstract type Operator end

isparametric(T::Type{<:Operator}) = !isnothing(parameters(T))
isparametric(::T) where {T<:Operator} = isparametric(T)

function parameters end
function lanes end

struct Gate{Op<:Operator,N}
    lanes::NTuple{N,Int}
    operator::Op
end

Base.length(::Type{Gate{Op,N}}) where {Op,N} = N
lanes(g::Gate) = g.lanes

operator(::Type{<:Gate{Op}}) where {Op} = Op
operator(g::Gate) = g.operator

parameters(::Type{<:Gate{Op}}) where {Op} = parameters(Op)
parameters(g::Gate) = parameters(operator(g))

Base.propertynames(::Type{<:Gate{Op}}) where {Op} = (keys(parameters(Op))...,)
Base.propertynames(::G) where {G<:Gate{Op}} where {Op} = propertynames(G)
Base.getproperty(g::Gate{Op}, i::Symbol) where {Op} = i ∈ propertynames(g) ? parameters(g)[i] : getfield(g, i)

Base.adjoint(::Type{Gate{Op,N}}) where {Op,N} = Gate{adjoint(Op),N}
Base.adjoint(g::Gate{Op,N}) where {Op,N} = Gate{Op',N}(g.lanes, g.operator')

function Base.summary(io::IO, gate::Gate)
    Op = operator(gate)
    Args = join(Iterators.flatten((lanes(gate), Iterators.map(x -> "$(x[1])=$(x[2])", pairs(parameters(gate))))), ",")
    print(io, "$Op($Args)")
end

Base.show(io::IO, ::MIME"text/plain", gate::Gate) = summary(io, gate)

macro gatedecl(name, opts...)
    @assert Meta.isidentifier(name)

    opts = collect(opts)
    params = !isempty(opts) && Meta.isexpr(opts[end], :block) ? pop!(opts) : Expr(:block)
    isparametric = !isempty(params.args)

    # options
    n = 1
    adjoint = isparametric ? :parametric : :hermitian

    for opt in opts
        @assert Meta.isexpr(opt, :(=))
        optname, optvalue = opt.args

        if optname == :n
            @assert optvalue isa Integer
            n = optvalue
        elseif optname == :adjoint
            @assert optvalue isa Symbol
            adjoint = optvalue
        else
            throw(ArgumentError("Invalid option $optname"))
        end
    end

    # parameters code
    code_parameters = if isparametric
        local params = Base.remove_linenums!(params)
        paramstuples = map(params.args) do arg
            arg = Meta.isexpr(arg, :(=)) ? arg.args[1] : arg
            @assert Meta.isexpr(arg, :(::))
            tuple(arg.args...)
        end
        paramstuple = :(NamedTuple{($(QuoteNode.(first.(paramstuples))...),),Tuple{$(last.(paramstuples)...)}})
        fieldaccesses = map(field -> :(op.$field), first.(paramstuples))
        quote
            $(esc(:(parameters(::Type{$name}) = $paramstuple)))
            $(esc(:(parameters(op::$name) = $paramstuple($(fieldaccesses...)))))
        end
    else
        quote
            $(esc(:(parameters(::Type{$name}) = nothing)))
            $(esc(:(parameters(::$name) = nothing)))
        end
    end

    # adjoint code
    code_adjoint = if adjoint == :hermitian
        esc(:(Base.adjoint(op::$name) = op))
    elseif adjoint == :parametric
        esc(:(Base.adjoint(op::$name) = parameters(op)))
    else
        esc(:(Base.adjoint(::Type{$name}) = $adjoint()))
    end

    return quote
        $(esc(:(Core.@__doc__ Base.@kwdef struct $name <: $Operator
            $(params.args...)
        end)))

        $(esc(:($name(lanes...; params...) = Gate{$name,$n}(lanes, $name(params...)))))

        $(esc(:(Base.length(::Type{$name}) = $n)))
        $(code_parameters.args...)
        $code_adjoint
    end
end

"""
    I(lane)

The ``σ_0`` Pauli matrix gate.

# Note

Due to name clashes with `LinearAlgebra.I`, `Quac.I` is not exported by default.
"""
@gatedecl I

"""
    X(lane)

The ``σ_1`` Pauli matrix gate.
"""
@gatedecl X

"""
    Y(lane)

The ``σ_2`` Pauli matrix gate.
"""
@gatedecl Y

"""
    Z(lane)

The ``σ_3`` Pauli matrix gate.
"""
@gatedecl Z

"""
    H(lane)

The Hadamard gate.
"""
@gatedecl H

"""
    S(lane)

The ``S`` gate or ``\\frac{π}{2}`` rotation around Z-axis.
"""
@gatedecl S adjoint = Sd

"""
    Sd(lane)

The ``S^\\dagger`` gate or ``-\\frac{π}{2}`` rotation around Z-axis.
"""
@gatedecl Sd adjoint = S

"""
    T(lane)

The ``T`` gate or ``\\frac{π}{4}`` rotation around Z-axis.
"""
@gatedecl T adjoint = Td

"""
    Td(lane)

The ``T^\\dagger`` gate or ``-\\frac{π}{4}`` rotation around Z-axis.
"""
@gatedecl Td adjoint = T

Base.sqrt(::Type{Z}) = S
Base.sqrt(::Type{S}) = T
Base.sqrt(::Type{Sd}) = Td

"""
    Rx(lane, θ)

The ``\\theta`` rotation around the X-axis gate.
"""
@gatedecl Rx begin
    θ::Float64 = 0
end

Base.sqrt(::X) = Rx(θ = π / 2)
Base.sqrt(op::Rx) = Rx(θ = op.θ / 2)

"""
    Rxx(lane1, lane2, θ)

The ``\\theta`` rotation around the XX-axis gate.
"""
@gatedecl Rxx n = 2 begin
    θ::Float64 = 0
end

"""
    Ry(lane, θ)

The ``\\theta`` rotation around the Y-axis gate.
"""
@gatedecl Ry begin
    θ::Float64 = 0
end

Base.sqrt(::Y) = Ry(θ = π / 2)
Base.sqrt(::Ry) = Ry(θ = op.θ / 2)

"""
    Ryy(lane1, lane2, θ)

The ``\\theta`` rotation around the YY-axis gate.
"""
@gatedecl Ryy n = 2 begin
    θ::Float64 = 0
end

"""
    Rz(lane, θ)

The ``\\theta`` rotation around the Z-axis gate.

# Notes

  - The `U1` gate is an alias of `Rz`.
"""
@gatedecl Rz begin
    θ::Float64 = 0
end

Base.sqrt(::Z) = S()
Base.sqrt(::S) = T()
Base.sqrt(::Sd) = Td()
Base.sqrt(op::Rz) = Rz(θ = op.θ / 2)
Base.sqrt(::T) = Rz(θ = π / 8)
Base.sqrt(::Td) = Rz(θ = -π / 8)

"""
    Rzz(lane1, lane2, θ)

The ``\\theta`` rotation around the ZZ-axis gate.
"""
@gatedecl Rzz n = 2 begin
    θ::Float64 = 0
end

const U1 = Rz

"""
    U2(lane, ϕ, λ)

The ``U2`` gate.
"""
@gatedecl U2 begin
    ϕ::Float64 = 0
    λ::Float64 = 0
end

"""
    U3(lane, θ, ϕ, λ)

The ``U3`` gate.
"""
@gatedecl U3 begin
    θ::Float64 = 0
    ϕ::Float64 = 0
    λ::Float64 = 0
end

"""
    Hz(lane, θ, ϕ)

The Hz (PhasedXPow) gate, equivalent to ``Z^ϕ X^θ Z^{-ϕ}``.
"""
@gatedecl Hz begin
    θ::Float64 = 0
    ϕ::Float64 = 0
end

"""
    Swap(lane1, lane2)

The SWAP gate.
"""
@gatedecl Swap n = 2

"""
    FSim(lane1, lane2, θ, ϕ)

The FSim (Fermionic Simulation) gate.
"""
@gatedecl FSim n = 2 begin
    θ::Float64 = 0
    ϕ::Float64 = 0
end

"""
    Control(lane, op::Gate)

A controlled gate.
"""
Base.@kwdef struct Control{Op<:Operator} <: Operator
    target::Op
end

Control{Op}(lanes...; params...) where {Op} = Gate{Control{Op},length(lanes)}(lanes...; params...)
Control(lane, op::Gate{Op,N}) where {Op,N} = Gate{Control{Op},N + 1}(lane, lanes(op)...; parameters(op)...)

Base.length(::Type{Control{T}}) where {T<:Operator} = 1 + length(T)
isparametric(::Type{<:Control{T}}) where {T<:Operator} = isparametric(T)
parameters(::Type{Control{Op}}) where {Op} = parameters(Op)

Base.adjoint(::Type{Control{Op}}) where {Op} = Control{adjoint(Op)}
Base.adjoint(op::Control{Op}) where {Op} =
    if isparametric(op)
        Control{Op'}([key => -val for (key, val) in pairs(parameters(op))]...)
    else
        Control{adjoint(Op)}()
    end

for Op in [:X, :Y, :Z, :Rx, :Ry, :Rz]
    @eval const $(Symbol("C" * String(Op))) = Control{$Op}
end

# Gate{Control}
targettype(::Type{Op}) where {Op<:Operator} = Op
targettype(::Type{Control{Op}}) where {Op} = Op
targettype(::Type{Control{Op}}) where {Op<:Control} = targettype(Op)
targettype(::Type{<:Gate{Op}}) where {Op} = targettype(Op)

control(g::G) where {G<:Gate{<:Control}} = lanes(g)[1:end-length(targettype(G))]
target(g::G) where {G<:Gate{<:Control}} = lanes(g)[end-length(targettype(G))+1:end]

"""
    SU{N}(lane_1, lane_2, ..., lane_N, array)

The `SU{N}` multi-qubit general unitary gate that can be used to represent any unitary matrix that acts on
`N` qubits. A new random `SU{N}` can be created with `rand(SU{N}, lanes...)`, where `N` is the dimension
of the unitary matrix and `lanes` are the qubit lanes on which the gate acts.

# Note

Unlike the general notation, `N` is not the dimension of the unitary matrix, but the number of qubits on which it acts.
The dimension of the unitary matrix is ``2^N \\times 2^N``.
"""
struct SU{N} <: Operator
    matrix::Matrix

    function SU{N}(; matrix) where {N}
        size(matrix) == (2^N, 2^N) || throw(ArgumentError("`matrix` ($(size(matrix))) must be a (2^N,2^N)-size matrix"))
        LinearAlgebra.det(matrix) ≈ 1 || throw(ArgumentError("`matrix` is not unitary"))
        new(matrix)
    end
end

function SU{N}(lanes...; params...) where {N}
    length(lanes) == N || throw(ArgumentError("SU{$N} requires $N lanes"))
    Gate{SU{N},N}(lanes, SU{N}(; params...))
end

Base.length(::Type{SU{N}}) where {N} = N
isparametric(::Type{<:SU}) = true
parameters(::Type{SU{N}}) where {N} = NamedTuple{(:matrix),Tuple{Matrix}}

Base.adjoint(::Type{SU{N}}) where {N} = SU{N}
Base.adjoint(op::SU{N}) where {N} = SU{N}(; matrix = op.matrix')

# operator sets
const Pauli = Union{I,X,Y,Z}
const Phase = Union{I,Z,S,Sd,T,Td,Rz}

randtuple(::Type{NamedTuple{N,T}}) where {N,T} = NamedTuple{N}(rand(type) for type in T.parameters)

Base.rand(::Type{Op}) where {Op<:Operator} = randtuple(parameters(Op))
Base.rand(::Type{Gate{Op}}, lanes::Integer...) where {Op} = Gate{Op}(lanes...; rand(Op)...)
