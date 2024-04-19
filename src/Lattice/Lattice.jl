using Graphs

abstract type Lattice <: AbstractGraph{CartesianIndex{2}} end

Graphs.edgetype(::Type{<:Lattice}) = NTuple{2,CartesianIndex{2}}

Graphs.is_directed(::Type{<:Lattice}) = false
Graphs.is_directed(::Lattice) = false

Base.getindex(lattice::Lattice, keys::Int...) = getindex(lattice, CartesianIndex(keys...))

function Graphs.neighbors(lattice::Lattice, v::CartesianIndex{2})
    has_vertex(lattice, v) || throw(ArgumentError("invalid vertex index"))

    map(Iterators.filter(∋(v), edges(lattice))) do e
        v == first(e) ? last(e) : first(e)
    end
end

include("Rectangular.jl")
include("Honeycomb.jl")
