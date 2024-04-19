using Graphs

abstract type Grid{N} <: AbstractGraph{CartesianIndex{N}} end

Graphs.edgetype(::Type{<:Grid{N}}) where {N} = NTuple{2,CartesianIndex{N}}

Graphs.is_directed(::Type{<:Grid}) = false
Graphs.is_directed(::Grid) = false

Base.getindex(grid::Grid, keys::Int...) = getindex(grid, CartesianIndex(keys...))

function Graphs.neighbors(grid::Grid, v::CartesianIndex)
    has_vertex(grid, v) || throw(ArgumentError("invalid vertex index"))

    map(Iterators.filter(∋(v), edges(grid))) do e
        v == first(e) ? last(e) : first(e)
    end
end

include("Chain.jl")
include("Rectangular.jl")
include("Honeycomb.jl")
