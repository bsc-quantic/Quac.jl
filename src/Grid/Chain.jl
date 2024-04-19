using Graphs

struct Chain <: Grid{1}
    length::Int

    function Chain(length::Int)
        length > 0 || throw(ArgumentError("length must be positive"))
        new(length)
    end
end

Graphs.vertices(grid::Chain) = (CartesianIndex(i) for i in 1:grid.length)
Graphs.nv(grid::Chain) = grid.length

function Graphs.has_vertex(grid::Chain, v::CartesianIndex{1})
    i = only(Tuple(v))
    1 <= i <= grid.length
end

Graphs.edges(grid::Chain) = ((CartesianIndex(i), CartesianIndex(i + 1)) for i in 1:grid.length-1)
Graphs.ne(grid::Chain) = grid.length - 1

function Graphs.has_edge(grid::Chain, e::Tuple{CartesianIndex{1},CartesianIndex{1}})
    i, j = minmax(e...)
    1 <= i <= grid.length && 1 <= j <= grid.length - 1
end

function Base.getindex(grid::Chain, keys::CartesianIndex{1})
    has_vertex(grid, keys) || throw(ArgumentError("invalid vertex index"))
    i
end
