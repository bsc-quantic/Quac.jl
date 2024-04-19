using Graphs

struct Rectangular <: Grid{2}
    rows::Int
    cols::Int

    function Rectangular(rows::Int, cols::Int)
        rows > 0 || throw(ArgumentError("rows must be positive"))
        cols > 0 || throw(ArgumentError("cols must be positive"))
        new(rows, cols)
    end
end

function Graphs.vertices(grid::Rectangular)
    Iterators.map(splat(CartesianIndex), Iterators.product(1:grid.rows, 1:grid.cols))
end

function Graphs.has_vertex(grid::Rectangular, v::CartesianIndex{2})
    i, j = Tuple(v)
    1 <= i <= grid.rows && 1 <= j <= grid.cols
end

function hedges(grid::Rectangular)
    ((CartesianIndex(i, j), CartesianIndex(i, j + 1)) for i in 1:grid.rows for j in 1:(grid.cols-1))
end

function vedges(grid::Rectangular)
    ((CartesianIndex(i, j), CartesianIndex(i + 1, j)) for i in 1:(grid.rows-1) for j in 1:grid.cols)
end

Graphs.edges(grid::Rectangular) = Iterators.flatten([hedges(grid), vedges(grid)])
Graphs.has_edge(grid::Rectangular, e::Tuple{CartesianIndex{2},CartesianIndex{2}}) = e ∈ edges(grid)

Graphs.nv(grid::Rectangular) = grid.rows * grid.cols
Graphs.ne(grid::Rectangular) = (grid.rows - 1) * grid.cols + grid.rows * (grid.cols - 1)

Base.getindex(grid::Rectangular, keys::CartesianIndex{2}) = keys[1] + (keys[2] - 1) * grid.rows
