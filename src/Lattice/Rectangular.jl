using Graphs

struct Rectangular <: Lattice
    rows::Int
    cols::Int

    function Rectangular(rows::Int, cols::Int)
        rows > 0 || throw(ArgumentError("rows must be positive"))
        cols > 0 || throw(ArgumentError("cols must be positive"))
        new(rows, cols)
    end
end

function Graphs.vertices(lattice::Rectangular)
    Iterators.map(splat(CartesianIndex), Iterators.product(1:lattice.rows, 1:lattice.cols))
end

function Graphs.has_vertex(lattice::Rectangular, v::CartesianIndex{2})
    i, j = Tuple(v)
    1 <= i <= lattice.rows && 1 <= j <= lattice.cols
end

function hedges(lattice::Rectangular)
    ((CartesianIndex(i, j), CartesianIndex(i, j + 1)) for i in 1:lattice.rows for j in 1:(lattice.cols-1))
end

function vedges(lattice::Rectangular)
    ((CartesianIndex(i, j), CartesianIndex(i + 1, j)) for i in 1:(lattice.rows-1) for j in 1:lattice.cols)
end

Graphs.edges(lattice::Rectangular) = Iterators.flatten([hedges(lattice), vedges(lattice)])
Graphs.has_edge(lattice::Rectangular, e::Tuple{CartesianIndex{2},CartesianIndex{2}}) = e ∈ edges(lattice)

Graphs.nv(lattice::Rectangular) = lattice.rows * lattice.cols
Graphs.ne(lattice::Rectangular) = (lattice.rows - 1) * lattice.cols + lattice.rows * (lattice.cols - 1)

Base.getindex(lattice::Rectangular, keys::CartesianIndex{2}) = keys[1] + (keys[2] - 1) * lattice.rows
