using Graphs

struct Honeycomb <: Lattice
    rows::Int
    cols::Int

    function Honeycomb(rows::Int, cols::Int)
        rows > 0 || throw(ArgumentError("rows must be positive"))
        cols > 0 || throw(ArgumentError("cols must be positive"))
        new(rows, cols)
    end
end

Graphs.vertices(lattice::Honeycomb) = Iterators.flatten([
    (CartesianIndex(1, j) for j in 1:(2*lattice.cols+1)), # first row
    (CartesianIndex(i, j) for i in 2:lattice.rows for j in 1:(2*lattice.cols+2)), # middle rows
    (CartesianIndex(lattice.rows + 1, j) for j in 2:(2*lattice.cols+2)), # last row
])

function Graphs.has_vertex(lattice::Honeycomb, v::CartesianIndex{2})
    i, j = Tuple(v)
    1 <= i <= lattice.rows + 1 && 1 <= j <= 2 * lattice.cols + 2 || return false

    if i == 1
        return j < 2 * lattice.cols + 2 ? true : false
    elseif i == lattice.rows + 1
        return j > 1 ? true : false
    else
        return true
    end
end

Graphs.edges(lattice::Honeycomb) = Iterators.flatten([aedges(lattice), dedges(lattice), vedges(lattice)])

"""
    vedges(lattice::Honeycomb)

Return an iterator of the vertical edges of `lattice`.
"""
vedges(lattice::Honeycomb) = (
    (CartesianIndex(i, j), CartesianIndex(i + 1, j)) for i in 1:lattice.rows for
    j in 1:(2*lattice.cols+2) if i % 2 == j % 2
)

# TODO test with other lattice configurations
"""
    aedges(lattice::Honeycomb)

Return an iterator of the antidiagonal edges of `lattice`.
"""
aedges(lattice::Honeycomb) = Iterators.flatten([
    ((CartesianIndex(1, j), CartesianIndex(1, j + 1)) for j in 1:2:lattice.cols*2), # first row
    ((CartesianIndex(i, j), CartesianIndex(i, j + 1)) for i in 2:lattice.rows for j in 1+(i+1)%2:2:lattice.cols*2+1), # middle rows
    (
        (CartesianIndex(lattice.rows + 1, j), CartesianIndex(lattice.rows + 1, j + 1)) for
        j in 2+(lattice.rows+1)%2:2:lattice.cols*2+2
    ), # last row
])

# TODO test with other lattice configurations
"""
    dedges(lattice::Honeycomb)

Return an iterator of the diagonal edges of `lattice`.
"""
dedges(lattice::Honeycomb) = Iterators.flatten([
    ((CartesianIndex(1, j), CartesianIndex(1, j + 1)) for j in 2:2:lattice.cols*2), # first row
    ((CartesianIndex(i, j), CartesianIndex(i, j + 1)) for i in 2:lattice.rows for j in 1+i%2:2:lattice.cols*2+1), # middle rows
    (
        (CartesianIndex(lattice.rows + 1, j), CartesianIndex(lattice.rows + 1, j + 1)) for
        j in 2+lattice.rows%2:2:lattice.cols*2+1
    ), # last row
])

Graphs.has_edge(lattice::Honeycomb, e::Tuple{CartesianIndex{2},CartesianIndex{2}}) = e ∈ edges(lattice)

Graphs.nv(lattice::Honeycomb) = (lattice.rows + 1) * (2 * lattice.cols + 2) - 2
function Graphs.ne(lattice::Honeycomb)
    nev = (lattice.cols + 1) * lattice.rows
    neh = (2 * lattice.cols + 1) * (lattice.rows - 1) + 2 * 2 * lattice.cols
    nev + neh
end

function Base.getindex(lattice::Honeycomb, key::CartesianIndex{2})
    has_vertex(lattice, key) || throw(ArgumentError("invalid vertex index"))
    i, j = Tuple(key)

    j == 1 && return i
    if j == lattice.cols * 2 + 2
        i -= 1
    end

    m = lattice.rows + 1
    offset = lattice.rows + (lattice.rows % 2)

    m * (j - 2) + i + offset
end
