using Graphs

struct Honeycomb <: Grid{2}
    rows::Int
    cols::Int

    function Honeycomb(rows::Int, cols::Int)
        rows > 0 || throw(ArgumentError("rows must be positive"))
        cols > 0 || throw(ArgumentError("cols must be positive"))
        new(rows, cols)
    end
end

Graphs.vertices(grid::Honeycomb) = Iterators.flatten([
    (CartesianIndex(1, j) for j in 1:(2*grid.cols+1)), # first row
    (CartesianIndex(i, j) for i in 2:grid.rows for j in 1:(2*grid.cols+2)), # middle rows
    (CartesianIndex(grid.rows + 1, j) for j in 2:(2*grid.cols+2)), # last row
])

function Graphs.has_vertex(grid::Honeycomb, v::CartesianIndex{2})
    i, j = Tuple(v)
    1 <= i <= grid.rows + 1 && 1 <= j <= 2 * grid.cols + 2 || return false

    if i == 1
        return j < 2 * grid.cols + 2 ? true : false
    elseif i == grid.rows + 1
        return j > 1 ? true : false
    else
        return true
    end
end

Graphs.edges(grid::Honeycomb) = Iterators.flatten([aedges(grid), dedges(grid), vedges(grid)])

"""
    vedges(grid::Honeycomb)

Return an iterator of the vertical edges of `grid`.
"""
vedges(grid::Honeycomb) =
    ((CartesianIndex(i, j), CartesianIndex(i + 1, j)) for i in 1:grid.rows for j in 1:(2*grid.cols+2) if i % 2 == j % 2)

# TODO test with other grid configurations
"""
    aedges(grid::Honeycomb)

Return an iterator of the antidiagonal edges of `grid`.
"""
aedges(grid::Honeycomb) = Iterators.flatten([
    ((CartesianIndex(1, j), CartesianIndex(1, j + 1)) for j in 1:2:grid.cols*2), # first row
    ((CartesianIndex(i, j), CartesianIndex(i, j + 1)) for i in 2:grid.rows for j in 1+(i+1)%2:2:grid.cols*2+1), # middle rows
    (
        (CartesianIndex(grid.rows + 1, j), CartesianIndex(grid.rows + 1, j + 1)) for
        j in 2+(grid.rows+1)%2:2:grid.cols*2+2
    ), # last row
])

# TODO test with other grid configurations
"""
    dedges(grid::Honeycomb)

Return an iterator of the diagonal edges of `grid`.
"""
dedges(grid::Honeycomb) = Iterators.flatten([
    ((CartesianIndex(1, j), CartesianIndex(1, j + 1)) for j in 2:2:grid.cols*2), # first row
    ((CartesianIndex(i, j), CartesianIndex(i, j + 1)) for i in 2:grid.rows for j in 1+i%2:2:grid.cols*2+1), # middle rows
    ((CartesianIndex(grid.rows + 1, j), CartesianIndex(grid.rows + 1, j + 1)) for j in 2+grid.rows%2:2:grid.cols*2+1), # last row
])

Graphs.has_edge(grid::Honeycomb, e::Tuple{CartesianIndex{2},CartesianIndex{2}}) = e ∈ edges(grid)

Graphs.nv(grid::Honeycomb) = (grid.rows + 1) * (2 * grid.cols + 2) - 2
function Graphs.ne(grid::Honeycomb)
    nev = (grid.cols + 1) * grid.rows
    neh = (2 * grid.cols + 1) * (grid.rows - 1) + 2 * 2 * grid.cols
    nev + neh
end

function Base.getindex(grid::Honeycomb, key::CartesianIndex{2})
    has_vertex(grid, key) || throw(ArgumentError("invalid vertex index"))
    i, j = Tuple(key)

    j == 1 && return i
    if j == grid.cols * 2 + 2
        i -= 1
    end

    m = grid.rows + 1
    offset = grid.rows + (grid.rows % 2)

    m * (j - 2) + i + offset
end
