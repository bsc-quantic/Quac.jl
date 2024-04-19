using Graphs

struct Sycamore <: Grid{2} end

Graphs.vertices(::Sycamore) = [
    CartesianIndex(1, 6),
    CartesianIndex(1, 7),
    CartesianIndex(2, 5),
    CartesianIndex(2, 6),
    CartesianIndex(2, 7),
    CartesianIndex(2, 8),
    CartesianIndex(3, 4),
    CartesianIndex(3, 5),
    CartesianIndex(3, 6),
    CartesianIndex(3, 7),
    CartesianIndex(3, 8),
    CartesianIndex(3, 9),
    CartesianIndex(4, 3),
    CartesianIndex(4, 4),
    CartesianIndex(4, 5),
    CartesianIndex(4, 6),
    CartesianIndex(4, 7),
    CartesianIndex(4, 8),
    CartesianIndex(4, 9),
    CartesianIndex(4, 10),
    CartesianIndex(5, 2),
    CartesianIndex(5, 3),
    CartesianIndex(5, 4),
    CartesianIndex(5, 5),
    CartesianIndex(5, 6),
    CartesianIndex(5, 7),
    CartesianIndex(5, 8),
    CartesianIndex(5, 9),
    CartesianIndex(5, 10),
    CartesianIndex(6, 1),
    CartesianIndex(6, 2),
    CartesianIndex(6, 3),
    CartesianIndex(6, 4),
    CartesianIndex(6, 5),
    CartesianIndex(6, 6),
    CartesianIndex(6, 7),
    CartesianIndex(6, 8),
    CartesianIndex(6, 9),
    CartesianIndex(7, 2),
    CartesianIndex(7, 3),
    CartesianIndex(7, 4),
    CartesianIndex(7, 5),
    CartesianIndex(7, 6),
    CartesianIndex(7, 7),
    CartesianIndex(7, 8),
    CartesianIndex(8, 3),
    CartesianIndex(8, 4),
    CartesianIndex(8, 5),
    CartesianIndex(8, 6),
    CartesianIndex(8, 7),
    CartesianIndex(9, 4),
    CartesianIndex(9, 5),
    CartesianIndex(9, 6),
    CartesianIndex(10, 5),
]

Graphs.has_vertex(grid::Sycamore, v::CartesianIndex{2}) = v in vertices(grid)
Graphs.nv(grid::Sycamore) = length(vertices(grid))

function hedges(grid::Sycamore)
    Iterators.filter([(CartesianIndex(i, j), CartesianIndex(i, j + 1)) for i in 1:10 for j in 1:9]) do (u, v)
        has_vertex(grid, u) && has_vertex(grid, v)
    end
end

function vedges(grid::Sycamore)
    Iterators.filter([(CartesianIndex(i, j + 1), CartesianIndex(i, j)) for i in 1:9 for j in 1:10]) do (u, v)
        has_vertex(grid, u) && has_vertex(grid, v)
    end
end

Graphs.edges(::Sycamore) = Iterators.flatten([vedges(Sycamore), hedges(Sycamore)])
Graphs.ne(grid::Sycamore) = length(edges(grid))

Base.getindex(grid::Sycamore, keys::CartesianIndex{2}) = findfirst(==(keys), vertices(grid))