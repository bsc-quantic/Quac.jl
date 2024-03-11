module QuacSimpleWeightedGraphsExt

using Quac
using SimpleWeightedGraphs

function Quac.connectivity(f, circuit::Circuit)
    connections = Iterators.map(Iterators.filter(f, circuit)) do gate
        n = length(lanes(gate))
        if n == 1
            src = dst = only(lanes(gate))
            return [[src, dst, 1]]
        else
            return vcat.(combinations(lanes(gate), 2), 1)
        end
    end |> Iterators.flatten |> collect

    src = [conn[1] for conn in connections]
    dst = [conn[2] for conn in connections]
    weights = [conn[3] for conn in connections]

    return SimpleWeightedGraph(src, dst, weights; combine = +)
end

end
