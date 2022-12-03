import Base: push!, length, iterate, IteratorSize, in, adjoint, getindex, show
using Base.Iterators: enumerate, filter
using SimpleWeightedGraphs
using Combinatorics

export Circuit
export lanes, connectivity, moments

"""
    Element{T}

One element of a queue, which contains an element of type `T` and multiple priority numbers.
"""
struct Element{T}
    data::T
    priority::Vector{Pair{Int,Int}}
end

Element(gate::AbstractGate, priority) = Element{AbstractGate}(gate, priority)

data(e::Element) = e.data

in(e::Element, head::Vector{Int}) = all(p == head[lane] for (lane, p) in e.priority)

"""
A quantum circuit implementation using multi-priority queues.

  - Queues are gate-buffers in qubit lanes.
  - Multi-priority numbers can be retrieved procedurally from gate-lanes encoded inside the gates of the queues.
"""
struct Circuit
    lanes::Vector{Vector{Element{AbstractGate}}}

    Circuit(n::Int) = new(fill(Element[], n))
    Circuit(lanes::Vector{Vector{Element{AbstractGate}}}) = new(lanes)
end

lanes(circ::Circuit) = length(circ.lanes)

Base.length(circ::Circuit) = sum(length(lane) for lane in circ.lanes)

Base.isempty(circ::Circuit) = all(isempty, circ.lanes)

"""
    push!(circ, gate)

Appends a gate to the circuit.
"""
Base.push!(circ::Circuit, gate::AbstractGate) = begin
    new_priority = [lane => length(circ.lanes[lane]) + 1 for lane in lanes(gate)]
    el = Element{AbstractGate}(gate, new_priority)

    for lane in lanes(gate)
        queue = circ.lanes[lane]
        push!(queue, el)
    end
end

Base.push!(circ::Circuit, gates::AbstractGate...) = foreach(g -> push!(circ, g), gates)

"""
    Base.iterate(circ::Circuit[, state])

Retrieves next gate from `state` by travelling through a topologically sorted path.

# Arguments

  - `circ::Circuit`
  - `state` (or head) should be a `NTuple{N, Int}` where `N` is the number of lanes. Each element is a pointer to the next gate on each lane.
"""
Base.iterate(circ::Circuit, state = fill(1, lanes(circ))) = begin
    # retrieve gates on the edge of the cut
    candidates =
        enumerate(state) |>
        (x -> filter(y -> ((lane, head) = y; head <= length(circ.lanes[lane])), x)) .|>
        (x -> begin
            lane, i = x
            circ.lanes[lane][i]
        end)
    if isempty(candidates)
        return nothing
    end

    # choose first valid gate
    winner = filter(x -> x ∈ state, candidates) |> first

    # update head by advancing cut on involved lanes
    for (lane, priority) in winner.priority
        state[lane] = priority + 1
    end

    (data(winner), state)
end

Base.IteratorSize(::Type{Circuit}) = Base.HasLength()

"""
    Base.adjoint(circ)

Retrieve the adjoint circuit which fulfills the following equality.

```julia
circ * circ' == I(n)
```

# Notes

If all gates are hermitian, then the following equality also holds.

```julia
circ * circ' == circ' * circ == I(n)
```
"""
Base.adjoint(circ::Circuit) = begin
    lanes = [
        [
            Element{AbstractGate}(
                data(el),
                [laneid => length(circ.lanes[laneid]) - p + 1 for (laneid, p) in el.priority],
            ) for el in lane
        ] for lane in circ.lanes
    ]

    Circuit(lanes)
end

"""
    Base.getindex(circ, lane, index)

Retrieve gate at lane `lane` and depth `index`.
"""
Base.getindex(circ::Circuit, lane, index) = data(circ.lanes[lane][index])

Base.show(io::IO, circ::Circuit) = print(io, "Circuit(#lanes=$(lanes(circ)), #gates=$(length(circ)))")

"""
    connectivity([f,] circuit)

Generate connectivity graph between qubits.

# Arguments

  - f: Function to filter gates from circuit.
  - circuit: Circuit.
"""
function connectivity(f, circ::Circuit)
    connections = Iterators.map(Iterators.filter(f, circ)) do gate
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

connectivity(circ::Circuit) = connectivity(gate -> length(lanes(gate)) >= 2, circ)

"""
    moments(circuit)

Return moments (lists of gates that _can_ execute at the same time) of the circuit.
"""
function moments(circ::Circuit)
    m = [AbstractGate[]]

    for gate in circ
        if isempty(last(m)) || isdisjoint(lanes(gate), ∪(lanes.(last(m))...))
            push!(last(m), gate)
        else
            push!(m, [gate])
        end
    end

    return m
end

"""
    hcat(circ::Circuit...)

Join circuits in the temporal dimension.
"""
function Base.hcat(circs::Circuit...)
    !allequal(lanes(circ) for circ in circs) && throw(DimensionMismatch("circuits must have same lanes"))

    reduce(map(x -> x.lanes, circs)) do acc, circ
        map(zip(acc, circ)) do (acclane, lane)
            offset = length(acclane)

            return [acclane..., map(lane) do el
                priority = [k => v + offset for (k, v) in el.priority]
                return Element{AbstractGate}(data(el), priority)
            end...]
        end
    end |> Circuit
end

"""
    vcat(circ::Circuit...)

Join circuits in the spatial dimension.
"""
function Base.vcat(circs::Circuit...)
    offsets = [0, cumsum(lanes.(circs))[1:end-1]...]

    mapreduce(vcat, zip(offsets, circs)) do (offset, circ)
        map(circ.lanes) do lane
            map(enumerate(lane)) do (i, el)
                priority = [k + offset => v for (k, v) in el.priority]
                gate = data(el)

                if isparametric(gate)
                    gate = typeof(gate)((lanes(gate) .+ offset)..., parameters(gate))
                else
                    gate = typeof(gate)((lanes(gate) .+ offset)...)
                end

                return Element{AbstractGate}(gate, priority)
            end
        end
    end |> Circuit
end