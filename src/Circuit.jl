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

Element(gate::Gate, priority) = Element{Gate}(gate, priority)

data(e::Element) = e.data

in(e::Element, head::Vector{Int}) = all(p == head[lane] for (lane, p) in e.priority)

"""
A quantum circuit implementation using multi-priority queues.

  - Queues are gate-buffers in qubit lanes.
  - Multi-priority numbers can be retrieved procedurally from gate-lanes encoded inside the gates of the queues.
"""
struct Circuit
    lanes::Vector{Vector{Element{Gate}}}

    Circuit(n::Int) = new(fill(Element[], n))
    Circuit(lanes::Vector{Vector{Element{Gate}}}) = new(lanes)
end

"""
    lanes(circuit)

Return the number of qubit lanes in a circuit.
"""
lanes(circuit::Circuit) = length(circuit.lanes)

"""
    length(circuit)

Return the number of gates in a circuit.
"""
Base.length(circuit::Circuit) = sum(length(lane) for lane in circuit.lanes; init = 0)

"""
    isempty(circuit)

Check whether the circuit contains any gate.
"""
Base.isempty(circuit::Circuit) = all(isempty, circuit.lanes)

"""
    push!(circuit, gate...)

Appends a gate to the circuit.
"""
Base.push!(circuit::Circuit, gate::Gate) = begin
    new_priority = [lane => length(circuit.lanes[lane]) + 1 for lane in lanes(gate)]
    el = Element{Gate}(gate, new_priority)

    for lane in lanes(gate)
        queue = circuit.lanes[lane]
        push!(queue, el)
    end
end

Base.push!(circuit::Circuit, gates::Gate...) = foreach(g -> push!(circuit, g), gates)

"""
    Base.iterate(circuit::Circuit[, state])

Retrieves next gate from `state` by travelling through a topologically sorted path.

# Arguments

  - `circuit::Circuit`
  - `state` (or head) should be a `NTuple{N, Int}` where `N` is the number of lanes. Each element is a pointer to the next gate on each lane.
"""
Base.iterate(circuit::Circuit, state = fill(1, lanes(circuit))) = begin
    # retrieve gates on the edge of the cut
    candidates =
        enumerate(state) |>
        (x -> filter(y -> ((lane, head) = y; head <= length(circuit.lanes[lane])), x)) .|>
        (x -> begin
            lane, i = x
            circuit.lanes[lane][i]
        end)
    if isempty(candidates)
        return nothing
    end

    # choose first valid gate
    winner = reduce(filter(x -> x ∈ state, candidates)) do a, b
        if mapreduce(x -> x[2], max, a.priority) <= mapreduce(x -> x[2], max, b.priority)
            return a
        else
            return b
        end
    end

    # update head by advancing cut on involved lanes
    for (lane, priority) in winner.priority
        state[lane] = priority + 1
    end

    (data(winner), state)
end

Base.IteratorSize(::Type{Circuit}) = Base.HasLength()

"""
    Base.adjoint(circuit)

Retrieve the adjoint circuit which fulfills the following equality.

```julia
circuit * circuit' == I(n)
```

# Notes

If all gates are hermitian, then the following equality also holds.

```julia
circuit * circuit' == circuit' * circuit == I(n)
```
"""
Base.adjoint(circuit::Circuit) = begin
    lanes = [
        reverse([
            Element{Gate}(data(el)', [laneid => length(circuit.lanes[laneid]) - p + 1 for (laneid, p) in el.priority]) for el in lane
        ]) for lane in circuit.lanes
    ]

    Circuit(lanes)
end

"""
    Base.getindex(circuit, lane, index)

Retrieve gate at lane `lane` and depth `index`.
"""
Base.getindex(circuit::Circuit, lane, index) = data(circuit.lanes[lane][index])

Base.show(io::IO, circuit::Circuit) = print(io, "Circuit(#lanes=$(lanes(circuit)), #gates=$(length(circuit)))")

"""
    connectivity([f,] circuit)

Generate connectivity graph between qubits.

# Arguments

  - f: Function to filter gates from circuit.
  - circuit: Circuit.
"""
function connectivity(f, circuit::Circuit)
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

connectivity(circuit::Circuit) = connectivity(gate -> length(lanes(gate)) >= 2, circuit)

"""
    moments(circuit)

Return moments (lists of gates that _can_ execute at the same time) of the circuit.
"""
function moments(circuit::Circuit)
    m = [Gate[]]

    for gate in circuit
        if isempty(last(m)) || isdisjoint(lanes(gate), ∪(lanes.(last(m))...))
            push!(last(m), gate)
        else
            push!(m, [gate])
        end
    end

    return m
end

"""
    hcat(circuits::Circuit...)

Join circuits in the temporal dimension.
"""
function Base.hcat(circuits::Circuit...)
    !allequal(lanes(circuit) for circuit in circuits) && throw(DimensionMismatch("circuits must have same lanes"))

    reduce(map(x -> x.lanes, circuits)) do acc, circuit
        map(zip(acc, circuit)) do (acclane, lane)
            offset = length(acclane)

            return [acclane..., map(lane) do el
                priority = [k => v + offset for (k, v) in el.priority]
                return Element{Gate}(data(el), priority)
            end...]
        end
    end |> Circuit
end

"""
    vcat(circuits::Circuit...)

Join circuits in the spatial dimension.
"""
function Base.vcat(circuits::Circuit...)
    offsets = [0, cumsum(lanes.(circuits))[1:end-1]...]

    mapreduce(vcat, zip(offsets, circuits)) do (offset, circuit)
        map(circuit.lanes) do lane
            map(enumerate(lane)) do (i, el)
                priority = [k + offset => v for (k, v) in el.priority]
                gate = data(el)

                if isparametric(gate)
                    gate = typeof(gate)((lanes(gate) .+ offset)..., parameters(gate))
                else
                    gate = typeof(gate)((lanes(gate) .+ offset)...)
                end

                return Element{Gate}(gate, priority)
            end
        end
    end |> Circuit
end