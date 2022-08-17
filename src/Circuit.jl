import Base: push!, length, iterate, IteratorSize

export Circuit
export lanes

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

isathead(e::Element, head::Vector{Int}) = all(p == head[lane] for (lane, p) in e.priority)

"""
A quantum circuit implementation using multi-priority queues.
- Queues are gate-buffers in qubit lanes.
- Multi-priority numbers can be retrieved procedurally from gate-lanes encoded inside the gates of the queues.
"""
struct Circuit
    lanes::Vector{Vector{Element{AbstractGate}}}

    Circuit(n::Int) = new(fill(Element[], n))
end

lanes(circ::Circuit) = length(circ.lanes)

Base.length(circ::Circuit) = sum(length(lane) for lane in circ.lanes)

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

"""
    Base.iterate(circ::Circuit[, state])

Retrieves next gate from `state` by travelling through a topologically sorted path.

# Arguments
- `circ::Circuit`
- `state` should be a `NTuple{N, Int}` where `N` is the number of lanes. Each element points to the head of the current cut.
"""
Base.iterate(circ::Circuit, state = fill(1, lanes(circ))) = begin
    # find winner
    candidates =
        enumerate(state) .|>
        (x -> begin
            lane, i = x
            get(circ.lanes[lane], i, nothing)
        end) |>
        (x -> filter(!isnothing, x))
    if isempty(candidates)
        return nothing
    end

    winner = filter(x -> isathead(x, state), candidates) |> first

    # update head
    for (lane, priority) in winner.priority
        state[lane] = priority + 1
    end

    (data(winner), state)
end

Base.IteratorSize(::Type{Circuit}) = Base.HasLength()