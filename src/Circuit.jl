import Base: push!, length

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
    new_priority = lanes(gate) .|> lane -> lane => circ.lanes[lane] + 1 |> collect
    el = Element(gate, new_priority)

    lanes(gate) .|> lane -> circ.lanes[lane] .|> queue -> push!(queue, el)
end
