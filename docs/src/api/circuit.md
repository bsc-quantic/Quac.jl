# Circuit

Quantum circuits can be seen as DAGs (Direct Acyclic Graphs) in which the width of the DAG is constant and equal to the number of qubits. Also the indgree and outdegree of quantum gates must always be equal.
But many state-of-art quantum circuit libraries use either (a) lists of moments of gates or (b) graphs for representing them, which do not exploit the sparsity of the circuit or already make decisions about their layout (thus not having a layout-independent representation).

Instead `Quac` uses **multi-priority queues** to store gates: there is a queue per qubit lane that stores the gates that act on it, and priorities are the order in which they are applied. If a gate acts on multiple qubits, it will contain a priority per qubit.
This data structure allows us to store gates in the most compact way while iterating on gates, reversing the circuit, ... are still efficient. It appears to be the perfect data structure for quantum circuits.

!!! question Was this really necessary?

    No, the bottleneck of quantum circuits is not on their representation but when reading the source code of other quantum circuit libraries, I wasn't convinced by their solutions: graphs, already laid out lists of lists, a serialized list of gates, ... It seemed like nobody could found the proper data structure for representing them. So I came up with multi-priority queues which seem like the perfect fit and as a consequence, the implementation is simple yet efficient.

```@docs
Circuit
```

## Methods

```@index
Pages = ["circuit.md"]
Order = [:function]
```

```@docs
Base.adjoint(::Circuit)
Base.hcat(::Circuit...)
Base.isempty(::Circuit)
Base.iterate(::Circuit)
Base.length(::Circuit)
Base.push!(::Circuit, ::Gate)
Base.vcat(::Circuit...)
connectivity
lanes(::Circuit)
moments
```