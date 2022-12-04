# Circuit

Circuits can be seen as DAGs (Direct Acyclic Graphs). In the case of quantum circuits, the width of the DAG is constant and equal to the number of qubits. Also the indgree and outdegree of quantum gates must be equal. Thus, using a graph for representing quantum circuits seems excesive because of its contraints.

Instead `Quac` uses **multi-priority queues** to store gates: there is a queue per qubit lane that stores the gates that act on it, and priorities are the order in which they are applied. If a gate acts on multiple qubits, it will contain a priority per qubit.
This data structure allows us to store gates in the most compact way while iterating on gates, create the reverse circuit, .... are still efficient. It seems to be the perfect fit for quantum circuits.

**Is this really necessary?** No, the bottleneck of quantum circuits is not on their representation but when reading the source code of other quantum circuit libraries, I wasn't convinced by their solutions: graphs, already laid out lists of lists, a serialized list of gates, ... It seemed like nobody could found the proper data structure for representing them. So I came up with multi-priority queues which seem like the perfect fit and as a consequence, the implementation is simple and efficient.

```@docs
Circuit
```