# Quac

`Quac` stands for _**Qua**ntum **c**ircuits_ and its a Julia library for quantum circuits with no assumptions about their use.

**_What does this means, you ask?_** Well, `Quac` is not a simulator, neither a controller of quantum computers. It just provides a `Circuit` data stracture, a set of gates and tools to manipulate them. Developers may use it as the core of their simulators or hardware controllers.

## Features
### Multiple representation of gates
TODO

## Example

### 4-qubit QFT
TODO

```julia
using Quac

circ = Circuit(4)
push!()
```

## Internals
### Gates
> ⚠️ Measurement gates are not currently supported as we are exploring how to fit non-unitary gates.


#### `AbstractGate` interface

If you want to create your own

### Circuit
Circuits can be seen as DAGs (Direct Acyclic Graphs).
In the case of quantum circuits, the width of the DAG is constant and equal to the number of qubits.

Instead `Quac` uses multi-priority queues to store gates.


## To do
- [ ] ZX-calculi
- [ ] Spatial layouts