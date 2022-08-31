# Quac

`Quac` stands for _**Qua**ntum **c**ircuits_ and its a Julia library for quantum circuits with no assumptions about their use.

**_What does this means, you ask?_** Well, `Quac` is not a simulator, neither a controller of quantum computers. It just provides a `Circuit` data stracture, a set of gates and tools to manipulate them. Developers may use it as the core of their simulators or hardware controllers.

> ⚠️ Measurement gates are not currently supported as we are exploring how to fit non-unitary gates.

## Features
### Multiple representation of gates

Gates are symbolic in the sense that they do not store their representation. In `Quac` a gate just stores the lane in which it acts, and parameters if it's a parametric gate. Thanks to Julia's multiple-dispatch different representations can be queried lazily.

For example, this is a $Z$ that acts on qubit 4.
```julia
> using Quac
> gate = Z(4)
```

Any gate can be represented by a dense matrix.
```julia
> Matrix(gate)
2×2 Matrix{ComplexF32}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im
```

You can even specify the `eltype`!
```julia
> Matrix{Int}(gate)
2×2 Matrix{Int64}:
 1   0
 0  -1
```

Furthermore, the $Z$ gate allows a `Diagonal` representation!
```julia
> using LinearAlgebra
> Diagonal{Float32}(gate)
2×2 Diagonal{Float32, Vector{Float32}}:
 1.0    ⋅
  ⋅   -1.0
```


## Example

### 4-qubit QFT
TODO

```julia
using Quac

circ = Circuit(4)
push!()
```

## Internals

### `AbstractGate` interface

Follow these instructions to implemente your own custom gate.

1. Set the parent abstract type to `AbstractGate`. Your struct should have a `lane` field of type `Int`.

```julia
struct CustomGate <: AbstractGate
	lane::Int
end
```
  - If your gate is a multi-qubit gate, then `lane` is of type `NTuple{N,Int}`.
  - If your gate is a parametric gate, then inherit from `AbstractParametricGate`.

2. Specify the type of the adjoint of your `CustomGate`. If your gate is hermitian, then it is itself.

```julia
Base.adjoint(::Type{CustomGate}) = CustomGate
```

3. Provide the representations of `CustomGate`. At least the `Matrix` representation should be provided, other possible representation is `Diagonal` if it allows it.

```julia
Matrix{T}(_::CustomGate) where {T} = Matrix{T}([...])
```

### Circuit
Circuits can be seen as DAGs (Direct Acyclic Graphs).
In the case of quantum circuits, the width of the DAG is constant and equal to the number of qubits.

Instead `Quac` uses multi-priority queues to store gates.


## To do
- [ ] ZX-calculi
- [ ] Spatial layouts