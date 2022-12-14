# Gates

```@meta
DocTestSetup = quote
    using Quac
    using LinearAlgebra
end
```

In `Quac`, gates are symbolic, i.e. they do not store their representation. A gate instance just stores the qubit lane in which it acts and its parameters if needed. Thanks to Julia's multiple-dispatch different representations can be queried lazily from type information.

For example, this is a $Z$ that acts on qubit 4.
```jldoctest z-gate
julia> gate = Z(4)
Z(4)
```

Any gate can be represented by a dense matrix.
```jldoctest z-gate
julia> Matrix(gate)
2×2 Matrix{ComplexF32}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im
```

You can even specify the `eltype`!
```jldoctest z-gate
julia> Matrix{Int}(gate)
2×2 Matrix{Int64}:
 1   0
 0  -1
```

Furthermore, the $Z$ gate allows a `Diagonal` representation!
```jldoctest z-gate
julia> Diagonal{Float32}(gate)
2×2 Diagonal{Float32, Vector{Float32}}:
 1.0    ⋅
  ⋅   -1.0
```

### `Gate` trait

All gates follow the `Gate` interface.

1. Set the parent abstract type to `Gate`. Your struct should have a `lane` field of type `Int`.

```julia
struct CustomGate <: Gate
	lane::Int
end
```
  - If your gate is a multi-qubit gate, then `lane` is of type `NTuple{N,Int}`.
  - If your gate is a parametric gate, then inherit from `ParametricGate`.

2. Specify the type of the adjoint of your `CustomGate`. If your gate is hermitian, then it is itself.

```julia
Base.adjoint(::Type{CustomGate}) = CustomGate
```

3. Provide the representations of `CustomGate`. At least the `Matrix` representation should be provided.

```julia
Matrix{T}(_::CustomGate) where {T} = Matrix{T}([...])
```

  - If the gate accepts other representations, you can implement them. For example, the $Z$ gate allows a `Diagonal`  representation.

```julia
Diagonal{T}(_::Z) where {T} = Diagonal{T}([1, -1])
```

```@docs
Gate
```

## Pauli gates
```@docs
Quac.I
X
Y
Z
```

## Hadamard gate
```@docs
H
```

## Phase gates
```@docs
S
Sd
T
Td
```

## Parametric gates
```@docs
ParametricGate
```

### Rotation gates
```@docs
Rx
Ry
Rz
```

### General U2, U3 gates
```@docs
U2
U3
```

## Controlled gates
```@docs
Control
```

## SWAP gate
```@docs
Swap
```