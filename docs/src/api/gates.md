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
Z() on 4
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

```@docs
Operator
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

### Rotation gates

```@docs
Rx
Ry
Rz
Rxx
Ryy
Rzz
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

## Special Unitary gate

!!! warn "Experimental interface"
    This interface is experimental and may change in the future.

```@docs
Quac.SU{N}
```
