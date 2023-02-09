# Quac

[![docs](https://img.shields.io/badge/docs-stable-blue)](https://bsc-quantic.github.io/Quac.jl)
[![CI](https://github.com/bsc-quantic/Quac.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/bsc-quantic/Quac.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/bsc-quantic/Quac.jl/branch/master/graph/badge.svg?token=D7KVE9SG9Z)](https://codecov.io/gh/bsc-quantic/Quac.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Registry](https://badgen.net/badge/registry/bsc-quantic/purple)](https://github.com/bsc-quantic/Registry)

`Quac` stands for _**Qua**ntum **c**ircuits_ and its a Julia library for quantum circuits with no assumptions about their use.

**_What does this means, you ask?_** Well, `Quac` is not a simulator, neither a controller of quantum computers. It just provides a `Circuit` data stracture, a set of gates and tools to manipulate them. Developers may use it as the core of their simulators or hardware controllers.

> ⚠️ Measurement gates are not currently supported as we are exploring how to fit non-unitary gates.

## Features

### Multiple representation of gates

Gates are symbolic in the sense that they do not store their representation. In `Quac` a gate just stores the lane in which it acts, and parameters if it's a parametric gate. Thanks to Julia's multiple-dispatch different representations can be queried lazily.

For example, this is a $Z$ that acts on qubit 4.

```julia
julia> using Quac
julia> gate = Z(4)
```

Any gate can be represented by a dense matrix.

```julia
julia> Matrix(gate)
2×2 Matrix{ComplexF32}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im
```

You can even specify the `eltype`!

```julia
julia> Matrix{Int}(gate)
2×2 Matrix{Int64}:
 1   0
 0  -1
```

Furthermore, the $Z$ gate allows a `Diagonal` representation!

```julia
julia> using LinearAlgebra
julia> Diagonal{Float32}(gate)
2×2 Diagonal{Float32, Vector{Float32}}:
 1.0    ⋅
  ⋅   -1.0
```

### Layout-agnostic `Circuit` representation

Quac uses multi-priority queues for representing `Circuit`s.

### SVG rendering of `Circuit`s

```julia
using Quac

circ = Quac.Algorithms.QFT(4)
draw(circ)
```

![Quantum Fourier Transform](docs/src/assets/qft.svg)

## Roadmap

- [ ] Gate decompositions
- [ ] ZX-calculus
- [ ] Spatial layouts
- [ ] Measurements
- [x] Visualization
- [ ] Support for _qudits_
