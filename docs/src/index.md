# Quac.jl

`Quac` stands for _**Qua**_ntum _**c**_ircuits and it's a library for quantum circuit representation in Julia. It serves as the core library of multiple libraries related to quantum computing: simulators, hardware controllers, ...

**Quac** is not a simulator or a controller, but it provides the tools to build them.

## Features

- **Multiple representation of gates** Gates are symbolic in Quac, and thanks to Julia's multiple-dispatch, they posess multiple representations like dense arrays, diagonal matrices, ...
- **(New) compact and efficient description of circuits** Unlike other libraries, Quac uses **multi-priority queues** as the underlying data-structure.
- **SVG rendering of quantum circuits**

<p align="center"><img alt="4-qubit Quantum Fourier Transform" src="assets/qft.svg"></p>

- **Extensibility** Every part of the package is extensible with new types or functionality as you wish.

## Contents

```@contents
```
