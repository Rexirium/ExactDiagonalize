# ExactDiagonalize.jl

A Julia package for exact diagonalization of quantum many-body systems. This package provides tools for constructing Hamiltonians, computing eigenspectra, and performing quantum dynamics simulations via time evolution.

## Overview

ExactDiagonalize enables computational studies of small to moderate-sized quantum systems through exact numerical methods. It supports:

- **Hamiltonian Construction**: Build quantum Hamiltonians from individual spin operators and composite operator sums
- **Spectrum Computation**: Diagonalize Hamiltonians to obtain eigenvalues and eigenvectors
- **Time Evolution**: Simulate quantum dynamics using the RK4 integration scheme
- **Observable Tracking**: Record time-dependent expectation values during evolution
- **Sparse Matrix Support**: Efficient sparse matrix representations of quantum operators

## Installation

Add ExactDiagonalize to your Julia project:

```julia
using Pkg
Pkg.add("ExactDiagonalize")
```

Or add directly from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/ExactDiagonalize.jl")
```

### Requirements

- Julia 1.12+
- Standard library: `LinearAlgebra`, `SparseArrays`
- Optional: `CairoMakie` for visualization

## Quick Start

### Basic XY Model Example

Simulate time evolution of an XY spin chain:

$$
    H = J \sum_{i = 1}^L \left(X_i X_{i+1} + Y_i Y_{i+1} + \Delta Z_i Z_{i+1} \right). 
$$

```julia
using ExactDiagonalize

# Set system type to spin
set_systype(:Spin)
L, N = 10, 1  # System size and particle number
Δ = 1.0       # Interaction strength

# Initial state: single excitation at site 1
init = NumState("1000000000")

# Build XY Hamiltonian
os = Tuple[]
for j in 1:L
    nj = mod1(j + 1, L)
    push!(os, (Δ, :Z, j, :Z, nj))        # ZZ interaction
    push!(os, (1.0, :X, j, :X, nj))      # XX coupling
    push!(os, (-1.0, :iY, j, :iY, nj))   # YY coupling
end
opsum = OpSum(os, Float64)

# Define observable and time points
obs = OperatorObserver((1.0, :Z, L), init.basis)
ts = 0.0:0.05:10.0

# Run time evolution
timeEvolve(opsum, init, ts, obs)
```

## Core API

### State Representation

- **`NumState`**: State represented by binary occupation string
- **`FullState`**: State defined on full dimension
- **`NumBasis`**: Basis generated from particle number occupation
- **`FullBasis`**: Basis defined on full dimension

### Operators

- **`SpinOp`**: Individual spin operators (`:X`, `:Y`, `:Z`, `:+`, `:-`, `:iY`)
- **`Operator`**: Multi-site operator products
- **`OpSum`**: Linear combinations of operators (Hamiltonian)

### Functions

- **`spectrum(opsum, basis)`**: Compute eigenvalues and eigenvectors
- **`exact(opsum, init_state, basis)`**: Exact diagonalization solution
- **`timeEvolve(opsum, init_state, times, observer)`**: RK4 time evolution
- **`record!(observer, state, time)`**: Record observable at given time
- **`makeHamiltonian(opsum, basis)`**: Convert OpSum to sparse matrix Hamiltonian

### System Configuration

- **`set_systype(type)`**: Set system type (`:Spin` or `:Fermion`) (`:Fermion` type is yet to be developed)
- **`get_systype()`**: Query current system type

## Examples

See the `examples/` directory for complete working examples:

- `xymodel.jl`: XY spin chain with time evolution and observable tracking

## Key Features

- **Efficient Sparse Representation**: Leverages sparse matrix formats for memory efficiency
- **Flexible Operator Syntax**: Intuitive specification of quantum Hamiltonians
- **Observable Recording**: Built-in framework for tracking time-dependent measurements
- **Hardware Acceleration**: Optional MKL support for accelerated linear algebra

## Project Structure

```
src/
├── ExactDiagonalize.jl   # Main module
├── exactdiag.jl          # Diagonalization core functions
├── operators.jl          # Operator and Hamiltonian construction
├── ode_solver.jl         # Time evolution (RK4)
├── sparsemat.jl          # Sparse matrix utilities
├── state_basis.jl        # State and basis definitions
└── utils.jl              # Helper utilities
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

See LICENSE file for details.

## Citation

If you use ExactDiagonalize in your research, please cite this repository.
