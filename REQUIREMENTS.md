# ExactDiagonalize.jl - Requirements & Setup

## System Requirements

- **Julia**: 1.12 or later
- **Operating System**: Linux, macOS, or Windows
- **RAM**: Minimum 2GB (more recommended for larger systems)

## Package Dependencies

### Core Dependencies
These are essential for the package to function:

| Package | Version | Purpose |
|---------|---------|---------|
| `LinearAlgebra` | ≥ 1.12.0 | Linear algebra operations (eigenvalue decomposition, matrix multiplication) |
| `SparseArrays` | ≥ 1.12.0 | Sparse matrix storage and operations for efficient memory usage |

### Optional Dependencies

| Package | Version | Purpose | Use Case |
|---------|---------|---------|----------|
| `CairoMakie` | ≥ 0.11 | Interactive visualization and plotting | Visualizing results from examples |
| `HDF5` | ≥ 0.17 | File I/O in HDF5 format | Saving/loading simulation data |
| `LaTeXStrings` | ≥ 1.3 | LaTeX formatting in plots | Publication-quality figures |
| `MKL` | ≥ 0.9.1 | Intel Math Kernel Library | Hardware acceleration (improves performance) |
| `StaticArrays` | ≥ 1.9 | Stack-allocated arrays | Performance optimization |
| `Revise` | ≥ 3.5 | Interactive development | Development workflow (not needed for users) |

## Installation

### Basic Installation

```julia
using Pkg
Pkg.add("ExactDiagonalize")
```

### Development Installation

```julia
using Pkg
Pkg.develop("ExactDiagonalize")
```

### Installing Dependencies Manually

```julia
using Pkg
Pkg.instantiate()  # Installs all dependencies from Project.toml
```

## Performance Optimization

For significantly faster linear algebra operations, install Intel MKL:

```julia
using Pkg
Pkg.add("MKL")
```

The package automatically uses MKL when available. Performance improvements are especially noticeable for large Hamiltonians and time evolution simulations.

## Checking Your Setup

Verify that all dependencies are properly installed:

```julia
using ExactDiagonalize
using LinearAlgebra, SparseArrays
# If no errors appear, your setup is correct
```

## Troubleshooting

### Missing Dependency Error
If you see an error like `Package X not found`, install it with:
```julia
using Pkg
Pkg.add("X")
```

### Performance Issues
- Ensure `MKL` is installed for hardware acceleration
- Check Julia version: `julia --version` (should be 1.12+)
- Consider using `--threads=auto` flag when starting Julia

### Julia Installation

If Julia is not installed, download it from: https://julialang.org/downloads/

## Development Environment

For development, additionally install:

```julia
using Pkg
Pkg.add("Revise")       # Hot code reloading
Pkg.add("TestEnv")      # Testing utilities
```

Then use interactive development:
```julia
using Revise
include("src/ExactDiagonalize.jl")
```

Changes to source files will reload automatically.

## Dependency Graph

```
ExactDiagonalize
├── LinearAlgebra [required]
├── SparseArrays [required]
│
├── CairoMakie [optional]
├── HDF5 [optional]
├── LaTeXStrings [optional]
├── MKL [optional - recommended]
├── StaticArrays [optional]
└── Revise [optional - development]
```

## Version Compatibility

- **LinearAlgebra, SparseArrays**: Standard library, included with Julia
- **MKL**: Recommended version 0.9.1+
- **CairoMakie**: Version 0.11+ recommended for visualization
- **Other packages**: Use latest stable versions from Julia General registry
