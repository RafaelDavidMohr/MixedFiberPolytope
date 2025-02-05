# MixedFiberPolytope

A small julia package to compute mixed fiber polytopes as defined
[here](https://link.springer.com/article/10.1007/s11853-008-0015-2).

# Installation

The package may be installed by cloning this repository to `~/your_directory`
and then adding it to your julia environment:

```julia
julia> using Pkg
julia> Pkg.add("~/your_directory/MixedFiberPolytope")
julia> using MixedFiberPolytope
```

# Usage

To compute mixed fiber polytopes of Newton polyhedra of polynomials,
polynomials may be defined either using
[Oscar](https://www.oscar-system.org/) or
[DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl). The
computation is done by calling the function `mixed_fiber_polytope`. We
refer to its docstring in `src/MixedFiberPolytope.jl` for further
details.
