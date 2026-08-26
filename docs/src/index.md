```@meta
CurrentModule = GaussianBasis
```

# GaussianBasis.jl

[GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl) evaluates
molecular integrals (and their nuclear-coordinate derivatives) over
contracted Gaussian-type atomic orbitals, using [libcint](https://github.com/sunqm/libcint)
as its default backend. It is the integral engine behind
[Fermi.jl](https://github.com/FermiQC/Fermi.jl).


## Installation

```julia
using Pkg
Pkg.add("GaussianBasis")
```

## Quickstart

A [`BasisSet`](@ref) pairs a set of atoms with a set of contracted Gaussian
basis functions read from a standard basis set library (e.g. `sto-3g`,
`cc-pvdz`):

```julia-repl
julia> using GaussianBasis

julia> bset = BasisSet("sto-3g", """
       O        0.000000000000     -0.143225816552      0.000000000000
       H        1.638036840407      1.136548822547     -0.000000000000
       H       -1.638036840407      1.136548822547     -0.000000000000
       """);

sto-3g Basis Set
Type: Spherical   Backend: Libcint

Number of shells: 5
Number of basis:  7

O: 1s 2s 1p
H: 1s
H: 1s
```

Geometries are Cartesian coordinates in Angstrom. From a `BasisSet`, you can compute:

  - one-electron integrals (overlap, kinetic, nuclear attraction) -- see
    [One-Electron Integrals](@ref)
  - two-electron integrals (2-, 3-, and 4-center ERIs) -- see
    [Two-Electron Integrals](@ref)
  - multipole integrals (dipole, quadrupole, ...) -- see
    [Multipole Integrals](@ref)
  - nuclear-coordinate gradients and Hessians of all of the above -- see
    [Gradients](@ref) and [Hessians](@ref)

```julia-repl
julia> S = overlap(bset);
7×7 Matrix{Float64}:
 1.0         0.236704    0.0        0.0        0.0  0.00410862   0.00410862
 0.236704    1.0         0.0        0.0        0.0  0.0644883    0.0644883
 0.0         0.0         1.0        0.0        0.0  0.0572785   -0.0572785
 0.0         0.0         0.0        1.0        0.0  0.0447509    0.0447509
 0.0         0.0         0.0        0.0        1.0  0.0          0.0
 0.00410862  0.0644883   0.0572785  0.0447509  0.0  1.0          0.0100209
 0.00410862  0.0644883  -0.0572785  0.0447509  0.0  0.0100209    1.0
```

Electron Repulsion Integrals (ERI) are returned as full $N^4$ arrays.

```julia-repl
julia> ERI_2e4c(bset) |> size
(7, 7, 7, 7)
```

To accommodate memory sensitive operations, such as integral direct (on-the-fly) methods, integrals can be computed in shell batches of size $(2l+1)^4$, where $l$ is the angular momentum of the shell.

> The formula above is only valid for Spherical basis set. For more information check [`CartesianShell`](@ref).

```julia-repl
julia> ERI_2e4c(bset, 4,4,4,3) # Note that those are shell indexes, not basis.
1×1×1×3 Array{Float64, 4}:
[:, :, 1, 1] =
 0.026307122405698283

[:, :, 1, 2] =
 0.02055337661033345

[:, :, 1, 3] =
 0.0
```
For an even more allocation friendly option, check out [`ERI_2e4c!`](@ref)

## Basis Set

A basis set is a collection of basis functions that, in turn, represent one-electron wave functions for an electron bound to a nucleus. A [`BasisSet`](@ref) object can be conveniently created using a string of molecular coordinates and the name of a standard basis set. The available basis sets can be found in `lib/`.

```julia-repl
julia> BasisSet("cc-pvdz", """
              F        0.000000000000     -0.143225816552      0.000000000000
              H        1.638036840407      1.136548822547     -0.000000000000
""")
cc-pvdz Basis Set
Type: Spherical   Backend: Libcint

Number of shells: 9
Number of basis:  19

F: 1s 2s 3s 1p 2p 1d 
H: 1s 2s 1p
```

You may add new basis sets by simply creating the `.gbs` file inside `lib/`. These files can be found in the [Basis Set Exchange](https://www.basissetexchange.org/). GaussianBasis.jl uses the Psi4 format. 

In GaussianBasis.jl, a [`BasisSet`](@ref) object is the fundamental object necessary to compute integrals. In simple terms, a basis set is simply a container of Shells (see [`SphericalShell`](@ref) and [`CartesianShell`](@ref)), which represent a group of similar wave functions with the same $l$ value but different $m_l$ values.

```julia-repl
julia> for b in bset
           println(b)
       end
S shell on Fluorine (1 basis function, 9 primitives)
S shell on Fluorine (1 basis function, 9 primitives)
S shell on Fluorine (1 basis function, 1 primitive)
P shell on Fluorine (3 basis functions, 4 primitives)
P shell on Fluorine (3 basis functions, 1 primitive)
D shell on Fluorine (5 basis functions, 1 primitive)
S shell on Hydrogen (1 basis function, 4 primitives)
S shell on Hydrogen (1 basis function, 1 primitive)
P shell on Hydrogen (3 basis functions, 1 primitive)
```

>There is no individual basis function object. That is because integrals are calculated in shell batches. It is important to pay attention to the difference between number of basis and number of shells. Likewise, the index of a shell and the index of a basis are fundamentally different.

```@docs
BasisSet
SphericalShell
CartesianShell
GaussianBasis.atomic_orbital_amplitude
GaussianBasis.Libcint
GaussianBasis.on_atom_flags
merge_basis
```

## Contributing

Contributions in the form of bug reports, feature requests, and pull
requests are all welcome.

### Reporting Issues and Requesting Features

If you run into a bug, are looking for a feature, or simply need some help
with GaussianBasis.jl, feel free to open an
[issue on GitHub](https://github.com/FermiQC/GaussianBasis.jl/issues). You
can also reach out directly to [Gustavo Aroeira](https://github.com/gustavojra).

### Development Setup

Clone the repository and instantiate its environment:

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Running Tests

From the package root:

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

or, from a Julia REPL started with `--project=.`:

```julia-repl
julia> ]
pkg> test
```

### Building the Docs

The docs live under `docs/` and are built with
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). To build and
preview them locally:

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

This generates `docs/build/`; open `docs/build/index.html` in a browser to
preview your changes.

### Submitting a Pull Request

  - Branch off `main`.
  - Make sure `Pkg.test()` passes locally -- CI runs the same suite on
    Julia 1.9 and the latest 1.x release, on Linux and macOS.
  - Keep PRs minimal, do not bundle several unrelated changes into a single PR.
  - There's no enforced code formatter for this repo, so just match the style of the surrounding code.
  - Use of AI is allowed, but you must have reviewed the entire code and assured its quality.
  - Your new code must be tested. Make sure to wire your test suite into `test/runtests.jl`.
  - At the minimum, you should provide docstrings for your user-facing functions.