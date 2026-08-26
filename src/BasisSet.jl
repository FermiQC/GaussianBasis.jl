using GaussianBasis
using Combinatorics: doublefactorial
import Base: getindex, show, iterate, length, eltype

include("Libs.jl")

@doc raw"""
    BasisSet

Object holding a set of ShellFunction objects associated with an array of atoms.

# Fields

| Name        | Type                               |   Description |
|:------------|:-----------------------------------|:-----------------------------------------------------------|
|`atoms`      | `Vector{Atom}`                     | An array of `Molecules.Atom` objects  |
|`name`       | `String`                           | String holding the basis set name  |
|`shells`     | `Vector{Vector{ShellFunction}}`    | An Array of arrays with ShellFunction          |
|`natoms`     | `Int32`                            | Number of atoms in the BasisSet |
|`nbas`       | `Int32`                            | Number of basis functions.  |
|`nshells`    | `Int32`                            | Number of shells, i.e. ShellFunction objects |
|`lc_atoms`   | `Array{Int32,1}`                   | Integer array mapping data to libcint |
|`lc_bas`     | `Array{Int32,1}`                   | Integer array mapping data to libcint |
|`lc_env`     | `Array{Float64,1}`                 | Float64 array mapping data to libcint |

# Example

Build a basis set from default options
```julia-repl
julia> water = \""" 
 O        1.2091536548      1.7664118189     -0.0171613972
 H        2.1984800075      1.7977100627      0.0121161719
 H        0.9197881882      2.4580185570      0.6297938832\"""

julia> bset = BasisSet("sto-3g", water)
sto-3g Basis Set
Type: Spherical   Backend: Libcint

Number of shells: 5
Number of basis:  7

O: 1s 2s 1p 
H: 1s 
H: 1s
```
The BasisSet object can be accessed as vector, returning the n-th shell.
```julia-repl
julia> bset[1] # Show the first shell
S shell on Oxygen at position [1.2091536548, 1.7664118189, -0.0171613972] Å
Contains 1 basis function built from 3 primitive gaussians

χ₀₀  =    0.8486970052⋅Y₀₀⋅exp(-5.033151319⋅r²)
     +    1.1352008076⋅Y₀₀⋅exp(-1.169596125⋅r²)
     +    0.8567529838⋅Y₀₀⋅exp(-0.38038896⋅r²)
```
You can also create your own crazy mix!
Let us create one S and one P basis functions for H
```julia-repl
julia> H1 = GaussianBasis.Atom(:H, [0.0, 0.0, 0.0]);
julia> H2 = GaussianBasis.Atom(:H, [0.0, 0.0, 0.7]);
julia> s = ShellFunction(0, [0.5215367271], [0.122], H1)
julia> p = ShellFunction(1, [1.9584045349], [0.727], H2)
```
The basis set is constructed with an array of atoms (`Vector{Atom}`) and a corresponding array of `Vector{ShellFunction}`
holding all basis functions for that particular atom. In this example, we consider an unequal treatment for the
two atoms in the H₂ molecule.
```julia-repl
julia> shells = [s, p]  # One s function on the first hydrogen
                        # One p function on the second hydrogen
2-element Vector{SphericalShell{Molecules.Atom{Float64, Float64}}}:
 S shell on Hydrogen (1 basis function, 1 primitive)
 P shell on Hydrogen (3 basis functions, 1 primitive)
julia> BasisSet("UnequalHydrogens", shells)
UnequalHydrogens Basis Set
Type: Spherical   Backend: Libcint

Number of shells: 2
Number of basis:  4

H: 1s 
H: 1p
```

> Note that the number of shells and the number of basis functions are different. Integrals are computed over 
> shells, while the number of basis functions is the actual number of functions in the basis set and dictates
> the size of your integral arrays.
"""
struct BasisSet{L<:IntLib,A<:Atom,B<:ShellFunction}
    name::String
    atoms::Vector{A}
    shells::Vector{B}
    basis_per_atom::Vector{Int}
    shells_per_atom::Vector{Int}
    natoms::Int
    nbas::Int
    nshells::Int
    lib::L

    function BasisSet(
        name::String,
        atoms::Vector{A},
        basis::Vector{B},
        basis_per_atom::Vector{Int},
        shells_per_atom::Vector{Int},
        natoms::Int,
        nbas::Int,
        nshells::Int,
        lib::L
    ) where {L<:IntLib,A<:Atom,B<:ShellFunction}
        new{L,A,B}(name, atoms, basis, basis_per_atom, shells_per_atom, natoms, nbas, nshells, lib)
    end

end


function BasisSet(name::String, str_atoms::String; spherical=true::Bool, lib=:libcint::Symbol)
    atoms = Molecules.parse_string(str_atoms)
    return build_basis_from_file(Val(spherical), name, atoms, lib)
end

build_basis_from_file(::Val{true}, name::String, atoms::Vector{A}, lib::Symbol) where {A<:Atom} =
    construct_basis_from_library(SphericalShell, name, atoms, Val(lib))

build_basis_from_file(::Val{false}, name::String, atoms::Vector{A}, lib::Symbol) where {A<:Atom} =
    construct_basis_from_library(CartesianShell, name, atoms, Val(lib))

function construct_basis_from_library(::Type{B}, name::String, atoms::Vector{A}, ::Val{:libcint}) where {B<:ShellFunction,A<:Atom}
    basis = B[]
    for atom in atoms
        append!(basis, read_basisset(B, name, atom))
    end
    BasisSet(name, atoms, basis, LCint(atoms, basis))
end

function construct_basis_from_library(::Type{B}, name::String, atoms::Vector{A}, ::Val{:acsint}) where {B<:ShellFunction,A<:Atom}
    B <: SphericalShell && throw(ArgumentError(
        "the ACSint backend (lib=:acsint) only supports Cartesian shells -- pass spherical=false"
    ))
    basis = B[]
    for atom in atoms
        append!(basis, read_basisset(B, name, atom))
    end
    BasisSet(name, atoms, basis, ACSint())
end

construct_basis_from_library(::Type{B}, name::String, atoms::Vector{A}, lib::Val) where {B<:ShellFunction,A<:Atom} =
    throw(ArgumentError("Unknown integral library: $lib"))


function BasisSet(name::String, atoms::Vector{A}, basis::Vector{B}, lib::L) where {A<:Atom,B<:ShellFunction,L<:IntLib}

    natm = length(atoms)

    # Number of shells per atom
    nshells = length(basis)

    # Number of basis per atom
    bpa = zeros(Int, natm)
    spa = zeros(Int, natm)
    for a in 1:natm
        for b in basis
            if atoms[a] == b.atom
                spa[a] += 1
                bpa[a] += num_basis(b)
            end
        end
    end
    nbas = sum(bpa)

    # `atoms`/`basis` may have an abstract element type here (e.g.
    # Molecules.parse_string builds `Atom[]`, and `construct_basis_from_library`
    # builds `B[]` from a bare `::Type{B}`) even though every element pushed
    # into them is concrete. Narrowing via `identity.` re-infers the tightest
    # concrete common element type from the actual values, so `BasisSet`'s own
    # type parameters -- and every subsequent `.atoms[i]`/`.shells[i]` access,
    # which is on the hot path of every integral routine -- are concrete
    # instead of boxed.
    BasisSet(name, identity.(atoms), identity.(basis), bpa, spa, natm, nbas, nshells, lib)
end

BasisSet(name::String, atoms::Vector{A}, basis::Vector{B}) where {A<:Atom,B<:ShellFunction} =
    BasisSet(name, atoms, basis, LCint(atoms, basis))

# The atoms list is recovered from the shells themselves (deduplicated,
# first-occurrence order) since each shell already carries its own atom.
BasisSet(name::String, basis::Vector{B}) where {B<:ShellFunction} =
    BasisSet(name, unique(b.atom for b in basis), basis)

BasisSet(name::String, atoms::Vector{A}; spherical::Bool=true, lib::Symbol=:libcint) where {A<:Atom} =
    build_basis_from_file(Val(spherical), name, atoms, lib)

function normalize_shell!(::Type{SphericalShell}, coef, exp, n)

    # Self-overlap of the contraction. Basis-set files give coefficients
    # against *normalized* primitives, and the overlap of two of those is
    # (2√(aᵢaⱼ)/(aᵢ+aⱼ))^(n+3/2), so this double sum is ⟨φ|φ⟩. It must be
    # accumulated here, before the loop below folds the primitive
    # normalization into `coef` -- afterwards these are no longer the
    # coefficients of normalized primitives and the sum means nothing.
    #
    # Skipping this rescaling is only harmless when the file's contraction is
    # already normalized, which holds for many basis-set files but not all.
    # Basis Set Exchange's "Optimize General Contractions" output, for one,
    # drops a primitive from a contracted block without rescaling what
    # remains, leaving ⟨φ|φ⟩ far from 1. The Cartesian branch below has
    # always applied this step; the spherical branch had not.
    normsq = zero(eltype(coef))
    @inbounds for i = eachindex(coef), j = eachindex(coef)
        ei = exp[i]; ej = exp[j]
        normsq += coef[i] * coef[j] * (2*sqrt(ei*ej) / (ei + ej))^(n + 1.5)
    end

    # Normalize primitives
    for i = eachindex(coef)
        a = exp[i]
        # normalization factor of function rⁿ exp(-ar²)
        s = 2^(2n+3) * factorial(n+1) * (2a)^(n+1.5) / (factorial(2n+2) * √π)
        coef[i] *= √s
    end

    coef ./= sqrt(normsq)
end

function normalize_shell!(::Type{CartesianShell}, coef, exp, l)
    df = (π^1.5) * (l == 0 ? 1 : doublefactorial(2 * l - 1))

    # Normalize primitives
    coef .*= [ sqrt.((2 * ei) ^ (l + 1.5) / df) for ei in exp ]

    # Normalize contractions
    normsq = (df / 2^l) * sum([ ci * cj / (ei + ej) ^ (l + 1.5) 
                                           for (ci,ei) in zip(coef, exp)
                                           for (cj,ej) in zip(coef, exp)])
    coef .*= 1 / sqrt(normsq)
end

function getindex(B::BasisSet, N::Int)
    return B.shells[N]
end

# Iterating a BasisSet runs through its shells, the same unit every
# shell-indexed integral routine (e.g. ERI_2e4c(bset, i, j, k, l)) uses.
iterate(B::BasisSet, state::Int=1) = state > B.nshells ? nothing : (B[state], state + 1)
length(B::BasisSet) = B.nshells
eltype(::Type{BasisSet{L,A,Bf}}) where {L,A,Bf} = Bf