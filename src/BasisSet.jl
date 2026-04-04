using GaussianBasis
using Combinatorics: doublefactorial
import Base: getindex, show

include("Libs.jl")

@doc raw"""
    BasisSet

Object holding a set of BasisFunction objects associated with an array of atoms.

# Fields

| Name        | Type                               |   Description |
|:------------|:-----------------------------------|:-----------------------------------------------------------|
|`atoms`      | `Vector{Atom}`                     | An array of `Molecules.Atom` objects  |
|`name`       | `String`                           | String holding the basis set name  |
|`basis`      | `Vector{Vector{BasisFunction}}`    | An Array of arrays with BasisFunction          |
|`natoms`     | `Int32`                            | Number of atoms in the BasisSet |
|`nbas`       | `Int32`                            | Number of basis functions (Note, not equal the number of BasisFunction objects) |
|`nshells`    | `Int32`                            | Number of shells, i.e. BasisFunction objects |
|`lc_atoms`   | `Array{Int32,1}`                   | Integer array maping data to libcint |
|`lc_bas`     | `Array{Int32,1}`                   | Integer array maping data to libcint |
|`lc_env`     | `Array{Float64,1}`                 | Float64 array maping data to libcint |

# Example

Build a basis set from default options
```julia
julia> water = Molecules.parse_string(
"O        1.2091536548      1.7664118189     -0.0171613972
 H        2.1984800075      1.7977100627      0.0121161719
 H        0.9197881882      2.4580185570      0.6297938832"
)
julia> bset = BasisSet("sto-3g", water)
sto-3g Basis Set
Number of shells: 5
Number of basis:  7

O: 1s 2s 1p 
H: 1s 
H: 1s
```
The BasisSet object can be accessed as two-dimensional array.
```julia
julia> bset[1,2] # Show the second basis for the first atom (O 2s)
S shell with 1 basis built from 3 primitive gaussians

χ₀₀  =    0.8486970052⋅Y₀₀⋅exp(-5.033151319⋅r²)
     +    1.1352008076⋅Y₀₀⋅exp(-1.169596125⋅r²)
     +    0.8567529838⋅Y₀₀⋅exp(-0.38038896⋅r²)
```
You can also create your own crazy mix!
Lets create one S and one P basis functions for H
```julia
julia> h2 = Molecules.parse_string(
   "H 0.0 0.0 0.0
    H 0.0 0.0 0.7"
)
julia> s = BasisFunction(0, [0.5215367271], [0.122])
julia> p = BasisFunction(1, [1.9584045349], [0.727])
```
The basis set is constructed with an array of atoms (Vector{Atom}) and a corresponding array of Vector{BasisFunction}
holding all basis functions for that particular atom. In this example, we consider an unequal treatment for the
two atoms in the H₂ molecule.
```
julia> shells = [
    [s],  # A s function on the first hydrogen
    [s,p] # One s and one p function on the second hydrogen
]
julia> BasisSet("UnequalHydrogens", h2, shells)
UnequalHydrogens Basis Set
Number of shells: 3
Number of basis:  5

H: 1s 
H: 1s 1p
```
"""
struct BasisSet{L<:IntLib,A<:Atom,B<:BasisFunction}
    name::String
    atoms::Vector{A}
    basis::Vector{B}
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
    ) where {L<:IntLib,A<:Atom,B<:BasisFunction}
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

function construct_basis_from_library(::Type{B}, name::String, atoms::Vector{A}, ::Val{:libcint}) where {B<:BasisFunction,A<:Atom}
    basis = B[]
    for atom in atoms
        append!(basis, read_basisset(B, name, atom))
    end
    BasisSet(name, atoms, basis, LCint(atoms, basis))
end

function construct_basis_from_library(::Type{B}, name::String, atoms::Vector{A}, ::Val{:acsint}) where {B<:BasisFunction,A<:Atom}
    basis = B[]
    for atom in atoms
        append!(basis, read_basisset(B, name, atom))
    end
    BasisSet(name, atoms, basis, ACSint())
end

construct_basis_from_library(::Type{B}, name::String, atoms::Vector{A}, lib::Val) where {B<:BasisFunction,A<:Atom} =
    throw(ArgumentError("Unknown integral library: $lib"))


function BasisSet(name::String, atoms::Vector{A}, basis::Vector{B}, lib::L) where {A<:Atom,B<:BasisFunction,L<:IntLib}

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

    BasisSet(name, atoms, basis, bpa, spa, natm, nbas, nshells, lib)
end

BasisSet(name::String, atoms::Vector{A}, basis::Vector{B}) where {A<:Atom,B<:BasisFunction} =
    BasisSet(name, atoms, basis, LCint(atoms, basis))

BasisSet(name::String, atoms::Vector{A}; spherical::Bool=true, lib::Symbol=:libcint) where {A<:Atom} =
    build_basis_from_file(Val(spherical), name, atoms, lib)

function normalize_shell!(::Type{SphericalShell}, coef, exp, n)
    for i = eachindex(coef)
        a = exp[i]
        # normalization factor of function rⁿ exp(-ar²)
        s = 2^(2n+3) * factorial(n+1) * (2a)^(n+1.5) / (factorial(2n+2) * √π)
        coef[i] *= √s 
    end
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
    return B.basis[N]
end