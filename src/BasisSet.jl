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
struct BasisSet{L<:IntLib, A<:Atom, B<:BasisFunction} 
    name::String
    atoms::Vector{A}
    basis::Vector{B}
    basis_per_atom::Vector{Int}
    shells_per_atom::Vector{Int}
    natoms::Int
    nbas::Int
    nshells::Int
    lib::L
end

function BasisSet(name::String, str_atoms::String; spherical=true, lib=:libcint)
    atoms = Molecules.parse_string(str_atoms)
    BasisSet(name, atoms, spherical=spherical, lib=lib)
end

function BasisSet(name::String, atoms::Vector{<:Atom}; spherical=true, lib=:libcint)

    basis = spherical ? SphericalShell[] : CartesianShell[]

    for i in eachindex(atoms)
        bfs = read_basisset(name, atoms[i]; spherical=spherical)
        push!(basis, bfs...)
    end
    BasisSet(name, atoms, basis, lib=lib)
end

function BasisSet(name::String, atoms::Vector{<:Atom}, basis::Vector{<:BasisFunction}; lib=:libcint)

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

    if lib == :acsint

        return BasisSet(name, atoms, basis, bpa, spa, natm, nbas, nshells, ACSint())

    elseif lib == :libcint

        return BasisSet(name, atoms, basis, bpa, spa, natm, nbas, nshells, LCint(atoms, basis))
    else
        throw(ArgumentError("invalid integral backend option: $(lib)"))
    end
end

function BasisSet(BS1::BasisSet, BS2::BasisSet)
    name = BS1.name*"+"*BS2.name
    atoms = vcat(BS1.atoms, BS2.atoms)
    basis = vcat(BS1.basis, BS2.basis)

    return BasisSet(name, atoms, basis)
end

"""
    GaussianBasis.normalize_basisfunction!(B::SphericalShell)

Modify the SphericalShell object by normalizing basis function.
"""
function normalize_basisfunction!(B::SphericalShell)
    for i = eachindex(B.coef)
        n = B.l
        a = B.exp[i]
        # normalization factor of function rⁿ exp(-ar²)
        s = 2^(2n+3) * factorial(n+1) * (2a)^(n+1.5) / (factorial(2n+2) * √π)
        B.coef[i] *= √s 
    end
end

function normalize_spherical!(coef, exp, n)
    for i = eachindex(coef)
        a = exp[i]
        # normalization factor of function rⁿ exp(-ar²)
        s = 2^(2n+3) * factorial(n+1) * (2a)^(n+1.5) / (factorial(2n+2) * √π)
        coef[i] *= √s 
    end
end

function normalize_cartesian!(coef, exp, l)
    df = (π^1.5) * (l == 0 ? 1 : doublefactorial(2 * l - 1))

    # Normalize primitives
    coef .*= [ sqrt.((2 * ei) ^ (l + 1.5) / df) for ei in exp ]

    # Normalize contractions
    normsq = (df / 2^l) * sum([ ci * cj / (ei + ej) ^ (l + 1.5) 
                                           for (ci,ei) in zip(coef, exp)
                                           for (cj,ej) in zip(coef, exp)])
    coef .*= 1 / sqrt(normsq)
end

"""
    GaussianBasis.normalize_basisfunction!(B::CartesianShell)

Modify the CartesianShell object by normalizing basis function.
"""
function normalize_basisfunction!(B::CartesianShell)
    L = B.l
    df = (π^1.5) * (L == 0 ? 1 : doublefactorial(2 * L - 1))

    # Normalize primitives
    B.coef .*= [ sqrt.((2 * ei) ^ (L + 1.5) / df) for ei in B.exp ]

    # Normalize contractions
    normsq = (df / 2^L) * sum([ ci * cj / (ei + ej) ^ (L + 1.5) 
                                           for (ci,ei) in zip(B.coef, B.exp)
                                           for (cj,ej) in zip(B.coef, B.exp)])
    B.coef .*= 1 / sqrt(normsq)
end

function get_shell(BS::BasisSet, N::Int)
    idx = 1
    for A in 1:BS.natoms
        for b in BS.basis[A]
            if idx == N
                return BS.atoms[A], b
            end
            idx += 1
        end
    end
end

function getindex(B::BasisSet, N::Int)
    return B.basis[N]
end