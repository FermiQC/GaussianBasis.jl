using GaussianBasis
using Combinatorics: doublefactorial
import Base: getindex, show

# struct with backend specific data for integral computation
abstract type IntLib end
# struct flag for computations using Asint.jl
struct ASint <: IntLib end
# struct flag for computations using Libcint.jl
struct LCint <: IntLib 
    atm::Vector{Cint}
    natm::Cint
    bas::Vector{Cint}
    nbas::Cint
    env::Vector{Cdouble}
end

function LCint(atoms::Vector{<:Atom}, basis::Vector{<:Vector{<:BasisFunction}})
    ATM_SLOTS = 6
    BAS_SLOTS = 8

    natm = length(atoms)

    nshells = 0
    nbas    = 0
    nexps   = 0
    nprims  = 0

    for i in eachindex(atoms)
        for b in basis[i]
            nshells += 1
            nbas    += num_basis(b)
            nexps   += length(b.exp)
            nprims  += length(b.coef)
        end
    end

    lc_atm = zeros(Cint, natm*ATM_SLOTS)
    lc_bas = zeros(Cint, nshells*BAS_SLOTS)
    env = zeros(Cdouble, 20+4*natm+nexps+nprims)

    # Prepare the lc_atom input 
    off = 20
    ib = 0 
    for i = eachindex(atoms)
        A = atoms[i]
        # lc_atom has ATM_SLOTS (6) "spaces" for each atom
        # The first one (Z_INDEX) is the atomic number
        lc_atm[1 + ATM_SLOTS*(i-1)] = A.Z
        # The second one is the env index address for xyz
        lc_atm[2 + ATM_SLOTS*(i-1)] = off
        env[off+1:off+3] .= A.xyz ./ Molecules.bohr_to_angstrom
        off += 4 # Skip an extra slot for the kappa (nuclear model parameter)
        # The remaining 4 slots are zero.

        for j = eachindex(basis[i])
            B = basis[i][j] 
            Ne = length(B.exp)
            Nc = length(B.coef)
            # lc_bas has BAS_SLOTS for each basis set
            # The first one is the index of the atom starting from 0
            lc_bas[1 + BAS_SLOTS*ib] = i-1
            # The second one is the angular momentum
            lc_bas[2 + BAS_SLOTS*ib] = B.l
            # The third is the number of primitive functions
            lc_bas[3 + BAS_SLOTS*ib] = Nc
            # The fourth is the number of contracted functions
            lc_bas[4 + BAS_SLOTS*ib] = 1
            # The fifth is a κ parameter
            lc_bas[5 + BAS_SLOTS*ib] = 0
            # Sixth is the env index address for exponents
            lc_bas[6 + BAS_SLOTS*ib] = off
            env[off+1:off+Ne] .= B.exp
            off += Ne
            # Seventh is the env index address for contraction coeff
            lc_bas[7 + BAS_SLOTS*ib] = off
            env[off+1:off+Nc] .= B.coef
            off += Nc
            # Eigth, nothing
            ib += 1
        end
    end
    return LCint(lc_atm, Cint(natm), lc_bas, Cint(nbas), env)
end

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
    basis::Vector{Vector{B}}
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

function BasisSet(name::String, atoms::Vector{Atom}; spherical=true, lib=:libcint)

    basis = spherical ? Vector{SphericalShell}[] : Vector{CartesianShell}[]

    # Cache entries so repeated atoms don't need to be looked up all the time
    cache = Dict()
    for i in eachindex(atoms)
        S = symbol(atoms[i])
        atom = atoms[i]
        if S in keys(cache)
            push!(basis, cache[S])
        else
            push!(basis, read_basisset(name, atom; spherical=spherical))
            cache[S] = basis[end]
        end
    end
    BasisSet(name, atoms, basis, lib=lib)
end

function BasisSet(name::String, atoms::Vector{Atom}, basis::Vector{<:Vector{<:BasisFunction}}; lib=:libcint)

    length(atoms) != length(basis) ? throw(ArgumentError("Error combining atoms and basis functions")) : nothing

    natm = length(atoms)

    # Number of shells per atom
    spa = length.(basis) 
    nshells = sum(spa)

    # Number of basis per atom
    bpa = zeros(Int, natm)
    for a in 1:natm
        bpa[a] = sum(num_basis(b) for b in basis[a])
    end
    nbas = sum(bpa)

    if lib == :asint

        return BasisSet(name, atoms, basis, bpa, spa, natm, nbas, nshells, ASint())

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