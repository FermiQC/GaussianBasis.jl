import Base: getindex, show

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
struct BasisSet
    name::String
    atoms::Vector{Atom}
    basis::Vector{Vector{BasisFunction}}
    ind_offset::Vector{Int}
    natoms::Cint
    nbas::Cint
    nshells::Cint
    lc_atoms::Array{Cint,1}
    lc_bas::Array{Cint, 1}
    lc_env::Array{Cdouble,1}
end

function BasisSet(name::String, str_atoms::String)
    atoms = Molecules.parse_string(str_atoms)
    BasisSet(name, atoms)
end

function BasisSet(name::String, atoms::Vector{T}) where T <: Atom

    basis = Vector{BasisFunction}[]

    # Cache entries so repeated atoms don't need to be looked up all the time
    cache = Dict()
    for i in eachindex(atoms)
        S = symbol(atoms[i])
        if S in keys(cache)
            push!(basis, cache[S])
        else
            push!(basis, read_basisset(name, S))
            cache[S] = basis[end]
        end
    end
    BasisSet(name, atoms, basis)
end

function BasisSet(name::String, atoms::Vector{T}, basis::Vector{Vector{BasisFunction}}) where T <: Atom

    length(atoms) != length(basis) ? throw(ArgumentError("Error combining atoms and basis functions")) : nothing

    ATM_SLOTS = 6
    BAS_SLOTS = 8

    natm = length(atoms)
    ind_offset = zeros(Int, natm)
    nbas = 0
    nshells = 0
    nexps = 0
    nprims = 0

    for i in eachindex(atoms)
        ind_offset[i] = nshells
        for b in basis[i]
            nshells += 1
            nbas += 2*b.l + 1
            nexps += length(b.exp)
            nprims += length(b.coef)
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

    return BasisSet(name, atoms, basis, ind_offset, natm, nbas, nshells, lc_atm, lc_bas, env)
end

function BasisSet(BS1::BasisSet, BS2::BasisSet)
    name = BS1.name*"+"*BS2.name
    atoms = vcat(BS1.atoms, BS2.atoms)
    basis = vcat(BS1.basis, BS2.basis)

    return BasisSet(name, atoms, basis)
end

"""
    GaussianBasis.gto_norm(n::Signed, a::AbstractFloat)

Function that returns the normalization factor for a gaussian basis function rⁿ⋅exp(-a⋅r²)
"""
function gto_norm(n, a)
   # normalization factor of function rⁿ exp(-ar²)
    s = 2^(2n+3) * factorial(n+1) * (2a)^(n+1.5) / (factorial(2n+2) * √π)
    return √s
end

"""
    GaussianBasis.normalize_basisfunction!(B::BasisFunction)

Modify the BasisFunction object by normalizing each primitive.
"""
function normalize_basisfunction!(B::BasisFunction)
    for i = eachindex(B.coef)
        B.coef[i] *= gto_norm(B.l, B.exp[i])
    end
end

function getindex(B::BasisSet, N::Int)
    return B.basis[N]
end

function getindex(B::BasisSet, i::Int, j::Int)
    return B.basis[i][j]
end

function string_repr(B::BasisFunction)
    # Generate Unicode symbol for sub number
    l_sub = Char(0x2080 + B.l)

    # Unicode for superscript is a bit messier, so gotta use control flow
    l_sup = B.l == 1 ? Char(0x00B9) :
            B.l in [2,3] ? Char(0x00B1 + B.l) :
            Char(0x2070 + B.l)

    nbas = 2*B.l + 1
    mvals = collect(-B.l:B.l)
    nprim = length(B.exp)

    # Reverse Dict(symbol=>num) to get Symbols from B.l
    Lsymbol = Dict(value => key for (key, value) in AMDict)[B.l]
    out = "$(Lsymbol) shell with $nbas basis built from $nprim primitive gaussians\n\n"
    for m in mvals
        # Add sub minus sign (0x208B) if necessary
        m_sub = m < 0 ? Char(0x208B)*Char(0x2080 - m) : Char(0x2080 + m)
        out *= format("{:<4s} = ","χ$(l_sub)$m_sub")
        for i in eachindex(B.coef)

            if i > 1
                out *= B.coef[i] > 0 ? "\n     + " : "\n     - "
            end

            #out *= "$(abs(B.coef[i]))⋅Y$(l_sub)$m_sub"
            out *= format("{:>15.10f}⋅Y$(l_sub)$m_sub", abs(B.coef[i]))

            if B.l != 0 
                out *= "⋅r$l_sup"
            end

            out *= "⋅exp(-$(B.exp[i])⋅r²)"
        end
        out *="\n\n"
    end
    return strip(out)
end

function string_repr(B::BasisSet)
    out  =  "$(B.name) Basis Set\n"
    out *= "Number of shells: $(B.nshells)\n"
    out *= "Number of basis:  $(B.nbas)\n\n"

    l_to_symbol = Dict(
        0 => "s",
        1 => "p",
        2 => "d",
        3 => "f",
        4 => "g",
        5 => "h",
        6 => "i",
    )

    for i in eachindex(B.atoms)
        A = B.atoms[i]
        # Count how many times s,p,d appears for numbering
        count = zeros(Int16, 7)
        out *= "$(symbol(A)): "
        for bfs in B.basis[i]
            L = bfs.l
            count[L+1] += 1
            out *= "$(count[L+1])$(l_to_symbol[L]) "
        end
        out *="\n"
    end

    return strip(out)
end

# Pretty printing
function show(io::IO, ::MIME"text/plain", X::T) where T<:Union{BasisFunction, BasisSet}
    print(io, string_repr(X))
end
