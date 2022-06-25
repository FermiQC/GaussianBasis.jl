module GaussianBasis

using Molecules
using Formatting
using StaticArrays
import Molecules: Atom, symbol, parse_file, parse_string

export BasisFunction, BasisSet, SphericalShell, CartesianShell, get_shell

abstract type BasisFunction end

@doc raw"""
     SphericalShell

Object holding basis function information for an Spherical Shell.

# Fields

| Name        |   Description |
|:------------|:-----------------------------------------------------------|
|`l`          | Angular momentum number |
|`coef`       | Array holding expansion coefficients for the primitive functions |
|`exp`        | Array holding exponents for primitive functions          |

# Examples

```julia
julia> χ = SphericalShell(1, [1/√2, 1/√2], [5.0, 1.2])
P shell with 3 basis built from 2 primitive gaussians

χ₁₋₁ =    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-1.2⋅r²)

χ₁₀  =    0.7071067812⋅Y₁₀⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₀⋅r¹⋅exp(-1.2⋅r²)

χ₁₁  =    0.7071067812⋅Y₁₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₁⋅r¹⋅exp(-1.2⋅r²)
```
"""
struct SphericalShell{A<:Atom, N, R} <: BasisFunction
    l::Int
    coef::SVector{N, R}
    exp::SVector{N, R}
    atom::A
end

@doc raw"""
     CartesianShell

Object holding basis function information for an Cartesian Shell.

# Fields

| Name        |   Description |
|:------------|:-----------------------------------------------------------|
|`l`          | Pseudo-angular momentum number l = a+b+c for xᵃyᵇzᶜ |
|`coef`       | Array holding expansion coefficients for the primitive functions |
|`exp`        | Array holding exponents for primitive functions          |

# Examples

```julia
julia> χ = CartesianShell(2, [1/√2], [5.0])
D shell with 6 basis built from 1 primitive gaussians

χ(x²) =    0.7071067812⋅x²⋅exp(-5.0⋅r²)

χ(xy) =    0.7071067812⋅xy⋅exp(-5.0⋅r²)

χ(xz) =    0.7071067812⋅xz⋅exp(-5.0⋅r²)

χ(y²) =    0.7071067812⋅y²⋅exp(-5.0⋅r²)

χ(yz) =    0.7071067812⋅yz⋅exp(-5.0⋅r²)

χ(z²) =    0.7071067812⋅z²⋅exp(-5.0⋅r²)
```
"""
struct CartesianShell{A<:Atom, N, R} <: BasisFunction
    l::Int
    coef::SVector{N, R}
    exp::SVector{N, R}
    atom::A
end

# Basis functions are created Spherical by default
BasisFunction(l, coef, exp, atom) = SphericalShell(l,coef,exp, atom)
# Number of basis in a shell
num_basis(B::CartesianShell) = ((B.l + 1) * (B.l + 2)) ÷ 2
num_basis(B::SphericalShell) = 2*B.l + 1

function index2(i, j)
    if i < j
        return (j * (j + 1)) >> 1 + i
    else
        return (i * (i + 1)) >> 1 + j
    end
end

function index4(i,j,k,l)
    return index2(index2(i,j), index2(k,l))
end

include("BasisParser.jl")
include("BasisSet.jl")
include("Misc.jl")
include("Libcint.jl")
include("Acsint.jl")
include("Integrals.jl")
include("Gradients.jl")

end # module
