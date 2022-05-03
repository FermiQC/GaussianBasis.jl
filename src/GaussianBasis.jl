module GaussianBasis

using Molecules
using Formatting
import Molecules: Atom, symbol, parse_file, parse_string

export BasisFunction, BasisSet, SphericalShell, CartesianShell, get_shell

abstract type BasisFunction end

@doc raw"""
     BasisFunction

Basis object holding basis function information.

# Fields

| Name        | Type                               |   Description |
|:------------|:-----------------------------------|:-----------------------------------------------------------|
|`l`          | `Cint`                             | Angular momentum number |
|`coef`       | `Vector{Cdouble}`                  | Array holding expansion coefficients for the primitives |
|`exp`        | `Vector{Cdouble}`                  | Array holding exponents for primitives          |

# Examples

```julia
julia> bf = BasisFunction(1, [1/√2, 1/√2], [5.0, 1.2])
P shell with 3 basis built from 2 primitive gaussians

χ₁₋₁ =    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-1.2⋅r²)

χ₁₀  =    0.7071067812⋅Y₁₀⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₀⋅r¹⋅exp(-1.2⋅r²)

χ₁₁  =    0.7071067812⋅Y₁₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₁⋅r¹⋅exp(-1.2⋅r²)
```
"""
struct SphericalShell <: BasisFunction
    l::Int
    coef::Vector{Float64}
    exp::Vector{Float64}
end

struct CartesianShell <: BasisFunction
    l::Int
    coef::Vector{Float64}
    exp::Vector{Float64}
end

# Basis functions are created Spherical by default
BasisFunction(l::Int, coef::Vector{Float64}, exp::Vector{Float64}) = SphericalShell(l,coef,exp)

include("BasisParser.jl")
include("BasisSet.jl")
include("Libcint.jl")
include("Integrals.jl")
include("Gradients.jl")

end # module
