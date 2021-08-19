module GaussianBasis

using Molecules
using Formatting
import Molecules: Atom, symbol

export BasisFunction, BasisSet

@doc raw"""
    GaussianBasis.BasisFunction

Object representing a shell of Gaussian basis functions composed of ``N`` primitives: 

```math
\chi_{l,m} = \sum_i^N C_i r^{l}e^{-\zeta_i r^2} Y_{l,m}(\theta,\phi)
```

# Fields

| Name    | Type | Description |
|:--------|:------|:-----------------------------------------------------------|
|`l`      |`Int32`            | Angular momentum number (e.g. 0, 1, 2 for S, P, D...)      |
|`coef`   |`Array{Float64,1}` | Array with coefficients (``C_i``) for each primitive       |
|`exp`    |`Array{Float64,1}` | Array with exponents (``\zeta_i``) for each primitive      |

# Examples

```julia
julia> bf = GaussianBasis.BasisFunction(1, [1/√2, 1/√2], [5.0, 1.2])
P shell with 3 basis built from 2 primitive gaussians

χ₁₋₁ =    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-1.2⋅r²)

χ₁₀  =    0.7071067812⋅Y₁₀⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₀⋅r¹⋅exp(-1.2⋅r²)

χ₁₁  =    0.7071067812⋅Y₁₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₁⋅r¹⋅exp(-1.2⋅r²)
```
"""
struct BasisFunction
    l::Cint
    coef::Array{Cdouble,1}
    exp::Array{Cdouble,1}
end

include("BasisParser.jl")
include("BasisSet.jl")
include("Libcint.jl")
include("Integrals.jl")
include("Gradients/FiniteDifferences.jl")
include("Gradients/OneElectronGrad.jl")

end # module
