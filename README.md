<p align="center">
  <img src="assets/gblogo.png" width="600" alt=""/>
</p>

<table align="center">
  <tr>
    <th>CI</th>
    <th>Coverage</th>
    <th>License</th>
  </tr>
  <tr>
    <td align="center">
      <a href=https://github.com/FermiQC/GaussianBasis.jl/actions/workflows/CI.yml>
      <img src=https://github.com/FermiQC/GaussianBasis.jl/actions/workflows/CI.yml/badge.svg>
      </a> 
    </td>
    <td align="center">
      <a href=https://codecov.io/gh/FermiQC/GaussianBasis.jl>
      <img src=https://codecov.io/gh/FermiQC/GaussianBasis.jl/branch/main/graph/badge.svg?token=JNouJPwoHm>
      </a> 
    </td>
    <td align="center">
      <a href=https://github.com/FermiQC/GaussianBasis.jl/blob/main/LICENSE>
      <img src=https://img.shields.io/badge/License-MIT-blue.svg>
      </a>
    </td>
  </tr>
</table>

GaussianBasis offers high-level utilities for molecular integral computations.

Current features include:

- Basis set parsing (`gbs` format)
- Standard basis set files from [BSE](https://www.basissetexchange.org/)
- One-electron integral (1e)
- Two-electron two-center integral (2e2c)
- Two-electrons three-center integral (2e3c)
- Two-electrons four-center integral (2e4c)
- Gradients (currrently under construction - *watch out!*)

Integral computations use by default the integral library [libcint](https://github.com/sunqm/libcint) *via* [libcint_jll.jl](https://github.com/JuliaBinaryWrappers/libcint_jll.jl). A simple Julia-written integral module `Acsint.jl` is also available, but it is significantly slower than the `libcint`.  

# Basic Usage

The simplest way to use the code is by first creating a `BasisSet` object. For example
```julia
julia> bset = BasisSet("sto-3g", """
              H        0.00      0.00     0.00                 
              H        0.76      0.00     0.00""")
sto-3g Basis Set
Type: Spherical   Backend: Libcint
Number of shells: 2
Number of basis:  2

H: 1s 
H: 1s
```
Next, call the desired integral function with the `BasisSet` object as the argument. Let's take the `overlap` function as an example:
```julia
julia> overlap(bset)
2×2 Matrix{Float64}:
 1.0       0.646804
 0.646804  1.0
```

| Function      | Description | Formula |
|---------------|-------------|:-------:|
| `overlap`       | Overlap between two basis functions | ![S](assets/ovlp.png)|
| `kinetic`       | Kinetic integral | ![T](assets/kin.png)|
| `nuclear`       | Nuclear attraction integral  | ![V](assets/nuc.png)|
| `ERI_2e4c`       | Electron repulsion integral - returns a full rank-4 tensor! | ![ERI](assets/4cERI.png)|
| `sparseERI_2e4c`       | Electron repulsion integral - returns non-zero elements along with a index tuple | ![sERI](assets/4cERI.png)|
| `ERI_2e3c`       | Electron repulsion integral over three centers. **Note:** this function requires another basis set as the second argument (that is the auxiliary basis set in [Density Fitting](http://vergil.chemistry.gatech.edu/notes/df.pdf)). It must be called as `ERI_2c3c(bset, aux)` | ![3cERI](assets/3cERI.png)|
| `ERI_2e2c`       | Electron repulsion integral over two centers  | ![2cERI](assets/2cERI.png)|

# Advanced Usage

## Basis Functions
`BasisFunction` object is the central data type within this package. Here, `BasisFunction` is an abstract type with two concrete structures: `SphericalShell` and `CartesianShell`. By default `SphericalShell` is created. In general a spherical basis function is

![BF](assets/bf.png)

where the sum goes over primitive functions. A `BasisFunction` object contains the data to reproduce the mathematical object, i.e. the angular momentum number (***l***), expansion coefficients (***c<sub>n</sub>***), and exponential factors (***&xi;<sub>n</sub>***). We can create a basis function by passing these arguments orderly:
```julia
julia> using StaticArrays
julia> atom = GaussianBasis.Atom(8, 16.0, [1.0, 0.0, 0.0])
julia> bf = BasisFunction(1, SVector(1/√2, 1/√2), SVector(5.0, 1.2), atom)
P shell with 3 basis built from 2 primitive gaussians

χ₁₋₁ =    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-1.2⋅r²)

χ₁₀  =    0.7071067812⋅Y₁₀⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₀⋅r¹⋅exp(-1.2⋅r²)

χ₁₁  =    0.7071067812⋅Y₁₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₁⋅r¹⋅exp(-1.2⋅r²)
```
We can now check the fields (attributes):
```julia
julia> bf.l
1

julia> bf.coef
2-element SVector{2, Float64} with indices SOneTo(2):
 0.7071067811865475
 0.7071067811865475

julia> bf.exp
2-element SVector{2, Float64} with indices SOneTo(2):
 5.0
 1.2
 ```
 Note that `exp` and `coef` are expected to be `SVector` from `StaticArrays`. 

 ## Basis Set

 The `BasisSet` object is the main ingredient for integrals. It can be created in a number of ways:

 - The highest level approach takes two strings as arguments, one for the basis set name and another for the XYZ file. See *Basic Usage*.

 - You can pass your vector of `Atom` structures instead of an XYZ string as the second argument. `GaussianBasis` uses the `Atom` structure from [Molecules.jl](https://github.com/FermiQC/GaussianBasis.jl).
  ```julia
atoms = GaussianBasis.parse_string("""
              H        0.00      0.00     0.00                 
              H        0.76      0.00     0.00""")
BasisSet("sto-3g", atoms)
```

 - Finally, instead of searching into `GaussianBasis/lib` for a basis set file matching the desired name, you can construct your own from scratch. We further discuss this approach below. 

 Basis sets are mainly composed of two arrays: a vector of atoms and a vector of basis functions objects. We can construct both manually for maximum flexibility: 
 ```julia
julia> h2 = GaussianBasis.parse_string(
   "H 0.0 0.0 0.0
    H 0.0 0.0 0.7"
)
2-element Vector{Atom{Int16, Float64}}:
 Atom{Int16, Float64}(1, 1.008, [0.0, 0.0, 0.0])
 Atom{Int16, Float64}(1, 1.008, [0.0, 0.0, 0.7])
 ```
Next, we create a vector of basis functions.
```julia
julia> shells = [BasisFunction(0, SVector(0.5215367271), SVector(0.122), h2[1]),
BasisFunction(0, SVector(0.5215367271), SVector(0.122), h2[2]),
BasisFunction(1, SVector(1.9584045349), SVector(0.727), h2[2])];
```
Finally, we create the basis set object. Note that, you got to make sure your procedure is consistent. The atoms used to construct the basis set object must be in the `atom` vector, otherwise unexpected results may arise. 
```julia
julia> bset = BasisSet("UnequalHydrogens", h2, shells)
UnequalHydrogens Basis Set
Type: Spherical{Molecules.Atom, 1, Float64}   Backend: Libcint
Number of shells: 3
Number of basis:  5

H: 1s 
H: 1s 1p
```
The most import fields here are:
```julia
julia> bset.name == "UnequalHydrogens"
true
julia> bset.basis == shells 
true
julia> bset.atoms == h2
true
```

### Integrals over different basis sets

Functions such as `ERI_2e3c` require two basis set as arguments. Looking at the corresponding equation
![3cERI](assets/3cERI.png) we see two basis set: ***&Chi;*** and ***P***. If your first basis set has 2 basis functions and the second has 4, your output array is a 2x2x4 tensor. For example
```julia
julia> b1 = BasisSet("sto-3g", """
              H        0.00      0.00     0.00                 
              H        0.76      0.00     0.00""")
julia> b2 = BasisSet("3-21g", """
              H        0.00      0.00     0.00                 
              H        0.76      0.00     0.00""")
julia> ERI_2e3c(b1,b2)
2×2×4 Array{Float64, 3}:
[:, :, 1] =
 3.26737  1.85666
 1.85666  2.44615

[:, :, 2] =
 6.18932  3.83049
 3.83049  5.60161

[:, :, 3] =
 2.44615  1.85666
 1.85666  3.26737

[:, :, 4] =
 5.60161  3.83049
 3.83049  6.18932
 ```
One electron integrals can also be employed with different basis set. 
```julia
julia> overlap(b1, b2)
2×4 Matrix{Float64}:
 0.914077  0.899458  0.473201  0.708339
 0.473201  0.708339  0.914077  0.899458

julia> kinetic(b1, b2)
2×4 Matrix{Float64}:
 1.03401  0.314867  0.20091  0.203163
 0.20091  0.203163  1.03401  0.314867
```
This can be useful when working with projections from one basis set onto another. 

### Computing integrals element-wise

For all integrals, you can get the full array by using the general syntax `integral(basisset)` (e.g. `overlap(bset)` or `ERI_2e4c(bset)`). Alternatively, you can specify a shell combination for which the integral must be computed
```julia
julia> ERI_2e4c(b1, 1,2,2,1)
1×1×1×1 Array{Float64, 4}:
[:, :, 1, 1] =
 0.2845189435761272

julia> kinetic(b1, 1,2)
1×1 Matrix{Float64}:
 0.2252049038643092
 ```
Mutating versions of the functions are also available 
```julia
julia> S = zeros(2,2);
julia> overlap!(S, b1)
julia> S
2×2 Matrix{Float64}:
 1.0       0.646804
 0.646804  1.0
 ```