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

Integral computations use by default the integral library [libcint](https://github.com/sunqm/libcint) *via* [libcint_jll.jl](https://github.com/JuliaBinaryWrappers/libcint_jll.jl). A simple Julia-written integral module `Acsint.jl` is also available, but it is significant slower than the `libcint`.  

# Basic Usage

The simplest way to use the code is by first creating a `BasisSet` object. For example
```julia
julia> bset = BasisSet("sto-3g", """
              H        0.00      0.00     0.00                 
              H        0.76      0.00     0.00""")
sto-3g Basis Set
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
`BasisFunction` object is the elementary data type within this package. In general, a basis function is 

![BF](assets/bf.png)

where the sum goes over primitive functions. A `BasisFunction` object contains the data to reproduce the mathematical object, i.e. the angular momentum number (***l***), expansion coefficients (***c<sub>n</sub>***), and exponential factors (***&xi;<sub>n</sub>***). We can create a basis function by passing these arguments orderly:
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
We can now check the fields (attributes):
```julia
julia> bf.l
1

julia> bf.coef
2-element Vector{Float64}:
 0.7071067811865475
 0.7071067811865475

julia> bf.exp
2-element Vector{Float64}:
 5.0
 1.2
 ```
 Note that, the `BasisFunction` object has no information about the atom it is attached to! This information is introduced into the `BasisSet` struct.

 ## Basis Set

 The `BasisSet` object is the main ingredient for integrals. It can be created in a number of ways:

 - The highest level approach takes two strings as arguments, one for the basis set name and another for the XYZ file. See *Basis Usage*.

 - Youcan pass your vector of `Atom` structures instead of an XYZ string as the second argument. `GaussianBasis` uses the `Atom` structure from [Molecules.jl](https://github.com/FermiQC/GaussianBasis.jl).
  ```julia
atoms = Molecules.parse_string("""
              H        0.00      0.00     0.00                 
              H        0.76      0.00     0.00""")
BasisSet("sto-3g", atoms)
```

 - Finally, instead of searching into `GaussianBasis/lib` for a basis set file matching the desired name, you can construct your own from scratch. We further discuss this approach below. 

 Basis sets are a map between atoms and their basis functions. Thus, the most important field here is a Vector (one dimensional array) of `BasisFunction` vectors (i.e. `Vector{Vector{BasisFunction}}`). First we create a Vector of `Atom` objects
 ```julia
julia> h2 = Molecules.parse_string(
   "H 0.0 0.0 0.0
    H 0.0 0.0 0.7"
)
2-element Vector{Atom{Int16, Float64}}:
 Atom{Int16, Float64}(1, 1.008, [0.0, 0.0, 0.0])
 Atom{Int16, Float64}(1, 1.008, [0.0, 0.0, 0.7])
 ```
Next, we create basis functions.
```julia
julia> s = BasisFunction(0, [0.5215367271], [0.122])
julia> p = BasisFunction(1, [1.9584045349], [0.727]);
```
Next, we create a map between atoms and basis function. In this case, for the sake of showing the flexibility here, we will do something unorthodox and attach one $s$ function to the first hydrogen and an $s$ and $p$ functions to the second:
```julia
julia> bmap = [
    [s],  # One s function on the first hydrogen
    [s,p] # One s and one p function on the second hydrogen
];
```
Note that the "mapping" is simply achieved by the ordering, the $n$-th entry in `bmap` is attributed to the $n$-th Atom in `h2`. Finally, we can create the basis set:
```julia
julia> bset = BasisSet("UnequalHydrogens", h2, shells)
UnequalHydrogens Basis Set
Number of shells: 3
Number of basis:  5

H: 1s 
H: 1s 1p
```
The most import fields here are:
```julia
julia> bset.name == "UnequalHydrogens"
true
julia> bset.basis == bmap
true
julia> bset.atoms == h2
true
```
Other fields (such as `bset.lc_env`) are mostly chewed up information for `libcint`. You can learn more about them [here](https://github.com/sunqm/libcint/blob/master/doc/program_ref.pdf).

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

### Calling raw libcint functions

If one desires, it is possible to directly call the simplest/lowest-level wrap over libcint functions. This can be useful if instead of the full array you only want a few specific shell pairs. The list of function currently exposed can be found in the `GaussianBasis.Libcint` module and the corresponding `GaussianBasis/src/Libcint.jl` file.
```julia
julia> names(GaussianBasis.Libcint)
12-element Vector{Symbol}:
 :Libcint
 :cint1e_ipkin_sph!
 :cint1e_ipnuc_sph!
 :cint1e_ipovlp_sph!
 :cint1e_kin_sph!
 :cint1e_nuc_sph!
 :cint1e_ovlp_sph!
 :cint1e_r_sph!
 :cint2c2e_sph!
 :cint2e_ip1_sph!
 :cint2e_sph!
 :cint3c2e_sph!
 ```
If you need to expose a new function you should do in that file, don't forget to export the function afterwards!

 > A note on conventions: original `libcint` functions do not have **!** at the end of their names, we add this to comply with Julia practices. For example, `cint1e_ovlp_sph` is exposed as `cint1e_ovlp_sph!`. The **!** indicates that the function is a [mutating function](https://docs.julialang.org/en/v1/manual/style-guide/#Use-naming-conventions-consistent-with-Julia-base/).

 All functions are called with the same type signature. We shall use the overlap functions as the example here:

`cint1e_ovlp_sph!(buf, shls, atm, natm, bas, nbas, env)`

A more detailed description of each argument in found in the [`libcint` documentation](https://github.com/sunqm/libcint/blob/master/doc/program_ref.pdf). Here, we offer an overview

- `buf` type: `Array{Cdouble}`

Array where integral results are written into, it is understood as a linear array (though, it does not need to be passed as such); hence, you may need reshape to get what you want. 

- `shls` type `Array{Cint}`

Indicates the shell pair for which the integral must be computed. See example below.

- `atom` type: `Array{Cint}`

Holds information about atoms, mapping onto values in `env`.

- `natom` type `Cint`

Number of atoms

- `bas` type: `Array{Cint}`

Holds information about the basis set, mapping onto values in `env`.

- `nbas` type: `Cint`

Number of basis functions.

> Note that, **basis functions** here refers to the *actual* number of functions, e.g. for a system with *s* and *p* functions, there are 4 basis functions (because the *p* shell contains 3 basis functions). This can be confusing since the `BasisFunction` object in `GaussianBasis` actually holds information of the whole shell. Such is life.

- `env` type: `Array{Cdouble}`

Holds all the `Cdouble` (`Float64`) data. The `atm` and `bas` arrays tell you how to read chunks of this array.

#### Example

The `BasisSet` object contains the information formatted as required by `libcint`. Given a `BasisSet` object named `bset`, `atm`, `bas`, and `env` can be fetched through the fields `bset.lc_atoms`, `bset.lc_bas`, and `bset.lc_env`, respectively. Thus, it is recommended that you just use this object to call these functions. Let us work out an example for water using sto-3g.
```julia
bset = BasisSet("sto-3g", """
     O     1.2091536548    1.7664118189   -0.0171613972
     H     2.1984800075    1.7977100627    0.0121161719
     H     0.9197881882    2.4580185570    0.6297938832
""")
sto-3g Basis Set
Number of shells: 5
Number of basis:  7

O: 1s 2s 1p 
H: 1s 
H: 1s
```
Notice how the **shells** are ordered: O 1s, O 2s, O 1p, H 1s, H 1s

Suppose we want the overlap over O 1p and H 1s. Since p contains 3 functions, we need a 3x1 output array. The argument `shls` is used to indicate our choice of shells, in this case `shls = Cint.([2,3])`

> Counting starts from zero since this is a C call.

```julia
julia> buf = zeros(3,1)
julia> GaussianBasis.Libcint.cint1e_ovlp_sph!(buf, Cint.([2,3]), 
bset.lc_atoms, bset.natoms, bset.lc_bas, 
bset.nbas, bset.lc_env)
julia> buf
3×1 Matrix{Float64}:
 0.3813773519418131
 0.012065221257168891
 0.011286267412344286
```
Compare that with the **S** matrix:
```julia
julia> S = overlap(bset)
julia> S[3:5, 5]
julia> S[3:5, 6] ≈ buf
true
julia> S[6, 3:5] ≈ buf
true
```
