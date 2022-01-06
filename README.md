# GaussianBasis.jl
[![CI](https://github.com/FermiQC/GaussianBasis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/FermiQC/GaussianBasis.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/FermiQC/GaussianBasis.jl/branch/main/graph/badge.svg?token=JNouJPwoHm)](https://codecov.io/gh/FermiQC/GaussianBasis.jl)

GaussianBasis offers a high-level wrap around the integral library [libcint](https://github.com/sunqm/libcint) *via* [libcint_jll.jl](https://github.com/FermiQC/Fermi.jl/issues/111). 

Install it using Julia's pkg manager:

```julia
pkg> add GaussianBasis
```

Current features include:

- Basis set parsing (`gbs` format)
- Standard basis set files from [BSE](https://www.basissetexchange.org/)
- One-electron integral (1e)
- Two-electron two-center integral (2e2c)
- Two-electrons three-center integral (2e3c)
- Two-electrons four-center integral (2e4c)
- Gradients (currrently under construction - *watch out!*)

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
| `overlap`       | Overlap between two basis functions | $\langle \chi_i \| \chi_j \rangle$|
| `kinetic`       | Kinetic integral | $\frac{1}{2}\langle \chi_i \| \vec{p} \cdot \vec{p}  \| \chi_j \rangle$|
| `nuclear`       | Nuclear attraction integral  | $\sum_A\langle \chi_i \|\frac{Z_A}{\|\vec{R}_A - \vec{r}\|} \|\chi_j \rangle$|
| `ERI_2e4c`       | Electron repulsion integral - returns a full rank-4 tensor! | $\left( \chi_i \chi_j \| \chi_k \chi_l\right)$|
| `sparseERI_2e4c`       | Electron repulsion integral - returns non-zero elements along with a index tuple | $\left( \chi_i \chi_j \| \chi_k \chi_l\right)$|
| `ERI_2e3c`       | Electron repulsion integral over three centers. **Note:** this function requires another basis set as the second argument (that is the auxiliary basis set in [Density Fitting](http://vergil.chemistry.gatech.edu/notes/df.pdf)). It must be called as `ERI_2c3c(bset, aux)` | $\left( \chi_i \chi_j \| P_k\right)$|
| `ERI_2e2c`       | Electron repulsion integral over two centers  | $\left( \chi_i \| \chi_j\right)$|

# Advanced Usage

## Basis Functions
`BasisFunction` object is the elementary data type within this package. In general, a basis function is $\chi_{m,l} = \sum_n c_n \cdot Y_{m,l} \cdot r^l \cdot e^{-\xi_n r^2}$
where the sum goes over primitive functions. A `BasisFunction` object contains the data to reproduce the mathematical object, i.e. angular momentum number ($l$), expansion coefficients ($c_n$), and exponential factors ($\xi_n$). We can create a basis function by passing these arguments orderly:
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

 `BasisSet` objects are the main ingredients for integrals. They can be created in a number of ways:

 - The "highest level" one: using two string arguments, one for the basis set name and another for the XYZ file. See *Basis Usage*.

 - Next, you can pass your vector of `Atom` structures instead of an XYZ string. `GaussianBasis` uses the `Atom` structure from [Molecules.jl](https://github.com/FermiQC/GaussianBasis.jl).
  ```julia
atoms = Molecules.parse_string("""
              H        0.00      0.00     0.00                 
              H        0.76      0.00     0.00""")
BasisSet("sto-3g", atoms)
```

 - Finally, instead of searching into `GaussianBasis/lib` for a basis set file matching the desired name, you can construct your own from scratch. We further discuss this approach below. 

 Basis sets are a map between atoms and their basis functions. Thus, the most important field here is a Vector (one dimensional array) of `BasisFunction` vectors (i.e. `Vector{Vector{BasisFunction}}`). Let us work out an example in the following lines.
 
 First we create a Vector of `Atom` objects

 ```julia
julia> h2 = Molecules.parse_string(
   "H 0.0 0.0 0.0
    H 0.0 0.0 0.7"
)
2-element Vector{Atom{Int16, Float64}}:
 Atom{Int16, Float64}(1, 1.008, [0.0, 0.0, 0.0])
 Atom{Int16, Float64}(1, 1.008, [0.0, 0.0, 0.7])
 ```
 This vector is ordered the same way as the input string. Next we create basis functions.
```julia
julia> s = BasisFunction(0, [0.5215367271], [0.122])
julia> p = BasisFunction(1, [1.9584045349], [0.727]);
```
Next, we create a map between atoms and basis function. In this case, for the sake of showing the flexibility here, let us do something unorthodox and attach one $s$ function to the first hydrogen and an $s$ and $p$ functions to the second:
```julia
julia> shells = [
    [s],  # A s function on the first hydrogen
    [s,p] # One s and one p function on the second hydrogen
];
```
Note that the "mapping" is simply done by the corresponding ordering, the $n$-th entry in `shells` is attributed to the $n$-th Atom in `h2`. Finally, we can create the basis set:
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
julia> bset.basis == shells
true
julia> bset.atoms == h2
true
```
Other fields (such as `bset.lc_env`) are mostly chewed up information for libcint. You can learn more about them [here](https://github.com/sunqm/libcint/blob/master/doc/program_ref.pdf).

