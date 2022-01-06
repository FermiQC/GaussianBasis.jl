# GaussianBasis.jl
[![CI](https://github.com/FermiQC/GaussianBasis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/FermiQC/GaussianBasis.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/FermiQC/GaussianBasis.jl/branch/master/graph/badge.svg?token=JNouJPwoHm)](https://codecov.io/gh/FermiQC/GaussianBasis.jl)

GaussianBasis offers a high-level wrap around the integral library [libcint](https://github.com/sunqm/libcint) *via* [libcint_jll.jl](https://github.com/FermiQC/Fermi.jl/issues/111). 

Install it using Julia's pkg manager:

```
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
```
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
```
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
