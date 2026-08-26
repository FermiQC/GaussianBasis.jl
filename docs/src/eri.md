```@meta
CurrentModule = GaussianBasis
```

# Two-Electron Integrals

## Four Centers

For the atomic orbitals ``\chi_\mu, \chi_\nu, \chi_\lambda, \chi_\sigma``, the two-electron four-center repulsion integral, in **chemist's notation**, is calculated as:

```math
(\mu\nu|\lambda\sigma) = \iint \chi_\mu(\mathbf{r}_1)\chi_\nu(\mathbf{r}_1)\,
\frac{1}{|\mathbf{r}_1-\mathbf{r}_2|}\, \chi_\lambda(\mathbf{r}_2)\chi_\sigma(\mathbf{r}_2)\, d\mathbf{r}_1 d\mathbf{r}_2
```

This integral obeys an 8-fold permutational symmetry,
``(\mu\nu|\lambda\sigma) = (\nu\mu|\lambda\sigma) = (\mu\nu|\sigma\lambda) =
(\lambda\sigma|\mu\nu) = \ldots``, and its full dense tensor scales as
``\mathcal{O}(n_\text{bas}^4)`` in memory -- prohibitive for anything but
small systems.

A full-dense tensor can be calculated directly from a basis set
```julia-repl
julia> bset = BasisSet("sto-3g", """
              O        0.000000000000     -0.143225816552      0.000000000000
              H        1.638036840407      1.136548822547     -0.000000000000
              H       -1.638036840407      1.136548822547     -0.000000000000
              """)
julia> I = ERI_2e4c(bset)

julia> size(I)
(7, 7, 7, 7)

julia> I[1, 1, 1, 1]   # (1s_O 1s_O | 1s_O 1s_O)
4.785065752279989

julia> a,b,c,d = 1, 3, 6, 4 # The array `I` has 8-fold symmetry
julia> I[a,b,c,d] == I[b,a,c,d] == I[a,b,d,c] == I[b,a,d,c] &&
       I[c,d,a,b] == I[d,c,a,b] == I[c,d,b,a] == I[d,c,b,a]
true
```
An in-place (mutating) version of the full-dense tensor can also be used:

```julia-repl
julia> out = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
julia> ERI_2e4c!(out, bset)
julia> out == I
true
```

For large basis sets, the full dense tensor quickly becomes intractable.
`sparseERI_2e4c` provides a middle ground: Cauchy-Schwarz shell-pair
screening discards quartets below a magnitude cutoff up front, and only the
permutationally-unique surviving elements are stored, rather than the full
dense tensor.

```julia-repl
julia> idx, vals = sparseERI_2e4c(bset)

julia> n = findfirst(==((1, 2, 2, 2)), idx) # For example reproducibility.
julia> vals[n]  # A non-zero integral value...
0.25663337137335623

julia> vals[n] ≈ I[idx[n]...]
true

julia> I[1,1,5,7] # Zero or near zero elements are not contained in `idx` and `vals`
0.0
julia> (1,1,5,7) in idx
false
julia> I[1,3,6,7]
-8.470329472543003e-22
julia> (1,3,6,7) in idx
false
julia> I[1,3,7,7] # All non zero elements (up to a threshold) are returned...
-0.0025845696440600086
julia> (1,3,7,7) in idx
true
julia> (3,1,7,7) in idx # But only one canonical permutation is stored.
false

julia> length(idx)   # number of unique, screened-surviving elements
228

julia> length(I) / length(idx) # Fraction of non-zero elements
10.530701754385966
```
The full ERI is thus ~10x larger than its sparse, screened counterpart.

> The ordering within `idx` and `vals` returned from `sparseERI_2e4c` is arbitrary.

In many cases, storing the full ERI array, even in sparse form, becomes impossible. In these cases, the per-shell-quartet option is also available:
```julia-repl
julia> I4443 = ERI_2e4c(bset, 4, 4, 4, 3) # 4 and 3 are shell indexes
1×1×1×3 Array{Float64, 4}:
[:, :, 1, 1] =
 0.026307122405698283

[:, :, 1, 2] =
 0.02055337661033345

[:, :, 1, 3] =
 0.0
```
For maximum efficiency, you may need to use the most primitive call, mutating and using a pre-allocated output
```julia-repl
julia> out = zeros(num_basis(bset[4]), num_basis(bset[4]), num_basis(bset[4]), num_basis(bset[3]))
julia> ERI_2e4c!(out, bset, 4, 4, 4, 3)
julia> out == I4443
true
```


```@docs
ERI_2e4c
ERI_2e4c!
sparseERI_2e4c
```


## Three Centers

For the atomic orbitals ``\chi_\mu, \chi_\nu`` of a "regular" orbital basis
and an auxiliary/fitting basis function ``\varphi_P``, the two-electron
three-center integral is calculated as:

```math
(\mu\nu|P) = \iint \chi_\mu(\mathbf{r}_1)\chi_\nu(\mathbf{r}_1) \frac{1}{|\mathbf{r}_1-\mathbf{r}_2|} \varphi_P(\mathbf{r}_2)\, d\mathbf{r}_1 d\mathbf{r}_2
```

This is one of the two building blocks for **density fitting** (also called
resolution-of-the-identity), which approximates the 4-center ERI by
expanding each electron's charge density in an auxiliary basis instead of
computing it exactly. It scales as ``\mathcal{O}(n_\text{bas}^2 n_\text{aux})``,
far more modestly than the ``\mathcal{O}(n_\text{bas}^4)`` full 4-center
tensor -- which is what makes density fitting an attractive approximation
for large systems.

A full-dense tensor can be calculated directly from an orbital basis set and
an auxiliary basis set
```julia-repl
julia> auxbset = BasisSet("cc-pvdz-rifit", """
              O        0.000000000000     -0.143225816552      0.000000000000
              H        1.638036840407      1.136548822547     -0.000000000000
              H       -1.638036840407      1.136548822547     -0.000000000000
              """)
julia> B = ERI_2e3c(bset, auxbset)

julia> size(B)
(7, 7, 84)

julia> B[1, 1, 1]
0.7085956990587061
```

Alternatively, an in-place (mutating) version can be used:
```julia-repl
julia> out = zeros(bset.nbas, bset.nbas, auxbset.nbas)
julia> ERI_2e3c!(out, bset, auxbset)
julia> out == B
true
```

A per-shell call is also available, but it requires a single basis set.
Hence, both your regular basis and your auxiliary basis must be merged, and
you must keep track of shell indexes carefully: shells `1:bset.nshells` of
the merged basis are the regular shells (unchanged), and shells beyond that
are the auxiliary ones, offset by `bset.nshells`.
```julia-repl
julia> merged = merge_basis(bset, auxbset)

julia> out = zeros(num_basis(bset[1]), num_basis(bset[1]), num_basis(auxbset[1]))
julia> ERI_2e3c!(out, merged, 1, 1, bset.nshells + 1) # shell 1 of auxbset is shell (bset.nshells + 1) of merged
julia> out[1, 1, 1] == B[1, 1, 1] # matches the full dense array from the example above
true
```


```@docs
ERI_2e3c
ERI_2e3c!
```


## Two Centers

For two auxiliary/fitting basis functions ``\varphi_P, \varphi_Q``, the
two-electron two-center integral is calculated as:

```math
(P|Q) = \iint \varphi_P(\mathbf{r}_1) \frac{1}{|\mathbf{r}_1-\mathbf{r}_2|} \varphi_Q(\mathbf{r}_2)\, d\mathbf{r}_1 d\mathbf{r}_2
```

This is the other **density fitting** building block: `(P|Q)` is the
Coulomb metric ``J_{PQ}`` used to fit the auxiliary-basis expansion
coefficients against the three-center integrals above. It only involves the
auxiliary basis, so it scales as ``\mathcal{O}(n_\text{aux}^2)``, far
cheaper than either the 3- or 4-center integral.

A full-dense matrix can be calculated directly from an auxiliary basis set
```julia-repl
julia> Jmetric = ERI_2e2c(auxbset)

julia> size(Jmetric)
(84, 84)

julia> Jmetric[1, 1]
0.1148022639511714
```

Alternatively, an in-place (mutating) version can be used:
```julia-repl
julia> out = zeros(auxbset.nbas, auxbset.nbas)
julia> ERI_2e2c!(out, auxbset)
julia> out == Jmetric
true
```
For maximum efficiency, you may need to use the most primitive call, mutating and using a pre-allocated output
```julia-repl
julia> out = zeros(num_basis(auxbset[3]), num_basis(auxbset[2]))
julia> ERI_2e2c!(out, auxbset, 3, 2) # 3 and 2 are shell indexes
1×1 Matrix{Float64}:
 0.7584111125059806
```


```@docs
ERI_2e2c
ERI_2e2c!
```
