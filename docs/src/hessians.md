```@meta
CurrentModule = GaussianBasis
```

# Hessians

Analytic nuclear Hessians of an integral are its second derivative with
respect to the Cartesian coordinates of two atoms ``A``, ``B`` (which may be
the same atom), e.g.

```math
\frac{\partial^2 S_{\mu\nu}}{\partial \mathbf{R}_A \partial \mathbf{R}_B}, \qquad
\frac{\partial^2 (\mu\nu|\lambda\sigma)}{\partial \mathbf{R}_A \partial \mathbf{R}_B}
```

returned as two extra trailing length-3 axes (``x,y,z`` for ``A``, then
``B``). These feed analytic vibrational frequency calculations and
second-order response (CPHF/CPKS) machinery.

!!! warning "Units: `R_A`/`R_B` are in bohr, not Angstrom"
    As with [Gradients](@ref), a [`BasisSet`](@ref) stores and accepts
    nuclear coordinates in Angstrom, but every Hessian here is with respect
    to bohr displacement of each nucleus. So `∇2overlap(bset, iA, iB)` gives
    `∂²S/∂R_iA∂R_iB` per bohr², *not* per Angstrom². To convert to a
    per-Angstrom² derivative, divide by `Molecules.bohr_to_angstrom^2`
    (≈0.529177 Å/bohr, squared).

Throughout this page the examples use

```julia-repl
julia> using GaussianBasis, Molecules

julia> water = Molecules.parse_string("""
       O        0.000000000000     -0.143225816552      0.000000000000
       H        1.638036840407      1.136548822547     -0.000000000000
       H       -1.638036840407      1.136548822547     -0.000000000000
       """);

julia> bset = BasisSet("sto-3g", water);   # 5 shells, 7 basis functions
```

## Overlap

```math
\frac{\partial^2 S_{\mu\nu}}{\partial \mathbf{R}_A \partial \mathbf{R}_B}
```

A full dense Hessian is computed directly from a basis set and two atom
indices:

```julia-repl
julia> hS = ∇2overlap(bset, 1, 2);   # ∂²S/∂R_O ∂R_H(first)

julia> size(hS)
(7, 7, 3, 3)
```

An in-place (mutating) version writes into a preallocated array:

```julia-repl
julia> out = zeros(bset.nbas, bset.nbas, 3, 3);

julia> ∇2overlap!(out, bset, 1, 2);

julia> out == hS
true
```

A per-shell-pair option is also available. This is usually what
integral-direct and CPHF code wants, since it never materializes the full
array:

```julia-repl
julia> b = ∇2overlap(bset, 1, 2, 4, 3)   # 4 and 3 are shell indexes
1×3×3×3 Array{Float64, 4}:
[:, :, 1, 1] =
 -0.00951483  -0.0358454  0.0
...
```

For maximum efficiency, use the mutating shell-pair form with a preallocated
output and a caller-owned `scratch` buffer, which makes the call
allocation-free:

```julia-repl
julia> Nmax = maximum(num_basis, bset.shells);

julia> scratch = Vector{Float64}(undef, 9*Nmax^2);

julia> out = zeros(num_basis(bset[4]), num_basis(bset[3]), 3, 3);

julia> ∇2overlap!(out, bset, 1, 2, 4, 3; scratch=scratch);

julia> out == b
true
```

```@docs
∇2overlap(::BasisSet, ::Any, ::Any)
∇2overlap!(::Any, ::BasisSet{LCint}, ::Any, ::Any, ::Int, ::Int)
```

## Kinetic Energy

```math
\frac{\partial^2 T_{\mu\nu}}{\partial \mathbf{R}_A \partial \mathbf{R}_B}
```

Identical structure to overlap -- `T` depends on the same two shell centers
and nothing else, so the same four levels are available:

```julia-repl
julia> hT = ∇2kinetic(bset, 1, 1);   # diagonal block, A = B = O

julia> hT[1, 1, :, :]
3×3 Matrix{Float64}:
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
```

That block is zero because shell 1 sits on atom 1, and translating a shell
together with its own atom leaves ``T`` unchanged.

```@docs
∇2kinetic(::BasisSet, ::Any, ::Any)
∇2kinetic!(::Any, ::BasisSet{LCint}, ::Any, ::Any, ::Int, ::Int)
```

## Nuclear-electron Attraction

```math
\frac{\partial^2 V_{\mu\nu}}{\partial \mathbf{R}_A \partial \mathbf{R}_B}
```

Nuclear attraction differs from overlap and kinetic in an important way: it
depends on **three** kinds of position -- the two shell centers *and* the
position of every nucleus supplying the potential. A derivative with respect
to atom ``A`` can therefore land on a shell center or on ``A``'s own nuclear
charge, so the second derivative splits into

1. both derivatives on shell centers, summed over all nuclei, and
2. at least one derivative on a moving nuclear charge, isolated one nucleus
   at a time via a repositionable point-charge operator.

`∇2nuclear` assembles both transparently:

```julia-repl
julia> hV = ∇2nuclear(bset, 1, 1);

julia> hV[1, 1, :, :]
3×3 Matrix{Float64}:
 -0.0284715   0.0         0.0
   0.0       -0.00452429  0.0
   0.0        0.0         0.0329958
```

Unlike overlap and kinetic, this block is **not** zero even though both
shells sit on atom 1 -- the potential's own center moves with the atom.

```@docs
∇2nuclear(::BasisSet, ::Any, ::Any)
∇2nuclear!
```

## Two-Electron Four Centers

```math
\frac{\partial^2 (\mu\nu|\lambda\sigma)}{\partial \mathbf{R}_A \partial \mathbf{R}_B}
```

There is deliberately **no dense form** here: the full tensor would be
``nbas^4 \times 3 \times 3``, which is not worth materializing for any
system where a Hessian is interesting. The shell-quartet form is the entry
point:

```julia-repl
julia> q = ∇2ERI_2e4c(bset, 1, 2, 3, 3, 1, 1);   # atoms 1,2; shells 3,3,1,1

julia> size(q)
(3, 3, 1, 1, 3, 3)
```

For a hot loop, use the mutating form with caller-owned scratch. The
`Xflag`/`Yflag` form skips the atom-membership lookup for callers that have
already screened:

```julia-repl
julia> using GaussianBasis: on_atom_flags

julia> X = on_atom_flags(bset, 1, 3, 3, 1, 1);   # which shells sit on atom 1

julia> Y = on_atom_flags(bset, 2, 3, 3, 1, 1);   # ... and on atom 2

julia> buf = Vector{Float64}(undef, 9*Nmax^4);

julia> out = zeros(size(q));

julia> ∇2ERI_2e4c!(out, bset, X, Y, 3, 3, 1, 1, buf);

julia> out == q
true
```

```@docs
∇2ERI_2e4c(::BasisSet, ::Int, ::Int, ::Int, ::Int, ::Int, ::Int)
∇2ERI_2e4c!(::Any, ::BasisSet, ::NTuple{4,Bool}, ::NTuple{4,Bool}, ::Int, ::Int, ::Int, ::Int, ::Vector{Float64})
```

## Two-Electron Three Centers

```math
\frac{\partial^2 (\mu\nu|P)}{\partial \mathbf{R}_A \partial \mathbf{R}_B}
```

Here ``\mu,\nu`` run over the regular basis and ``P`` over an
auxiliary/fitting basis (density fitting). A dense form exists but is
usually the wrong shape -- prefer the shell-triple form:

```julia-repl
julia> aux = BasisSet("cc-pvqz-jkfit", water);

julia> t = ∇2ERI_2e3c(bset, aux, 1, 2, 3, 3, 1);   # atoms 1,2; shells 3,3 and aux shell 1

julia> size(t)
(3, 3, 1, 3, 3)
```

In a loop, hoist both the scratch buffer and the merged basis (which depends
only on the two basis sets, never on the atoms or shells):

```julia-repl
julia> using GaussianBasis: merge_basis

julia> Bm = merge_basis(bset, aux);

julia> scratch = Vector{Float64}(undef, 9*Nmax^2*maximum(num_basis, aux.shells));

julia> out = zeros(size(t));

julia> ∇2ERI_2e3c!(out, bset, aux, 1, 2, 3, 3, 1; scratch=scratch, Bmerged=Bm);

julia> out == t
true
```

```@docs
∇2ERI_2e3c(::BasisSet, ::BasisSet, ::Any, ::Any)
∇2ERI_2e3c!(::Any, ::BasisSet, ::BasisSet, ::Any, ::Any, ::Int, ::Int, ::Int)
```

## Two-Electron Two Centers

```math
\frac{\partial^2 (P|Q)}{\partial \mathbf{R}_A \partial \mathbf{R}_B}
```

The density-fitting Coulomb metric. It depends on exactly two shell centers,
so it shares the overlap/kinetic structure -- only the basis set (auxiliary)
and the kernel family differ:

```julia-repl
julia> hJ = ∇2ERI_2e2c(aux, 1, 2);

julia> size(hJ)
(208, 208, 3, 3)
```

```@docs
∇2ERI_2e2c(::BasisSet, ::Any, ::Any)
∇2ERI_2e2c!(::Any, ::BasisSet{LCint}, ::Any, ::Any, ::Int, ::Int)
```

## Choosing a level

Every Hessian here is available at the same levels as the
[Gradients](@ref), from most convenient to fastest. All of them go through
the same libcint primitives, so they agree to the last bit -- they differ
only in how much bookkeeping is done for you.

| level | example | notes |
|:------|:--------|:------|
| dense, allocating | `∇2overlap(bs, A, B)` | fresh `nbas × nbas × 3 × 3` array; not offered for the 4-center ERI |
| dense, preallocated | `∇2overlap!(out, bs, A, B)` | reuses `out`; zeroes it for you |
| shell pair/triple/quartet | `∇2overlap!(out, bs, A, B, i, j)` | validates sizes, resolves shell membership, returns the free zero when the block vanishes |
| bare primitive | `∇2overlap_μμ!(out, bs, i, j)` | one libcint call, no checks at all |

Two things are worth knowing before dropping a level.

**A Hessian needs two primitives per integral, not one.** Where the
gradients express every case as the single `_μ!` primitive with permuted
arguments, a second derivative can put both derivatives on the same center
(`_μμ!`, libcint's `ipip` family) or one on each (`_μν!`, the `ipXip`
family). These are genuinely different kernels; neither is an argument
permutation of the other. The 3- and 4-center ERIs need more still, because
their centers fall into groups -- see the tables in their primitive
docstrings.

**Hoistable state should be hoisted.** Passing it in makes the call
allocation-free:

| routine | keyword | build it with |
|:--------|:--------|:--------------|
| `∇2overlap!`, `∇2kinetic!`, `∇2ERI_2e2c!` (shell pair) | `scratch` | `Vector{Float64}(undef, 9*Nmax^2)` |
| `∇2ERI_2e3c!` (shell triple) | `scratch`, `Bmerged` | [`merge_basis`](@ref) |
| `∇2ERI_2e4c!` (shell quartet) | `buf` (positional) | `Vector{Float64}(undef, 9*Nmax^4)` |

### Bare libcint primitives

The lowest level. Each does exactly one libcint call -- no bounds checking,
no zero-block skipping, no output-size validation. Unlike the gradient
primitives they apply **no sign flip**: libcint differentiates with respect
to the electron coordinate, each derivative contributes a factor of ``-1``,
and a second derivative applies two, which cancel.

!!! warning "Derivative-axis order"
    The cross primitives (`_μν!`, `_μλ!`, `_μP!`) return their two
    derivative-component axes in **(q,p)** order, reversed from the naive
    expectation. For same-center kernels this is invisible, since mixed
    partials of one point commute, but it is *not* invisible in general.
    Each primitive's docstring says which convention it uses.

```@docs
∇2overlap_μμ!
∇2overlap_μν!
∇2nuclear_rinv_μμ!
∇2ERI_2e2c_μμ!
∇2ERI_2e3c_μμ!
∇2ERI_2e4c_μμ!
```
