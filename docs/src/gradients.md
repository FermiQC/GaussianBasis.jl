```@meta
CurrentModule = GaussianBasis
```

# Gradients

## Overview

Analytic nuclear gradients of an integral are its first derivative with
respect to the Cartesian coordinates of a single atom `A`, e.g.
```math
\frac{\partial S_{\mu\nu}}{\partial \mathbf{R}_A}, \qquad
\frac{\partial (\mu\nu|\lambda\sigma)}{\partial \mathbf{R}_A}
```
returned as an extra length-3 trailing axis (``x,y,z``). For example:

```julia-repl
julia> bset = BasisSet("sto-3g", """
              H 0.0 0.0 0.0
              H 0.7 0.0 0.0   
              """)  # Note that the `x` is the bonding axis
julia> ∇overlap(bset, 1) # Derivative w.r.t the position of the first atom
2×2×3 Array{Float64, 3}:
[:, :, 1] =
 0.0       0.347146
 0.347146  0.0

[:, :, 2] =
 0.0  0.0
 0.0  0.0

[:, :, 3] =
 0.0  0.0
 0.0  0.0
```
The first sub-array (`[:, :, 1]`) is the derivative with respect to the `x` position
of the requested atom. In this simple example, displacements along `y` and `z` are rotations and do not affect 
overlap. Note also that the `x` sub-array is zero for shells on the same atom (such as diagonals).

!!! warning "Units: `R_A` is in bohr, not Angstrom"
    A [`BasisSet`](@ref) stores and accepts nuclear coordinates in Angstrom
    (see [Quickstart](@ref)), but every gradient/Hessian in this package is
    computed and returned in atomic units (bohr). To convert to a per-Angstrom
    derivative, divide by `Molecules.bohr_to_angstrom` (≈0.529177 Å/bohr),
    i.e. multiply by ≈1.8897 bohr/Å.

## Overlap

The derivative of the overlap integral is defined as follows:

```math
\frac{\partial S_{\mu\nu}}{\partial \mathbf{R}_A} = \frac{\partial\langle \chi_\mu | \chi_\nu \rangle}{\partial \mathbf{R}_A} 
```

If the shells are on the same atoms ($\mu, \nu \in A$) or if neither shell is on the atom ($\mu, \nu \notin A$), the derivative is exacly zero. 

Functions follow the same pattern found in [One-Electron Integrals](@ref), except you must include a mandatory argument `A` which indicates the atom for each derivatives are evaluated.
```julia-repl
julia> dr = ∇overlap(bset, 2, 1, 2) # Second atom, first shell, second shell
julia> dr[1,1,:] # Shows derivatives along x,y, and z
3-element Vector{Float64}:
 -0.3471463589514361
  0.0
  0.0
```

```@docs
∇overlap
∇overlap!
```

## Kinetic

The derivative of the kinetic energy integral is defined as follows:

```math
\frac{\partial T_{\mu\nu}}{\partial \mathbf{R}_A} = \langle \frac{\partial \chi_\mu}{\partial \mathbf{R}_A} | -\frac{1}{2}\nabla^2 | \chi_\nu \rangle + \langle \chi_\mu | -\frac{1}{2}\nabla^2 | \frac{\partial \chi_\nu}{\partial \mathbf{R}_A} \rangle
```

As with overlap, if the shells are on the same atom ($\mu, \nu \in A$) or if neither shell is on the atom ($\mu, \nu \notin A$), the derivative is exactly zero.

Functions follow the same pattern found in [One-Electron Integrals](@ref), except you must include a mandatory argument `A` which indicates the atom for which derivatives are evaluated.
```julia-repl
julia> dr = ∇kinetic(bset, 2, 1, 2) # Second atom, first shell, second shell
julia> dropdims(dr, dims=(1,2)) # Shows derivatives along x, y, and z
3-element Vector{Float64}:
 -0.33446292736778443
 -0.0
 -0.0
```

```@docs
∇kinetic
∇kinetic!
```

## Nuclear

The nuclear attraction integral sums the potential of every nucleus $C$ in
the molecule:

```math
V_{\mu\nu} = -\sum_C Z_C \langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \chi_\nu \rangle
```

Its derivative w.r.t. atom $C$ has **three** contributing terms instead of
two: the usual pair of shell-center ("bra"/"ket") derivatives from moving
$\chi_\mu$/$\chi_nu$ summed over every nucleus `C` plus an extra term from the operator itself, since the `C = A` piece of the
potential, $Z_A/|r-R_A|$, depends on $R_A$ too:

```math
\frac{\partial V_{\mu\nu}}{\partial \mathbf{R}_A} =
-\sum_C Z_C \left[ \langle \frac{\partial \chi_\mu}{\partial \mathbf{R}_A} | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \chi_\nu \rangle
+ \langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \frac{\partial \chi_\nu}{\partial \mathbf{R}_A} \rangle \right]
- Z_A \langle \chi_\mu | \frac{\partial}{\partial \mathbf{R}_A}\frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle
```

#### ${\mu,\nu} \in A$
For the case in which both shells are on the atom, we have
```math
\frac{\partial V_{\mu\nu}}{\partial \mathbf{R}_A} =
-\sum_{C\notin A} Z_C \left[ \langle \frac{\partial \chi_\mu}{\partial \mathbf{R}_A} | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \chi_\nu \rangle
+ \langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \frac{\partial \chi_\nu}{\partial \mathbf{R}_A} \rangle \right]
- Z_A \frac{\partial}{\partial \mathbf{R}_A}\langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle
```
but translation invariance requires
```math
\frac{\partial}{\partial \mathbf{R}_A}\langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle = 0
```
Hence,
```math
\frac{\partial V_{\mu\nu}}{\partial \mathbf{R}_A} =
-\sum_{C\notin A} Z_C \left[ \langle \frac{\partial \chi_\mu}{\partial \mathbf{R}_A} | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \chi_\nu \rangle
+ \langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \frac{\partial \chi_\nu}{\partial \mathbf{R}_A} \rangle \right]
```
#### ${\mu,\nu} \notin A$

When no shell is at the atom, the formula simplifies to
```math
\frac{\partial V_{\mu\nu}}{\partial \mathbf{R}_A} =
- Z_A \langle \chi_\mu | \frac{\partial}{\partial \mathbf{R}_A}\frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle
```
Translational invariance still applies, but with a caveat:
```math
\frac{\partial}{\partial \mathbf{R}_A}\langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle = 0 \\
 \langle \frac{\partial \chi'_\mu}{\partial \mathbf{R}_A} | \frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle
+ \langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \frac{\partial \chi'_\nu}{\partial \mathbf{R}_A} \rangle
+ \langle \chi_\mu | \frac{\partial}{\partial \mathbf{R}_A}\frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle = 0
```
The trick here is that $\chi_\mu'$ and $\chi_\nu'$ are moving with the atom $A$. This is artificial, but computable. Hence:
```math
\frac{\partial V_{\mu\nu}}{\partial \mathbf{R}_A} =
Z_A \left[  \langle \frac{\partial \chi'_\mu}{\partial \mathbf{R}_A} | \frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle
+ \langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \frac{\partial \chi'_\nu}{\partial \mathbf{R}_A} \rangle\right]
```

#### $\mu \in A$ and $\nu \; \notin A$

```math
\frac{\partial V_{\mu\nu}}{\partial \mathbf{R}_A} =
-\sum_C Z_C \langle \frac{\partial \chi_\mu}{\partial \mathbf{R}_A} | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \chi_\nu \rangle
- Z_A \langle \chi_\mu | \frac{\partial}{\partial \mathbf{R}_A}\frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle
```


Unlike overlap and kinetic, this is **never** exactly zero. The first two
(shell-derivative) terms still vanish unless `μ`/`ν` sits on atom `A` --
same translational-invariance argument as before -- but the third
(operator-derivative) term survives regardless of which atoms `μ`/`ν` sit
on, since nucleus `A`'s potential reaches every AO pair. For example, with
shells 1 and 1 both sitting on atom 1 (not atom 2):

```julia-repl
julia> dropdims(∇overlap(bset, 2, 1, 1), dims=(1,2))  # shell-derivative terms only -- exactly zero
3-element Vector{Float64}:
 0.0
 0.0
 0.0

julia> dropdims(∇nuclear(bset, 2, 1, 1), dims=(1,2))  # operator-derivative term survives
3-element Vector{Float64}:
 0.36265343354910295
 0.0
 0.0
```

###


### Reducing it to two libcint calls

The operator-derivative term above isn't something libcint can compute
directly -- `cint1e_ipnuc_sph!` only differentiates a *shell*, never the
potential's own center. So instead the code rewrites the whole
three-term expression as a sum of two shell-derivative pieces, split by
whether each nucleus's charge is `A`'s or not.

**Split the sum at `C = A`:**

```math
V_{\mu\nu} = V^A_{\mu\nu} + V^{\neg A}_{\mu\nu}, \qquad
V^A_{\mu\nu} \equiv -Z_A \langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_A|} | \chi_\nu \rangle, \qquad
V^{\neg A}_{\mu\nu} \equiv -\sum_{C \notin A} Z_C \langle \chi_\mu | \frac{1}{|\mathbf{r}-\mathbf{R}_C|} | \chi_\nu \rangle
```

**Define the two libcint-computable pieces.** For a shell pair `(p,q)` and a
charge set (either "just `A`" or "every nucleus except `A`"), one
`cint1e_ipnuc_sph!` call (using a charge array with every *other* nucleus's
charge zeroed out) gives the bra-shell electron-coordinate derivative of the
*very same* quantity `cint1e_nuc_sph!` itself would return for that charge
set -- `V^A_{pq}` or `V^{¬A}_{pq}` from above, nothing separately signed:

```math
K_A(p,q) \equiv \frac{\partial V^A_{pq}}{\partial \mathbf{r}_p}, \qquad
K_{\neg A}(p,q) \equiv \frac{\partial V^{\neg A}_{pq}}{\partial \mathbf{r}_p}
```

This is the one rule every `ip`-prefixed libcint kernel follows, uniformly
across `ovlp`, `kin`, and `nuc`: `cint1e_ipXXX_sph!` is the plain
electron-coordinate derivative of whatever `cint1e_XXX_sph!` itself
returns, no hidden extra sign. The translation identity
`∂χ/∂R_p = -∂χ/∂r` then converts this to the nuclear-coordinate derivative,
exactly the same single flip used for overlap and kinetic:

```math
\frac{\partial V^A_{pq}}{\partial \mathbf{R}_p} = -K_A(p,q), \qquad
\frac{\partial V^{\neg A}_{pq}}{\partial \mathbf{R}_p} = -K_{\neg A}(p,q)
```

`∇overlap_μ!`/`∇kinetic_μ!` apply that same flip by hand (`.*= -1.0`)
because `S`/`T` carry no leading sign of their own to fold it into;
`∇nuclear_μ!` doesn't need to, only because `V^A`/`V^{¬A}`
*already* carry a leading minus (from `Z` being stored positive but the
potential being attractive) for the flip to land inside of -- not because
`cint1e_ipnuc_sph!` is doing anything libcint's other kernels don't.

**Piece 1 -- `V^{¬A}`.** `R_A` never appears in this operator, only possibly
in the shell centers, so plain chain rule plus the relation above gives

```math
\frac{\partial V^{\neg A}_{\mu\nu}}{\partial \mathbf{R}_A} =
\frac{\partial V^{\neg A}_{\mu\nu}}{\partial \mathbf{R}_\mu}\Big|_{\mu \in A}
+ \frac{\partial V^{\neg A}_{\mu\nu}}{\partial \mathbf{R}_\nu}\Big|_{\nu \in A}
= -K_{\neg A}(\mu,\nu)\Big|_{\mu \in A} \;-\; K_{\neg A}(\nu,\mu)^\top\Big|_{\nu \in A}
```
(each term present only if that shell sits on `A`.)

**Piece 2 -- `V^A`.** Here `R_A` can appear in up to *three* places: the
operator itself, and (if the shells sit there) `R_μ`, `R_ν`. Treat
`V^A(R_μ,R_ν,R_A) = -Z_A⟨χ_μ|1/|r-R_A||χ_ν⟩` as a function of three
independent positions -- shifting all three together leaves the integral
unchanged (same argument as overlap's two-position version, one position
more), so

```math
\frac{\partial V^A_{\mu\nu}}{\partial \mathbf{R}_\mu} + \frac{\partial V^A_{\mu\nu}}{\partial \mathbf{R}_\nu} + \frac{\partial V^A_{\mu\nu}}{\partial \mathbf{R}_A} = 0
```

The full derivative with respect to atom `A`'s position -- dragging along
whichever shells are attached to it, plus the operator's own dependence,
substituted from the invariance above:

```math
\frac{dV^A_{\mu\nu}}{d\mathbf{R}_A} =
\underbrace{\frac{\partial V^A_{\mu\nu}}{\partial \mathbf{R}_\mu}\Big|_{\mu \in A}}_{\text{shell-drag, if any}}
+ \underbrace{\frac{\partial V^A_{\mu\nu}}{\partial \mathbf{R}_\nu}\Big|_{\nu \in A}}_{\text{shell-drag, if any}}
+ \underbrace{\left(-\frac{\partial V^A_{\mu\nu}}{\partial \mathbf{R}_\mu}-\frac{\partial V^A_{\mu\nu}}{\partial \mathbf{R}_\nu}\right)}_{\partial V^A_{\mu\nu}/\partial \mathbf{R}_A\text{, substituted}}
```

The shell-drag term and the substituted term cancel exactly whenever that
shell *is* on `A`, and survive (sign-flipped) whenever it *isn't*:

```math
\frac{dV^A_{\mu\nu}}{d\mathbf{R}_A} = K_A(\mu,\nu)\Big|_{\mu \notin A} \;+\; K_A(\nu,\mu)^\top\Big|_{\nu \notin A}
```

**Add the two pieces.** Writing them side by side gives the closed form
`∇nuclear!` actually evaluates:

```math
\frac{\partial V_{\mu\nu}}{\partial \mathbf{R}_A} = f(\mu) + f(\nu)^\top,
\qquad
f(p) = \begin{cases} +K_A(p,q) & p \notin A \\ -K_{\neg A}(p,q) & p \in A \end{cases}
```

which needs exactly two libcint calls per shell pair (`K_A` and `K_{¬A}`,
each also computed with `p,q` swapped and transposed for the `f(ν)`
piece) -- never three, and never a dedicated operator-derivative kernel.

Functions follow the same pattern found in [One-Electron Integrals](@ref), except you must include a mandatory argument `A` which indicates the atom for which derivatives are evaluated.
```julia-repl
julia> dr = ∇nuclear(bset, 2, 1, 2) # Second atom, first shell, second shell
julia> dropdims(dr, dims=(1,2)) # Shows derivatives along x, y, and z
3-element Vector{Float64}:
 0.9707070249006159
 0.0
 0.0
```

```@docs
∇nuclear
∇nuclear!
```

## Two-electron four centers

```@docs
∇ERI_2e4c(::BasisSet, ::Any)
∇ERI_2e4c!(::Any, ::BasisSet, ::Any)
∇sparseERI_2e4c
```

## Two-electron three centers

```@docs
∇ERI_2e3c(::BasisSet, ::BasisSet, ::Any)
∇ERI_2e3c!
```

## Two-electron two centers

```@docs
∇ERI_2e2c(::BasisSet, ::Any)
∇ERI_2e2c!(::Any, ::BasisSet, ::Any)
```

## Choosing a level

Every gradient here comes in the same four levels, from most convenient to
fastest. All of them ultimately go through the same libcint primitive, so
they agree to the last bit -- they differ only in how much bookkeeping they
do for you.

| level | example | notes |
|:------|:--------|:------|
| dense, allocating | `∇overlap(bs, A)` | returns a fresh `nbas × nbas × 3` array |
| dense, preallocated | `∇overlap!(out, bs, A)` | reuses `out`; zeroes it for you |
| shell pair/quartet | `∇overlap!(out, bs, A, i, j)` | validates sizes, resolves shell membership, returns the free zero when the block vanishes |
| bare primitive | `∇overlap_μ!(out, bs, i, j)` | one libcint call, no checks at all |

Two things are worth knowing before dropping a level.

**The dense drivers exploit symmetry that a per-pair loop must reproduce
itself.** `S`, `T`, `V`, `(P|Q)` and `(μν|P)` are all symmetric under the
bra-shell swap, so the dense builders visit only `i <= j` and write the
mirror from the same block. A loop that calls the shell-pair form for every
ordered `(i,j)` gets the right answer at twice the libcint cost, with
nothing to warn it. Loop `i <= j` and mirror, and the shell-pair form costs
essentially nothing over the dense driver (measured ~1.0x).

**Hoistable state should be hoisted.** Several routines rebuild per-call
state that depends only on the basis and the atom. Pass it in and the call
becomes allocation-free:

| routine | keyword | build it with |
|:--------|:--------|:--------------|
| `∇overlap!`, `∇kinetic!`, `∇nuclear!` (shell pair) | `scratch` | `Vector{Cdouble}(undef, 3*Nmax^2)` |
| `∇nuclear!` (shell pair) | `charges` | [`nuclear_charge_sets`](@ref) |
| `∇ERI_2e3c!` | `Bmerged` | merged `BasisSet` of the two bases |
| `∇sparseERI_2e4c` | `ij_vals`, `σvals` | [`schwarz_bounds`](@ref) |

### Shell-pair/shell-quartet-level forms

For callers building an integral-direct gradient or CPHF loop that never
wants the full dense array materialized:

```@docs
∇ERI_2e4c(::BasisSet, ::Int, ::Int, ::Int, ::Int, ::Int)
∇ERI_2e4c!(::Any, ::BasisSet, ::NTuple{4,Bool}, ::Int, ::Int, ::Int, ::Int, ::Vector{Float64})
```

### Bare libcint primitives

The lowest level. Each does exactly one libcint call and applies the sign
flip to the nuclear-coordinate convention -- no bounds checking, no
zero-block skipping, no output-size validation. Every one of them
differentiates its **first** shell argument, so any other center is reached
by passing that shell first and permuting the result yourself.

```@docs
∇overlap_μ!
∇kinetic_μ!
∇nuclear_μ!
∇ERI_2e2c_μ!
∇ERI_2e3c_μ!
∇ERI_2e4c_μ!
nuclear_charge_sets
schwarz_bounds
```

## Example

```julia-repl
julia> using GaussianBasis, Molecules

julia> water = Molecules.parse_string("""
       O        0.000000000000     -0.143225816552      0.000000000000
       H        1.638036840407      1.136548822547     -0.000000000000
       H       -1.638036840407      1.136548822547     -0.000000000000
       """);

julia> bset = BasisSet("sto-3g", water);

julia> gS = ∇overlap(bset, 2);  # ∂S/∂R for the first H atom

julia> size(gS)
(7, 7, 3)

julia> gS[2, 6, :]   # ⟨2s_O| shell | 1s_H(2)⟩ block, (x, y, z) components
3-element Vector{Float64}:
 -0.05625663685821103
 -0.0439525017729466
  0.0

julia> gERI = ∇ERI_2e4c(bset, 2);

julia> size(gERI)
(7, 7, 7, 7, 3)
```

For repeated calls, the mutating `∇overlap!`, `∇kinetic!`, `∇nuclear!`,
`∇ERI_2e4c!`, `∇ERI_2e3c!`, and `∇ERI_2e2c!` forms write into a
preallocated output array instead of reallocating.
