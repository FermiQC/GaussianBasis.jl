# --- Shell-quartet-level dense 4-center ERI Hessian ---
#
# GaussianBasis-idiomatic light bridge to libcint's second-derivative 2e
# kernels, mirroring the shell-quartet gradient primitive in
# `Gradients/TwoElectronGrad.jl` one derivative order up: no Schwarz
# screening here (that's Fermi's own loop's job, same split the gradient
# primitive already uses -- see that file's header comment), only the free,
# exact (not threshold-based) structural zero checks. Everything below this
# comment (the placement algebra below) is ported
# near-verbatim from Fermi.jl's `Hessians/TwoElectronHess.jl`, which already
# did the hard, error-prone work of confirming libcint's kernel-to-shell-pair
# mapping and axis orientation against `cint_funcs.h` and finite difference --
# moved here rather than re-derived, and no longer needed at the Fermi.jl
# call site once Fermi switches to calling this instead.
#
# Which kernel gives which shell-pair's second derivative is NOT what the
# "ip1", "ip2" naming suggests at first glance:
#
#   cint2e_ipip1  "(NABLA NABLA i j|R12|k l)"   both derivatives on shell i
#                                                (same shell, position 1)
#   cint2e_ipvip1 "(NABLA i NABLA j|R12|k l)"    one derivative each on shells
#                                                i,j (SAME electron pair --
#                                                i,j both "bra", positions 1,2)
#   cint2e_ip1ip2 "(NABLA i j|R12|NABLA k l)"    one derivative each on shells
#                                                i,k (OPPOSITE electron pairs
#                                                -- i is "bra" position 1, k
#                                                is "ket" position 3)
#
# ip1ip2 giving the opposite-pair (i.e. i vs k, not i vs j) cross derivative
# directly is the crucial piece -- with it, all 10 distinct shell-pair second
# derivatives of a quartet (same-shell x4, same-pair-cross x2, opposite-pair-
# cross x4) are directly available via these 3 kernels plus shell-argument
# reordering (exploiting the standard 8-fold ERI permutation symmetry
# (pq|rs)=(qp|rs)=(pq|sr)=(rs|pq)=...).
#
# Every raw kernel whose two derivatives land on DIFFERENT shells (ipvip1,
# ip1ip2) has its two derivative-component axes in (position-2-or-4,
# position-1-or-3) order in the raw buffer -- reversed from the naive
# expectation, same pattern found for the one-electron mixed kernels
# (cint1e_iprinvip_sph!). ipip1 (same shell, mixed partials of the same
# point commute) needs no such correction.

"""
    ∇2ERI_2e4c_μμ!(out, BS::BasisSet{LCint}, i, j, k, l)
    ∇2ERI_2e4c_μν!(out, BS::BasisSet{LCint}, i, j, k, l)
    ∇2ERI_2e4c_μλ!(out, BS::BasisSet{LCint}, i, j, k, l)

Bare libcint calls for the 4-center `(ij|kl)` Hessian. Three kernels cover
every placement, because the four centers form two bra/ket pairs and a second
derivative can land within one pair or across them:

| primitive | libcint | placement |
|:----------|:--------|:----------|
| `_μμ!` | `ipip1` | both derivatives on shell-argument position 1 |
| `_μν!` | `ipvip1` | one on position 1, one on position 2 (same pair) |
| `_μλ!` | `ip1ip2` | one on position 1, one on position 3 (across pairs) |

Every other placement is reached by permuting the shell arguments -- e.g.
both derivatives on the third shell is `_μμ!(out, BS, k, l, i, j)` -- and
transposing the AO axes of the result back. `∇2ERI_2e4c!` does exactly that,
folding each permutation into the index expression of its accumulate.

!!! warning "Derivative-axis order"
    All three emit their two derivative-component axes in **(q,p)** order,
    reversed from the naive expectation. For `_μμ!` this is invisible (mixed
    partials of the same point commute), but not for the cross kernels.

`out` is a raw `(Na,Nb,Nc,Nd,3,3)` block and may be a contiguous view.

> No bounds checking, no zero-block skipping, no output-size validation.
"""
function ∇2ERI_2e4c_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int, l::Int)
    lib = BS.lib
    cint2e_ipip1_sph!(out, @SVector(Cint[i-1, j-1, k-1, l-1]),
                      lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2ERI_2e4c_μμ!)
function ∇2ERI_2e4c_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int, l::Int)
    lib = BS.lib
    cint2e_ipvip1_sph!(out, @SVector(Cint[i-1, j-1, k-1, l-1]),
                       lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2ERI_2e4c_μμ!)
function ∇2ERI_2e4c_μλ!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int, l::Int)
    lib = BS.lib
    cint2e_ip1ip2_sph!(out, @SVector(Cint[i-1, j-1, k-1, l-1]),
                       lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

# Accumulate a raw libcint block into `out`, folding both permutations into
# the index expression rather than materializing a transposed copy.
#
#   W      compile-time 4-tuple: raw axis n carries block axis W[n]. So a
#          kernel called with shells (q,p,r,s) has W = (2,1,3,4).
#   dswap  final derivative-axis orientation. libcint always emits (q,p);
#          `false` corrects it to (p,q), `true` leaves it (used when
#          posA > posB, where mixed partials commute).
@inline function _acc4c!(out, bufv, N::NTuple{4,Int}, ::Val{W}, dswap::Bool) where W
    D1 = N[W[1]]; D2 = N[W[2]]; D3 = N[W[3]]; D4 = N[W[4]]
    n = D1*D2*D3*D4
    s2 = D1; s3 = D1*D2; s4 = D1*D2*D3
    @inbounds for q = 1:3, p = 1:3
        x = dswap ? p : q
        y = dswap ? q : p
        od = n*(x-1) + 3n*(y-1)
        for d = 1:N[4], c = 1:N[3], b = 1:N[2], a = 1:N[1]
            t = (a, b, c, d)
            li = t[W[1]] + s2*(t[W[2]]-1) + s3*(t[W[3]]-1) + s4*(t[W[4]]-1)
            out[a,b,c,d,p,q] += bufv[od + li]
        end
    end
    return out
end

# One (posA,posB) placement accumulated into `out`. Positions 1..4 map to the
# four shell arguments (i,j,k,l). Replaces the former
# former same/cross helper pair *and* their permutedims!-based scratch twins:
# one implementation now serves both the allocating and preallocated entry
# points, so the placement algebra exists in exactly one place.
@inline function _∇24c_place!(out, bufv, BS, i, j, k, l, N::NTuple{4,Int},
                              posA::Int, posB::Int)
    if posA == posB
        if posA == 1
            ∇2ERI_2e4c_μμ!(bufv, BS, i, j, k, l); _acc4c!(out, bufv, N, Val((1,2,3,4)), false)
        elseif posA == 2
            ∇2ERI_2e4c_μμ!(bufv, BS, j, i, k, l); _acc4c!(out, bufv, N, Val((2,1,3,4)), false)
        elseif posA == 3
            ∇2ERI_2e4c_μμ!(bufv, BS, k, l, i, j); _acc4c!(out, bufv, N, Val((3,4,1,2)), false)
        else
            ∇2ERI_2e4c_μμ!(bufv, BS, l, k, i, j); _acc4c!(out, bufv, N, Val((4,3,1,2)), false)
        end
    else
        swp = posA > posB
        a, b = swp ? (posB, posA) : (posA, posB)
        if (a, b) == (1, 2)
            ∇2ERI_2e4c_μν!(bufv, BS, i, j, k, l); _acc4c!(out, bufv, N, Val((1,2,3,4)), swp)
        elseif (a, b) == (3, 4)
            ∇2ERI_2e4c_μν!(bufv, BS, k, l, i, j); _acc4c!(out, bufv, N, Val((3,4,1,2)), swp)
        elseif (a, b) == (1, 3)
            ∇2ERI_2e4c_μλ!(bufv, BS, i, j, k, l); _acc4c!(out, bufv, N, Val((1,2,3,4)), swp)
        elseif (a, b) == (1, 4)
            ∇2ERI_2e4c_μλ!(bufv, BS, i, j, l, k); _acc4c!(out, bufv, N, Val((1,2,4,3)), swp)
        elseif (a, b) == (2, 3)
            ∇2ERI_2e4c_μλ!(bufv, BS, j, i, k, l); _acc4c!(out, bufv, N, Val((2,1,3,4)), swp)
        else # (2,4)
            ∇2ERI_2e4c_μλ!(bufv, BS, j, i, l, k); _acc4c!(out, bufv, N, Val((2,1,4,3)), swp)
        end
    end
    return out
end

"""
    ∇2ERI_2e4c(BS::BasisSet, iA::Int, iB::Int, i::Int, j::Int, k::Int, l::Int)
    ∇2ERI_2e4c(BS::BasisSet, Xflag::NTuple{4,Bool}, Yflag::NTuple{4,Bool}, i::Int, j::Int, k::Int, l::Int)

Shell-quartet-level dense 4-center ERI Hessian: `∂²(ij|kl)/∂R_iA∂R_iB` for
shells `i,j,k,l` w.r.t. atoms `iA,iB`'s three Cartesian directions each,
with `R_iA`/`R_iB` in bohr (see [Hessians](@ref) for units), as an
`(Ni,Nj,Nk,Nl,3,3)` block. No Schwarz screening either way -- callers
wanting that should screen before calling.

Two forms, same relationship as `∇ERI_2e4c`'s: the `iA::Int,iB::Int` form
computes `Xflag`/`Yflag` (which of `i,j,k,l` sit on `iA`/`iB`) and returns
the free zero (no libcint call) whenever no shell touches `iA`, no shell
touches `iB`, or `iA==iB` with all four shells on that one atom
(translational invariance). The `Xflag,Yflag::NTuple{4,Bool}` form is the
unchecked core: no membership computation, no zero-skip -- for callers
that have already screened at a higher level and precomputed the flags
once outside a hot loop.
"""
function ∇2ERI_2e4c(BS::BasisSet, iA::Int, iB::Int, i::Int, j::Int, k::Int, l::Int)
    Xflag = on_atom_flags(BS, iA, i, j, k, l)
    Yflag = on_atom_flags(BS, iB, i, j, k, l)
    Ni = num_basis(BS.shells[i]); Nj = num_basis(BS.shells[j])
    Nk = num_basis(BS.shells[k]); Nl = num_basis(BS.shells[l])
    out = zeros(Ni, Nj, Nk, Nl, 3, 3)
    (!any(Xflag) || !any(Yflag) || (iA == iB && all(Xflag))) && return out
    return ∇2ERI_2e4c!(out, BS, Xflag, Yflag, i, j, k, l)
end

function ∇2ERI_2e4c!(out, BS::BasisSet, iA::Int, iB::Int, i::Int, j::Int, k::Int, l::Int)
    Xflag = on_atom_flags(BS, iA, i, j, k, l)
    Yflag = on_atom_flags(BS, iB, i, j, k, l)
    if !any(Xflag) || !any(Yflag) || (iA == iB && all(Xflag))
        out .= 0.0
        return out
    end
    return ∇2ERI_2e4c!(out, BS, Xflag, Yflag, i, j, k, l)
end

function ∇2ERI_2e4c(BS::BasisSet, Xflag::NTuple{4,Bool}, Yflag::NTuple{4,Bool}, i::Int, j::Int, k::Int, l::Int)
    Ni = num_basis(BS.shells[i]); Nj = num_basis(BS.shells[j])
    Nk = num_basis(BS.shells[k]); Nl = num_basis(BS.shells[l])
    out = zeros(Ni, Nj, Nk, Nl, 3, 3)
    return ∇2ERI_2e4c!(out, BS, Xflag, Yflag, i, j, k, l)
end

function ∇2ERI_2e4c!(out, BS::BasisSet, Xflag::NTuple{4,Bool}, Yflag::NTuple{4,Bool}, i::Int, j::Int, k::Int, l::Int)
    Ni = num_basis(BS.shells[i]); Nj = num_basis(BS.shells[j])
    Nk = num_basis(BS.shells[k]); Nl = num_basis(BS.shells[l])
    buf = Vector{Cdouble}(undef, 9 * Ni * Nj * Nk * Nl)
    return ∇2ERI_2e4c!(out, BS, Xflag, Yflag, i, j, k, l, buf)
end

"""
    ∇2ERI_2e4c!(out, BS::BasisSet, Xflag, Yflag, i, j, k, l, buf)

Scratch-buffer-accepting core: `buf` (sized `>= 9*Nmax^4`, `Nmax` = the
largest `num_basis` over any shell the caller will ever pass) is
caller-owned and reused across every call instead of allocated fresh --
zero-allocation, for callers in a hot per-quartet loop. Not thread-safe to
share: each concurrent caller (e.g. each worker task) needs its own `buf`.
"""
function ∇2ERI_2e4c!(out, BS::BasisSet, Xflag::NTuple{4,Bool}, Yflag::NTuple{4,Bool}, i::Int, j::Int, k::Int, l::Int,
                      buf::Vector{Cdouble})
    # Up to 16 placements per shell quartet (4x4 posA/posB, against the
    # gradient's 4 branches), so anything allocated per placement is
    # multiplied accordingly -- this sits in Fermi.jl's Hessian inner loop.
    Ni = num_basis(BS.shells[i]); Nj = num_basis(BS.shells[j])
    Nk = num_basis(BS.shells[k]); Nl = num_basis(BS.shells[l])
    if size(out) != (Ni, Nj, Nk, Nl, 3, 3)
        throw(DimensionMismatch("Size of the output array needs to be ($Ni, $Nj, $Nk, $Nl, 3, 3)."))
    end

    fill!(out, 0.0)
    N = (Ni, Nj, Nk, Nl)
    # `buf` is passed whole rather than as a `view(buf, 1:n)`: libcint writes
    # only the first n entries and `_acc4c!` indexes linearly, so the view
    # bought nothing and cost a SubArray allocation per call.
    for posA in 1:4
        Xflag[posA] || continue
        for posB in 1:4
            Yflag[posB] || continue
            _∇24c_place!(out, buf, BS, i, j, k, l, N, posA, posB)
        end
    end

    return out
end
