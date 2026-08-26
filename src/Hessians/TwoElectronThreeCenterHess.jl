# Second derivative of the 3-center two-electron integral (μν|P) w.r.t. two
# nuclear positions (atoms iA, iB). Unlike overlap/kinetic (2 shell
# positions) or the 4-center ERI (4 shell positions but bra/ket-paired),
# this integral has 3 shell positions with no pairing symmetry between them:
# μ, ν (both "regular" basis shells, symmetric under μ<->ν swap since
# (μν|P)=(νμ|P)) and P (the auxiliary/fitting shell). That gives 3
# same-shell placements (μμ, νν, PP) and 3 distinct cross placements (μν,
# μP, νP) -- 6 total instead of the 2-shell case's 2 (same, cross).
#
# Kernel-to-placement mapping (confirmed against libcint's own cint_funcs.h
# doc comments, mirroring exactly how Fermi.jl's 4-center ERI Hessian
# (TwoElectronHess.jl) already worked this out for cint2e_ipvip1/ip1ip2 --
# same naming logic, one derivative order down in shell count):
#
#   cint3c2e_ipip1  both derivatives on shell-argument-position-1 (μ, by
#                   default -- shell-argument order can be permuted to put
#                   ν there instead, same trick as the first-derivative
#                   ∇ERI_2e3c! already uses)
#   cint3c2e_ipip2  both derivatives on the auxiliary shell (libcint's own
#                   "position 2" for 3c2e, matching how the existing
#                   first-derivative cint3c2e_ip1/ip2 already split
#                   "regular" vs "auxiliary" rather than literal shell-arg
#                   position 1 vs 2)
#   cint3c2e_ipvip1 one derivative each on μ and ν -- the "same group"
#                   cross term (both regular shells, no aux involved)
#   cint3c2e_ip1ip2 one derivative each on a regular shell and the
#                   auxiliary shell -- the "opposite group" cross term
#
# Same as the 4-center case's ipvip1/ip1ip2, the two derivative-component
# axes of a CROSS kernel's raw output come out reversed from the naive
# (position-1-component, position-2-component) expectation -- ipip-type
# (same-shell) kernels need no such correction, since mixed partials of the
# same point commute. The primitives below hand back libcint's raw (q,p)
# axis order unchanged -- see their docstring -- and `_acc3c!` applies the
# correction as it accumulates, so the placement logic never has to
# re-derive it per call site.
#
# Validated: ∇2FD_ERI_2e3c (central difference of the already-trusted
# ∇ERI_2e3c), translational invariance (sum over all iB of the Hessian
# block, for fixed iA, is zero -- checked to machine precision), and an
# independent cross-check of the μν cross term against the value implied by
# translational invariance applied to the already-validated μμ/νν/μP/νP
# terms (d2f/dμdν = -d2f/dμ² - d2f/dμdP, from differentiating the
# first-derivative translational-invariance identity dμ+dν+dP=0 once more
# w.r.t. μ) -- an independent check on cint3c2e_ipvip1 specifically that
# doesn't depend on finite difference at all.

"""
    ∇2ERI_2e3c_μμ!(out, BS::BasisSet{LCint}, i, j, k)
    ∇2ERI_2e3c_PP!(out, BS::BasisSet{LCint}, i, j, k)
    ∇2ERI_2e3c_μν!(out, BS::BasisSet{LCint}, i, j, k)
    ∇2ERI_2e3c_μP!(out, BS::BasisSet{LCint}, i, j, k)

Bare libcint calls for the 3-center `(μν|P)` Hessian. `BS` is the **merged**
basis (regular shells then auxiliary), so an auxiliary shell `k` is addressed
as `k + BS1.nshells` -- same convention as `∇ERI_2e3c_μ!`.

Four kernels are needed because the three centers fall into two groups (the
two bra shells, and the lone ket/auxiliary shell), so a second derivative can
land in four distinct ways:

| primitive | libcint | both/one derivative on |
|:----------|:--------|:-----------------------|
| `_μμ!` | `ipip1` | both on the FIRST shell argument |
| `_PP!` | `ipip2` | both on the auxiliary shell |
| `_μν!` | `ipvip1` | one on each of the two bra shells |
| `_μP!` | `ip1ip2` | one on the first bra shell, one on the auxiliary |

`ν` is reached by passing it first (`_μμ!(out, BS, j, i, k)` etc.) and
transposing the AO axes back.

!!! warning "Derivative-axis order"
    All four return their two derivative-component axes in **(q,p)** order,
    i.e. `out[a,b,c,q,p]` is `∂²/∂(first)_p ∂(second)_q`. Callers must
    orient this deliberately; `∇2ERI_2e3c!` folds the correction into its
    scatter.

`out` is a raw `(Na,Nb,Nc,3,3)` block and may be a contiguous view.

> No bounds checking, no zero-block skipping, no output-size validation.
"""
function ∇2ERI_2e3c_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int)
    lib = BS.lib
    cint3c2e_ipip1_sph!(out, @SVector(Cint[i-1, j-1, k-1]),
                        lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2ERI_2e3c_μμ!)
function ∇2ERI_2e3c_PP!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int)
    lib = BS.lib
    cint3c2e_ipip2_sph!(out, @SVector(Cint[i-1, j-1, k-1]),
                        lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2ERI_2e3c_μμ!)
function ∇2ERI_2e3c_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int)
    lib = BS.lib
    cint3c2e_ipvip1_sph!(out, @SVector(Cint[i-1, j-1, k-1]),
                         lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2ERI_2e3c_μμ!)
function ∇2ERI_2e3c_μP!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int)
    lib = BS.lib
    cint3c2e_ip1ip2_sph!(out, @SVector(Cint[i-1, j-1, k-1]),
                         lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

# Accumulate one raw libcint block into `block`, folding both possible
# permutations into the index expression:
#   aoswap : the kernel was called with the two bra shells swapped, so the raw
#            block is (Nj,Ni,Nk,...) and its first two AO axes transpose back
#   dswap  : final derivative-axis orientation. libcint always emits (q,p);
#            `false` corrects that to (p,q), `true` leaves it (used when
#            posA > posB, where mixed partials commute).
@inline function _acc3c!(block, bufv, Ni, Nj, Nk, aoswap::Bool, dswap::Bool)
    Nijk = Ni*Nj*Nk
    @inbounds for q = 1:3, p = 1:3
        x = dswap ? p : q
        y = dswap ? q : p
        od = Nijk*(x-1) + 3*Nijk*(y-1)
        for c = 1:Nk, b = 1:Nj, a = 1:Ni
            idx = aoswap ? (b + Nj*(a-1) + Nj*Ni*(c-1)) : (a + Ni*(b-1) + Ni*Nj*(c-1))
            block[a,b,c,p,q] += bufv[od + idx]
        end
    end
    return block
end

# One (posA,posB) placement: pick the kernel, the shell order and the two
# permutation flags, evaluate, and accumulate. Positions are 1=μ, 2=ν, 3=aux.
@inline function _∇23c_place!(block, bufv, Bm, i, j, kshell, Ni, Nj, Nk, posA::Int, posB::Int)
    if posA == posB
        if posA == 1
            ∇2ERI_2e3c_μμ!(bufv, Bm, i, j, kshell);  _acc3c!(block, bufv, Ni,Nj,Nk, false, false)
        elseif posA == 2
            ∇2ERI_2e3c_μμ!(bufv, Bm, j, i, kshell);  _acc3c!(block, bufv, Ni,Nj,Nk, true,  false)
        else
            ∇2ERI_2e3c_PP!(bufv, Bm, i, j, kshell);  _acc3c!(block, bufv, Ni,Nj,Nk, false, false)
        end
    else
        swp = posA > posB
        a, b = swp ? (posB, posA) : (posA, posB)
        if (a, b) == (1, 2)
            ∇2ERI_2e3c_μν!(bufv, Bm, i, j, kshell);  _acc3c!(block, bufv, Ni,Nj,Nk, false, swp)
        elseif (a, b) == (1, 3)
            ∇2ERI_2e3c_μP!(bufv, Bm, i, j, kshell);  _acc3c!(block, bufv, Ni,Nj,Nk, false, swp)
        else # (2,3)
            ∇2ERI_2e3c_μP!(bufv, Bm, j, i, kshell);  _acc3c!(block, bufv, Ni,Nj,Nk, true,  swp)
        end
    end
    return block
end

"""
    ∇2ERI_2e3c(BS1::BasisSet, BS2::BasisSet, iA, iB)

Second derivative (Hessian, atoms `iA`,`iB`, with `R_iA`/`R_iB` in bohr --
see [Hessians](@ref) for units) of the 3-center two-electron integral
`(μν|P)` (`BS1`=regular basis, `BS2`=auxiliary/fitting basis -- density
fitting). Output `(BS1.nbas,BS1.nbas,BS2.nbas,3,3)`. See this file's header
comment for the shell-position combinatorics and kernel mapping.
"""
function ∇2ERI_2e3c(BS1::BasisSet, BS2::BasisSet, iA, iB)
    out = zeros(BS1.nbas, BS1.nbas, BS2.nbas, 3, 3)
    return ∇2ERI_2e3c!(out, BS1, BS2, iA, iB)
end

function ∇2ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet, iA, iB; Bmerged::Union{Nothing,BasisSet}=nothing)

    if size(out) != (BS1.nbas, BS1.nbas, BS2.nbas, 3, 3)
        throw(DimensionMismatch("Size of the output array needs to be (N1, N1, N2, 3, 3)."))
    end
    out .= 0.0

    # Bmerged depends only on BS1/BS2, never on iA/iB -- callers looping
    # over atom pairs (e.g. Fermi.jl's DF-Hessian, O(natm^2) calls) can
    # build it once and pass it in via this keyword. See ∇ERI_2e3c!'s
    # identical pattern (Gradients/TwoElectronGrad.jl).
    if Bmerged === nothing
        Bmerged = merge_basis(BS1, BS2)
    end

    Aat = BS1.atoms[iA]
    Bat = BS1.atoms[iB]

    Nvals1 = num_basis.(BS1.shells)
    ao_offset1 = cumsum(Nvals1) .- Nvals1

    Nvals2 = num_basis.(BS2.shells)
    ao_offset2 = cumsum(Nvals2) .- Nvals2

    Nmax = 9*maximum(Nvals1)^2*maximum(Nvals2)
    buf = Vector{Cdouble}(undef, Nmax)
    blkbuf = Vector{Cdouble}(undef, Nmax)

    for i = 1:BS1.nshells
        atom_i = BS1.shells[i].atom
        for j = i:BS1.nshells
            atom_j = BS1.shells[j].atom
            for k = 1:BS2.nshells
                atom_k = BS2.shells[k].atom

                Xflag = (atom_i == Aat, atom_j == Aat, atom_k == Aat)
                Yflag = (atom_i == Bat, atom_j == Bat, atom_k == Bat)
                (any(Xflag) && any(Yflag)) || continue

                Ni = Nvals1[i]
                Nj = Nvals1[j]
                Nk = Nvals2[k]

                ioff = ao_offset1[i]
                joff = ao_offset1[j]
                koff = ao_offset2[k]

                I = (ioff+1):(ioff+Ni)
                J = (joff+1):(joff+Nj)
                K = (koff+1):(koff+Nk)
                kshell = k + BS1.nshells

                blk = reshape(view(blkbuf, 1:9*Ni*Nj*Nk), Ni, Nj, Nk, 3, 3)
                fill!(blk, 0.0)
                bufv = view(buf, 1:9*Ni*Nj*Nk)
                for posA in 1:3, posB in 1:3
                    (Xflag[posA] && Yflag[posB]) || continue
                    _∇23c_place!(blk, bufv, Bmerged, i, j, kshell, Ni, Nj, Nk, posA, posB)
                end

                # (μν|P) is symmetric under μ<->ν, so mirror the block; the
                # transpose is just the swapped index expression.
                @inbounds for q = 1:3, p = 1:3, c = 1:Nk, b = 1:Nj, a = 1:Ni
                    v = blk[a,b,c,p,q]
                    out[ioff+a, joff+b, koff+c, p, q] += v
                    i != j && (out[joff+b, ioff+a, koff+c, p, q] = v)
                end
            end
        end
    end

    return out
end


"""
    ∇2ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet, iA, iB, i::Int, j::Int, k::Int;
                scratch=nothing, Bmerged=nothing)
    ∇2ERI_2e3c(BS1::BasisSet, BS2::BasisSet, iA, iB, i::Int, j::Int, k::Int)

Shell-triple Hessian block of `(μν|P)`: `∂²(ij|k)/∂R_iA∂R_iB` for regular
shells `i,j` of `BS1` and auxiliary shell `k` of `BS2` (shell indices, not AO
indices). `out` must be `(Ni,Nj,Nk,3,3)` and is overwritten, so a reused
buffer is safe.

This is the level integral-direct and CPHF code wants: the dense form is a
`BS1.nbas² × BS2.nbas × 3 × 3` tensor, which is rarely worth materializing.

Returns the free zero -- no libcint call -- unless at least one of the three
shells sits on `iA` *and* at least one sits on `iB`.

Pass `scratch` (at least `9*Ni*Nj*Nk` elements) and `Bmerged` (from the two
bases; it depends only on them, never on the atoms or shells) to make the
call allocation-free in a loop.
"""
function ∇2ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet, iA, iB,
                     i::Int, j::Int, k::Int;
                     scratch=nothing, Bmerged::Union{Nothing,BasisSet}=nothing)
    Ni = num_basis(BS1.shells[i]); Nj = num_basis(BS1.shells[j])
    Nk = num_basis(BS2.shells[k])
    if size(out) != (Ni, Nj, Nk, 3, 3)
        throw(DimensionMismatch("Size of the output array needs to be ($Ni, $Nj, $Nk, 3, 3)"))
    end
    fill!(out, 0.0)

    Aat = BS1.atoms[iA]; Bat = BS1.atoms[iB]
    ai = BS1.shells[i].atom; aj = BS1.shells[j].atom; ak = BS2.shells[k].atom
    Xflag = (ai == Aat, aj == Aat, ak == Aat)
    Yflag = (ai == Bat, aj == Bat, ak == Bat)
    (any(Xflag) && any(Yflag)) || return out

    Bm = Bmerged === nothing ? merge_basis(BS1, BS2) : Bmerged
    n = 9*Ni*Nj*Nk
    buf = scratch === nothing ? Vector{Cdouble}(undef, n) : scratch
    length(buf) >= n ||
        throw(DimensionMismatch("scratch must hold at least $n elements"))
    bufv = view(buf, 1:n)

    kshell = k + BS1.nshells
    for posA in 1:3, posB in 1:3
        (Xflag[posA] && Yflag[posB]) || continue
        _∇23c_place!(out, bufv, Bm, i, j, kshell, Ni, Nj, Nk, posA, posB)
    end
    return out
end

function ∇2ERI_2e3c(BS1::BasisSet, BS2::BasisSet, iA, iB, i::Int, j::Int, k::Int)
    out = zeros(num_basis(BS1.shells[i]), num_basis(BS1.shells[j]),
                num_basis(BS2.shells[k]), 3, 3)
    return ∇2ERI_2e3c!(out, BS1, BS2, iA, iB, i, j, k)
end
