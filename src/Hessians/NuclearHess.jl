# Nuclear attraction Hessian, ∂²(sum_C -Z_C/|r-R_C|)_ij / ∂A∂B.
#
# Unlike overlap/kinetic, this integral depends on THREE kinds of positions:
# shell i's center ("bra"), shell j's center ("ket"), and the position of
# whichever nucleus C is providing the potential ("op"). Differentiating
# twice w.r.t. two *atoms* A, B therefore has two parts:
#
#   1. Both derivatives land on shell centers (bra/ket only, "op" untouched).
#      This part doesn't care which nucleus supplies the potential, so it's
#      just the whole-molecule, all-nuclei-summed kernels ipipnuc/ipnucip --
#      handled by ∇21e!(..., "nuclear_shellonly", ...), identical in
#      structure to the overlap/kinetic case.
#
#   2. At least one derivative lands on a nuclear charge position. A
#      derivative w.r.t. atom A can only hit nucleus C's "op" role if C==A
#      (moving a nucleus that ISN'T A doesn't change anything when we
#      differentiate w.r.t. A). So this part only ever needs nucleus C = A
#      or C = B, isolated one at a time via the repositionable rinv operator
#      (env[PTR_RINV_ORIG]), each contribution scaled by that nucleus's
#      charge. This is the piece the existing gradient code gets via
#      charge-masking (lc_atoms_A/lc_atoms_woA); the rinv-repositioning
#      approach isolates a single nucleus's contribution directly instead,
#      which is much less bookkeeping one derivative order up.
#
# For a fixed shell pair (i,j) and isolated nucleus C (charge Z_C, position
# R_C), translational invariance (bra + ket + op directions sum to zero,
# component-wise, at every derivative order since it holds identically as a
# function) gives the second-derivative building blocks in terms of three
# DIRECTLY available kernels -- ipiprinv (both derivatives on shell i,
# symmetric in its own two derivative-component axes since both land on the
# same point), iprinvip (one derivative on each shell -- *not* symmetric in
# its own two axes, since bra/ket/op are three independent points once a
# third center (R_C) is in the picture; its raw buffer layout is also
# (ket,bra), not (bra,ket) -- both facts verified against direct finite
# differences of cint1e_iprinv_sph!), and ipiprinv with shells swapped and
# the AO axes transposed back (both derivatives on shell j). Writing
# bk[p,q] := d2(rinv_ij)/dbra_p dket_q (i.e. already corrected to (bra,ket)
# axis order) and bkT := its (3,3)-axis transpose:
#
#   d2/dbra_p dbra_q = -Z_C * ipiprinv(i,j)[p,q]                      (symmetric in p,q)
#   d2/dbra_p dket_q = -Z_C * bk[p,q]
#   d2/dket_p dket_q = -Z_C * ipiprinv(j,i)^T[p,q]                    (symmetric in p,q)
#   d2/dbra_p dop_q  = +Z_C * (bb + bk)[p,q]
#   d2/dop_p  dbra_q = +Z_C * (bb + bkT)[p,q]      (= d2/dbra_q dop_p, mixed partials commute)
#   d2/dket_p dop_q  = +Z_C * (bkT + kk)[p,q]
#   d2/dop_p  dket_q = +Z_C * (bk + kk)[p,q]       (= d2/dket_q dop_p)
#   d2/dop_p  dop_q  = -Z_C * (bb + bk + bkT + kk)[p,q]               (symmetric in p,q)
#
# where bb=ipiprinv(i,j) and kk=ipiprinv(j,i) transposed-back, both raw
# (unscaled). Note d2/dbra dop and d2/dop dbra are NOT interchangeable
# lookups into the same tensor here (bk isn't symmetric under its own axis
# swap) -- each of the two placements below uses its own orientation. All of
# this was validated (sign and axis convention, at every step) against the
# whole-molecule kernels via cint1e_ipipnuc == -sum_C Z_C*cint1e_ipiprinv(R_C)
# (and the ipnucip/iprinvip analog) to machine precision, and the full
# assembly against finite differences of the existing, trusted ∇nuclear
# gradient to ~1e-9, plus an exact (machine-precision) translational
# invariance check (sum over all B of the full Hessian block, for fixed A,
# is zero). The dbra-dbra/dbra-dket/dket-dket combinations are exactly the
# whole-molecule piece handled by part 1 above (summed over ALL nuclei
# automatically) and must be excluded here to avoid double-counting -- only
# combinations touching "op" at least once are added below.

function rinv_env(BS::BasisSet, R)
    env = copy(BS.lib.env)
    env[5:7] .= R
    return env
end

"""
    ∇2nuclear_rinv_μμ!(out, BS::BasisSet{LCint}, env, i::Int, j::Int)
    ∇2nuclear_rinv_μν!(out, BS::BasisSet{LCint}, env, i::Int, j::Int)

Libcint calls for the `1/|r-C|` (single-nucleus) Hessian pieces of the
nuclear-attraction integral, for shells `i,j` of `BS`. `env` selects WHICH
nucleus: it is a copy of `BS.lib.env` with the rinv origin overwritten (see
`rinv_env`), the same way [`∇nuclear_μ!`](@ref) takes a charge-fudged `atm`
array to select which potential it differentiates.

`_μμ!` places both derivatives on shell `i`; `_μν!` places one on shell `i`
and one on shell `j`. Both return a raw `(Ni,Nj,3,3)` block and accept a
contiguous view.

!!! warning "Derivative-axis order"
    `_μν!` returns its two derivative axes in **(ket, bra)** order. Unlike
    the two-center `ipXip` kernels, a mixed bra/ket derivative here involves
    *three* independent points (bra center, ket center, rinv origin), so it
    has no reason to be symmetric under swapping those axes and the caller
    must orient it explicitly.

> No bounds checking, no zero-block skipping, no output-size validation.
"""
function ∇2nuclear_rinv_μμ!(out, BS::BasisSet{LCint}, env, i::Int, j::Int)
    lib = BS.lib
    cint1e_ipiprinv_sph!(out, @SVector(Cint[i-1, j-1]),
                         lib.atm, lib.natm, lib.bas, lib.nbas, env)
    return out
end

@doc (@doc ∇2nuclear_rinv_μμ!)
function ∇2nuclear_rinv_μν!(out, BS::BasisSet{LCint}, env, i::Int, j::Int)
    lib = BS.lib
    cint1e_iprinvip_sph!(out, @SVector(Cint[i-1, j-1]),
                         lib.atm, lib.natm, lib.bas, lib.nbas, env)
    return out
end

"""
    ∇2nuclear(BS::BasisSet, iA, iB) -> Array{Float64,4}

Second derivative (Hessian) of the AO nuclear attraction matrix `V` w.r.t.
atoms `iA`,`iB`'s three Cartesian coordinates each, `∂²V/∂R_iA∂R_iB`
(accounting for both shell-center and nuclear-charge-position derivatives),
with `R_iA`/`R_iB` in bohr (see [Hessians](@ref) for units). Returns a
dense `nbas × nbas × 3 × 3` array. For repeated calls, see `∇2nuclear!`.
"""
function ∇2nuclear(BS::BasisSet, iA, iB)
    out = zeros(BS.nbas, BS.nbas, 3, 3)
    return ∇2nuclear!(out, BS, iA, iB)
end

"""
    ∇2nuclear!(out, BS::BasisSet, iA, iB)

Mutating counterpart of [`∇2nuclear`](@ref): writes the dense
`nbas × nbas × 3 × 3` Hessian into `out` instead of allocating it. `out` is
zeroed first, so a reused buffer is safe.

Assembles both contributions -- the piece with both derivatives on shell
centers (summed over every nucleus) and the piece with at least one
derivative on a moving nuclear charge (isolated one nucleus at a time
through a repositionable point-charge operator, see
[`∇2nuclear_rinv_μμ!`](@ref)).
"""
function ∇2nuclear!(out, BS::BasisSet, iA, iB)

    if size(out) != (BS.nbas, BS.nbas, 3, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3, 3)"))
    end
    out .= 0.0

    # Part 1: both derivatives on shell centers, whole-molecule potential.
    ∇21e!(out, BS, "nuclear_shellonly", iA, iB)

    # Part 2: at least one derivative on a nuclear charge position, isolated
    # per relevant nucleus (only A and/or B can ever be the moving nucleus).
    Aat = BS.atoms[iA]
    Bat = BS.atoms[iB]

    Nvals = num_basis.(BS.shells)
    ao_offset = cumsum(Nvals) .- Nvals
    Nmax = maximum(Nvals)
    buf = zeros(Cdouble, 9*Nmax^2)
    # Reused across every (nucleus, shell pair) below instead of allocating
    # bb/bk/kk/bkT/block fresh each time -- same rationale as
    # OneElectronHess.jl's fix. bo/bo_T/ko/ko_T/oo (below) are handled via
    # fused broadcast (`block .+= Zc .* (bb .+ bk)`, etc.) instead of naming
    # an intermediate array, so no separate buffers are needed for those.
    bb_buf = zeros(Cdouble, Nmax, Nmax, 3, 3)
    bk_buf = zeros(Cdouble, Nmax, Nmax, 3, 3)
    kk_buf = zeros(Cdouble, Nmax, Nmax, 3, 3)
    bkT_buf = zeros(Cdouble, Nmax, Nmax, 3, 3)
    block_buf = zeros(Cdouble, Nmax, Nmax, 3, 3)

    nuclei = Aat == Bat ? (Aat,) : (Aat, Bat)

    @inbounds for c in nuclei
        Zc = c.Z
        Rc = collect(c.xyz) ./ Molecules.bohr_to_angstrom
        env_c = rinv_env(BS, Rc)
        Ca = (c == Aat)
        Cb = (c == Bat)

        for i in 1:BS.nshells
            atom_i = BS.shells[i].atom
            X_i = (atom_i == Aat)
            Y_i = (atom_i == Bat)

            Ni = Nvals[i]
            ioff = ao_offset[i]
            I = (ioff+1):(ioff+Ni)

            for j in 1:BS.nshells
                atom_j = BS.shells[j].atom
                X_j = (atom_j == Aat)
                Y_j = (atom_j == Bat)

                (X_i || X_j || Ca) || continue
                (Y_i || Y_j || Cb) || continue

                # Only combinations touching "op" (this nucleus c) at least
                # once are needed here -- pure bra/ket combinations are
                # already fully accounted for in Part 1.
                need_bo = (X_i && Cb) || (Ca && Y_i)
                need_ko = (X_j && Cb) || (Ca && Y_j)
                need_oo = Ca && Cb
                (need_bo || need_ko || need_oo) || continue

                Nj = Nvals[j]
                joff = ao_offset[j]
                J = (joff+1):(joff+Nj)
                Nij = Ni*Nj
                lib = BS.lib

                bb = reshape(view(bb_buf, 1:Ni, 1:Nj, :, :), Ni, Nj, 3, 3)
                bk = reshape(view(bk_buf, 1:Ni, 1:Nj, :, :), Ni, Nj, 3, 3)
                kk = reshape(view(kk_buf, 1:Ni, 1:Nj, :, :), Ni, Nj, 3, 3)
                bkT = reshape(view(bkT_buf, 1:Ni, 1:Nj, :, :), Ni, Nj, 3, 3)
                block = view(block_buf, 1:Ni, 1:Nj, :, :)
                block .= 0.0

                bufv = view(buf, 1:9*Nij)

                # The three permutations below used to go through
                # permutedims!, which allocates ~384 B per call even as the
                # in-place form -- 3 of them per (nucleus, shell pair) was the
                # bulk of this routine's allocation. Folded into scalar fills.
                ∇2nuclear_rinv_μμ!(bufv, BS, env_c, i, j)          # raw (Ni,Nj,3,3)
                @inbounds for q = 1:3, p = 1:3, b = 1:Nj, a = 1:Ni
                    bb[a,b,p,q] = bufv[Nij*(p-1) + 3*Nij*(q-1) + a + Ni*(b-1)]
                end

                # _μν! returns its derivative axes in (ket,bra) order -- see
                # its docstring -- so bk[a,b,p,q] = raw[a,b,q,p].
                ∇2nuclear_rinv_μν!(bufv, BS, env_c, i, j)          # raw (Ni,Nj,3,3)
                @inbounds for q = 1:3, p = 1:3, b = 1:Nj, a = 1:Ni
                    bb_idx = Nij*(q-1) + 3*Nij*(p-1) + a + Ni*(b-1)
                    bk[a,b,p,q] = bufv[bb_idx]
                end

                ∇2nuclear_rinv_μμ!(bufv, BS, env_c, j, i)          # raw (Nj,Ni,3,3)
                @inbounds for q = 1:3, p = 1:3, b = 1:Nj, a = 1:Ni
                    kk[a,b,p,q] = bufv[Nij*(p-1) + 3*Nij*(q-1) + b + Nj*(a-1)]
                end

                # bk is NOT symmetric under its own (3,3)-axis swap (unlike
                # bb/kk, which are same-point mixed partials and always are),
                # so d2/dbra_p dop_q = Zc*(bb+bk)[p,q] and
                # d2/dket_p dop_q = Zc*(bk^T+kk)[p,q] genuinely differ from
                # their (p,q)<->(q,k) counterparts -- each of the two
                # activation sites below needs its own orientation.
                @inbounds for q = 1:3, p = 1:3, b = 1:Nj, a = 1:Ni
                    bkT[a,b,p,q] = bk[a,b,q,p]
                end

                if need_bo
                    (X_i && Cb) && (block .+= Zc .* (bb .+ bk))    # d2/dbra_m dop_k
                    (Ca && Y_i) && (block .+= Zc .* (bb .+ bkT))   # d2/dop_m dbra_k
                end

                if need_ko
                    (X_j && Cb) && (block .+= Zc .* (bkT .+ kk))  # d2/dket_m dop_k
                    (Ca && Y_j) && (block .+= Zc .* (bk .+ kk))   # d2/dop_m dket_k
                end

                if need_oo
                    block .+= (-Zc) .* (bb .+ bk .+ bkT .+ kk)
                end

                out[I,J,:,:] .+= block
            end
        end
    end

    return out
end
