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

function ∇2nuclear(BS::BasisSet, iA, iB)
    out = zeros(BS.nbas, BS.nbas, 3, 3)
    return ∇2nuclear!(out, BS, iA, iB)
end

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

    Nvals = num_basis.(BS.basis)
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:BS.nshells]
    Nmax = maximum(Nvals)
    buf = zeros(Cdouble, 9*Nmax^2)

    nuclei = Aat == Bat ? (Aat,) : (Aat, Bat)

    @inbounds for c in nuclei
        Zc = c.Z
        Rc = collect(c.xyz) ./ Molecules.bohr_to_angstrom
        env_c = rinv_env(BS, Rc)
        Ca = (c == Aat)
        Cb = (c == Bat)

        for i in 1:BS.nshells
            atom_i = BS.basis[i].atom
            X_i = (atom_i == Aat)
            Y_i = (atom_i == Bat)

            Ni = Nvals[i]
            ioff = ao_offset[i]
            I = (ioff+1):(ioff+Ni)

            for j in 1:BS.nshells
                atom_j = BS.basis[j].atom
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

                cint1e_ipiprinv_sph!(buf, Cint.([i,j].-1), BS.lib.atm, BS.lib.natm, BS.lib.bas, BS.lib.nbas, env_c)
                bb = reshape(copy(view(buf, 1:9*Nij)), Ni, Nj, 3, 3)

                # cint1e_iprinvip_sph!'s two derivative-component axes come out
                # in (ket,bra) order, not (bra,ket) -- unlike the same-shell
                # ipXip kernels (ipovlpip/ipkinip/ipnucip and ipiprinv here),
                # a mixed bra/ket derivative of a *3-center* integral (bra,
                # ket, and the rinv origin are three independent points) has
                # no reason to be symmetric under swapping those two axes, so
                # the raw layout must be corrected explicitly (verified
                # against a direct finite difference of cint1e_iprinv_sph!).
                cint1e_iprinvip_sph!(buf, Cint.([i,j].-1), BS.lib.atm, BS.lib.natm, BS.lib.bas, BS.lib.nbas, env_c)
                bk = permutedims(reshape(copy(view(buf, 1:9*Nij)), Ni, Nj, 3, 3), (1,2,4,3))

                cint1e_ipiprinv_sph!(buf, Cint.([j,i].-1), BS.lib.atm, BS.lib.natm, BS.lib.bas, BS.lib.nbas, env_c)
                kk = permutedims(reshape(copy(view(buf, 1:9*Nij)), Nj, Ni, 3, 3), (2,1,3,4))

                block = zeros(Cdouble, Ni, Nj, 3, 3)

                # bk is NOT symmetric under its own (3,3)-axis swap (unlike
                # bb/kk, which are same-point mixed partials and always are),
                # so d2/dbra_p dop_q = Zc*(bb+bk)[p,q] and
                # d2/dket_p dop_q = Zc*(bk^T+kk)[p,q] genuinely differ from
                # their (p,q)<->(q,k) counterparts -- each of the two
                # activation sites below needs its own orientation.
                bkT = permutedims(bk, (1,2,4,3))

                if need_bo
                    bo = Zc .* (bb .+ bk)      # d2/dbra_m dop_k
                    bo_T = Zc .* (bb .+ bkT)   # d2/dop_m dbra_k = (d2/dbra_k dop_m)
                    (X_i && Cb) && (block .+= bo)
                    (Ca && Y_i) && (block .+= bo_T)
                end

                if need_ko
                    ko = Zc .* (bk .+ kk)      # d2/dop_m dket_k = (d2/dket_k dop_m)
                    ko_T = Zc .* (bkT .+ kk)   # d2/dket_m dop_k
                    (X_j && Cb) && (block .+= ko_T)
                    (Ca && Y_j) && (block .+= ko)
                end

                if need_oo
                    oo = -Zc .* (bb .+ bk .+ bkT .+ kk)
                    block .+= oo
                end

                out[I,J,:,:] .+= block
            end
        end
    end

    return out
end
