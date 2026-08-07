function ∇21e(BS::BasisSet, compute::String, iA, iB)
    out = zeros(BS.nbas, BS.nbas, 3, 3)
    return ∇21e!(out, BS, compute, iA, iB)
end

# Second derivative of overlap/kinetic w.r.t. two nuclear positions (atoms iA,
# iB). Unlike nuclear attraction, these integrals only depend on two
# positions total (the two shell centers), so -- for a given shell pair (i,j)
# with atom(i)=Ai, atom(j)=Aj -- the only combinations of (dR_iA, dR_iB) that
# can be nonzero are built from whichever of shell i's/shell j's centers
# happen to coincide with atom iA and/or atom iB:
#   - both derivatives land on shell i (Ai==iA and Ai==iB, i.e. iA==iB==Ai):
#     same-center kernel (ipip) on shell i, i.e. cint1e_ipip*_sph!([i,j],...)
#   - one derivative on shell i, one on shell j (cross terms): cross-center
#     kernel (ipXip) on [i,j]
#   - both derivatives land on shell j: same-center kernel on shell j,
#     obtained by calling it with shells reversed ([j,i]) and transposing
#     back into (i,j) AO ordering
# When iA==iB and both shells sit on that atom, all four combinations apply
# simultaneously and add together (chain rule for a doubled variable), which
# fixes the size of the diagonal Hessian block.
function ∇21e!(out, BS::BasisSet, compute::String, iA, iB)

    if size(out) != (BS.nbas, BS.nbas, 3, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3, 3)"))
    end
    out .= 0.0

    if compute == "overlap"
        ipip! = cint1e_ipipovlp_sph!
        ipXip! = cint1e_ipovlpip_sph!
    elseif compute == "kinetic"
        ipip! = cint1e_ipipkin_sph!
        ipXip! = cint1e_ipkinip_sph!
    elseif compute == "nuclear_shellonly"
        # Whole-molecule (all nuclei summed, unweighted by *which* atom pair we
        # target) piece of the nuclear-attraction Hessian where both
        # derivatives land on shell centers. This is only PART of
        # ∇2nuclear! -- the remaining piece, where at least one derivative
        # lands on a nuclear *charge* position, is added separately in
        # NuclearHess.jl via isolated per-nucleus rinv terms (to avoid
        # double-counting this shell-only contribution).
        ipip! = cint1e_ipipnuc_sph!
        ipXip! = cint1e_ipnucip_sph!
    else
        throw(ArgumentError("compute must be \"overlap\", \"kinetic\", or \"nuclear_shellonly\" (use ∇2nuclear! for the full nuclear attraction Hessian)"))
    end

    Aat = BS.atoms[iA]
    Bat = BS.atoms[iB]

    Nvals = num_basis.(BS.basis)
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:BS.nshells]
    Nmax = maximum(Nvals)
    buf = zeros(Cdouble, 9*Nmax^2)

    @inbounds for i in 1:BS.nshells
        atom_i = BS.basis[i].atom
        X_i = (atom_i == Aat)
        Y_i = (atom_i == Bat)
        # No early-exit on X_i||Y_i alone here: shell i not touching A or B is
        # still relevant if shell j does (e.g. i on some other atom, j on A --
        # contributes the single-shell "d2_jj" term below).

        Ni = Nvals[i]
        ioff = ao_offset[i]
        I = (ioff+1):(ioff+Ni)

        for j in 1:BS.nshells
            atom_j = BS.basis[j].atom
            X_j = (atom_j == Aat)
            Y_j = (atom_j == Bat)

            (X_i || X_j) || continue
            (Y_i || Y_j) || continue

            Nj = Nvals[j]
            joff = ao_offset[j]
            J = (joff+1):(joff+Nj)
            Nij = Ni*Nj

            block = zeros(Cdouble, Ni, Nj, 3, 3)

            if X_i && Y_i
                ipip!(buf, [i,j], BS.lib)
                block .+= reshape(view(buf, 1:9*Nij), Ni, Nj, 3, 3)
            end

            if (X_i && Y_j) || (X_j && Y_i)
                # cint1e_ip*ip_sph!'s two derivative-component axes come out
                # in (ket,bra) order, not (bra,ket): raw[:,:,q,p] is
                # d2X/d(shell i)_p d(shell j)_q. For overlap/kinetic this is
                # invisible (those integrals depend only on the shell-shell
                # separation, so their cross-Hessian is provably symmetric
                # under the p<->q swap, making raw and its transpose
                # numerically identical) -- but not in general (e.g. nuclear
                # attraction's shell-only piece, which also has a fixed
                # background charge distribution breaking that symmetry), so
                # the two placements need genuinely different orientations:
                # X_i&&Y_j wants d2/d(shell i=A)_m d(shell j=B)_k = raw[:,:,k,m]
                # X_j&&Y_i wants d2/d(shell j=A)_m d(shell i=B)_k = raw[:,:,m,k]
                ipXip!(buf, [i,j], BS.lib)
                raw = reshape(view(buf, 1:9*Nij), Ni, Nj, 3, 3)
                X_i && Y_j && (block .+= permutedims(raw, (1,2,4,3)))
                X_j && Y_i && (block .+= raw)
            end

            if X_j && Y_j
                ipip!(buf, [j,i], BS.lib)
                blockT = reshape(view(buf, 1:9*Nij), Nj, Ni, 3, 3)
                block .+= permutedims(blockT, (2,1,3,4))
            end

            out[I,J,:,:] .+= block
        end
    end

    return out
end

∇2overlap(BS::BasisSet, iA, iB) = ∇21e(BS, "overlap", iA, iB)
∇2overlap!(out, BS::BasisSet, iA, iB) = ∇21e!(out, BS, "overlap", iA, iB)
∇2kinetic(BS::BasisSet, iA, iB) = ∇21e(BS, "kinetic", iA, iB)
∇2kinetic!(out, BS::BasisSet, iA, iB) = ∇21e!(out, BS, "kinetic", iA, iB)
