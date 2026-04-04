function ∇ERI_2e4c(BS::BasisSet, iA)
    # Pre allocate output
    out = zeros(BS.nbas, BS.nbas, BS.nbas, BS.nbas, 3)
    return ∇ERI_2e4c!(out, BS, iA)
end

function ∇ERI_2e4c!(out, BS::BasisSet, iA)

    if size(out) != (BS.nbas, BS.nbas, BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (N, N, N, N, 3)."))
    end

    A = BS.atoms[iA]

    # Shell indexes for basis in the atom A
    Ashells = Int[]
    notAshells = Int[]
    for i in 1:BS.nshells
        b = BS.basis[i]
        if b.atom == A
            push!(Ashells, i)
        else
            push!(notAshells, i)
        end
    end

    Nvals = num_basis.(BS.basis)
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:BS.nshells]
    Nmax = maximum(Nvals)
    #buf_arrays = [zeros(Cdouble, 3*Nmax^4) for _ = 1:Threads.nthreads()]


    # Find unique (i,j,k,l) combinations given permutational symmetry
    unique_idx = Tuple{Int, Int, Int, Int}[]
    for i = 1:BS.nshells
        for j = i:BS.nshells # i <= j
            for k = 1:BS.nshells
                for l = k:BS.nshells # k <= l
                    if index2(i-1,j-1) < index2(k-1,l-1)
                        continue
                    end
                    push!(unique_idx, (i,j,k,l))
                end
            end
        end
    end

    allocate(body) = body(zeros(Cdouble, 3*Nmax^4))
    workerpool(allocate, unique_idx; chunksize=10) do (i,j,k,l), buf
        x_in_A = map(in(Ashells), (i,j,k,l))

        # If no basis is centered on A, skip
        # If all basis are centered on A, skip
        if !any(x_in_A) || all(x_in_A)
            return
        end

        Ni = Nvals[i]
        Nj = Nvals[j]
        Nk = Nvals[k]
        Nl = Nvals[l]
        Nijkl = Ni*Nj*Nk*Nl

        ioff = ao_offset[i]
        joff = ao_offset[j]
        koff = ao_offset[k]
        loff = ao_offset[l]

        I = (ioff+1):(ioff+Ni)
        J = (joff+1):(joff+Nj)
        K = (koff+1):(koff+Nk)
        L = (loff+1):(loff+Nl)

        # [i'j|kl]
        if x_in_A[1]
            cint2e_ip1_sph!(buf, @SVector([i,j,k,l]), BS.lib)
            ∇q = reshape(view(buf, 1:3*Nijkl), Ni, Nj, Nk, Nl, 3)
            for q in 1:3
                @inbounds for (rl, nl) in enumerate(L), (rk, nk) in enumerate(K), (rj, nj) in enumerate(J), (ri, ni) in enumerate(I)
                    out[ni, nj, nk, nl, q] -= ∇q[ri, rj, rk, rl, q]
                end
            end
        end

        # [ij'|kl]
        if x_in_A[2]
            cint2e_ip1_sph!(buf, @SVector([j,i,k,l]), BS.lib)
            ∇q = reshape(view(buf, 1:3*Nijkl), Nj, Ni, Nk, Nl, 3)
            for q in 1:3
                @inbounds for (rl, nl) in enumerate(L), (rk, nk) in enumerate(K), (rj, nj) in enumerate(J), (ri, ni) in enumerate(I)
                    out[ni, nj, nk, nl, q] -= ∇q[rj, ri, rk, rl, q]
                end
            end
        end

        # [ij|k'l]
        if x_in_A[3]
            cint2e_ip1_sph!(buf, @SVector([k,l,i,j]), BS.lib)
            ∇q = reshape(view(buf, 1:3*Nijkl), Nk, Nl, Ni, Nj, 3)
            for q in 1:3
                @inbounds for (rl, nl) in enumerate(L), (rk, nk) in enumerate(K), (rj, nj) in enumerate(J), (ri, ni) in enumerate(I)
                    out[ni, nj, nk, nl, q] -= ∇q[rk, rl, ri, rj, q]
                end
            end
        end

        # [ij|kl']
        if x_in_A[4]
            cint2e_ip1_sph!(buf, @SVector([l,k,i,j]), BS.lib)
            ∇q = reshape(view(buf, 1:3*Nijkl), Nl, Nk, Ni, Nj, 3)
            for q in 1:3
                @inbounds for (rl, nl) in enumerate(L), (rk, nk) in enumerate(K), (rj, nj) in enumerate(J), (ri, ni) in enumerate(I)
                    out[ni, nj, nk, nl, q] -= ∇q[rl, rk, ri, rj, q]
                end
            end
        end

        for q in 1:3
            if i != j && k != l && index2(i,j) != index2(k,l)
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[nj, ni, nk, nl, q] = out[ni, nj, nk, nl, q]
                    out[ni, nj, nl, nk, q] = out[ni, nj, nk, nl, q]
                    out[nj, ni, nl, nk, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, ni, nj, q] = out[ni, nj, nk, nl, q]
                    out[nl, nk, ni, nj, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, nj, ni, q] = out[ni, nj, nk, nl, q]
                    out[nl, nk, nj, ni, q] = out[ni, nj, nk, nl, q]
                end
            elseif k != l && index2(i,j) != index2(k,l)
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[ni, nj, nl, nk, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, ni, nj, q] = out[ni, nj, nk, nl, q]
                    out[nl, nk, ni, nj, q] = out[ni, nj, nk, nl, q]
                end
            elseif i != j && index2(i,j) != index2(k,l)
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[nj, ni, nk, nl, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, ni, nj, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, nj, ni, q] = out[ni, nj, nk, nl, q]
                end
            elseif i != j && k != l 
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[nj, ni, nk, nl, q] = out[ni, nj, nk, nl, q]
                    out[ni, nj, nl, nk, q] = out[ni, nj, nk, nl, q]
                    out[nj, ni, nl, nk, q] = out[ni, nj, nk, nl, q]
                end
            elseif index2(i,j) != index2(k,l) 
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[nk, nl, ni, nj, q] = out[ni, nj, nk, nl, q]
                end
            end
        end
    end

    return out
end

function ∇sparseERI_2e4c(BS::BasisSet, iA)

    A = BS.atoms[iA]

    # Number of unique integral elements
    N = ((BS.nbas^2 - BS.nbas) ÷ 2) + BS.nbas
    N = ((N^2 - N) ÷ 2) + N

    # Pre allocate output
    ∇x = zeros(N)
    ∇y = zeros(N)
    ∇z = zeros(N)

    indexes = Array{NTuple{4,Int16}}(undef, N)
    mask = repeat([false], N) 

    # Shell indexes for basis in the atom A
    Ashells = Int[]
    notAshells = Int[]
    for i in 1:BS.nshells
        b = BS.basis[i]
        if b.atom == A
            push!(Ashells, i)
        else
            push!(notAshells, i)
        end
    end

    # Pre compute a list of number of basis for each shell (2l +1)
    Nvals = num_basis.(BS.basis)
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:BS.nshells]
    Nmax = maximum(Nvals)

    # Unique shell pairs with i < j
    num_ij = Int((BS.nshells^2 - BS.nshells)/2) + BS.nshells

    # Save ij pairs
    ij_vals = Array{NTuple{2,Int32}}(undef, num_ij)
    lim = Int32(BS.nshells - 1)
    for i = 1:BS.nshells
        for j = i:BS.nshells
            idx = index2(i-1,j-1) + 1
            ij_vals[idx] = (i,j)
        end
    end

    buf = zeros(Cdouble, 3*Nmax^4)
    
    # i,j,k,l => Shell indexes starting at zero
    # I, J, K, L => AO indexes starting at one
    for ij in eachindex(ij_vals)
    #@sync for ij in eachindex(ij_vals)
        #Threads.@spawn begin
        #@inbounds begin
            i,j = ij_vals[ij]
            Ni, Nj = Nvals[i], Nvals[j]
            Nij = Ni*Nj
            ioff = ao_offset[i]
            joff = ao_offset[j]
            for kl in ij:num_ij
                k,l = ij_vals[kl]

                # If all basis are centered at A, or none is, the derivative is zero
                x_in_A = [x in Ashells for x = (i,j,k,l)]
                if !any(x_in_A) || all(x_in_A)
                    continue
                end

                Nk, Nl = Nvals[k], Nvals[l]
                Nijk = Nij*Nk
                Nijkl = Nijk*Nl
                koff = ao_offset[k]
                loff = ao_offset[l]

                # NOTE: Using loops instead of array operations could make this more efficient
                # Compute ERI
                bufx = zeros(Cdouble, Int(Ni), Int(Nj), Int(Nk), Int(Nl))
                bufy = zeros(Cdouble, Int(Ni), Int(Nj), Int(Nk), Int(Nl))
                bufz = zeros(Cdouble, Int(Ni), Int(Nj), Int(Nk), Int(Nl))
                if x_in_A[1]
                    cint2e_ip1_sph!(buf, [i,j,k,l], BS.lib)
                    bufx += reshape(buf[1:Nijkl], Int(Ni), Int(Nj), Int(Nk), Int(Nl))
                    bufy += reshape(buf[Nijkl+1:2*Nijkl], Int(Ni), Int(Nj), Int(Nk), Int(Nl))
                    bufz += reshape(buf[2*Nijkl+1:3*Nijkl], Int(Ni), Int(Nj), Int(Nk), Int(Nl))
                end

                if x_in_A[2]
                    cint2e_ip1_sph!(buf, [j,i,k,l], BS.lib)
                    bufx += permutedims(reshape(buf[1:Nijkl],           Int(Nj), Int(Ni), Int(Nk), Int(Nl)), (2,1,3,4))
                    bufy += permutedims(reshape(buf[Nijkl+1:2*Nijkl],   Int(Nj), Int(Ni), Int(Nk), Int(Nl)), (2,1,3,4))
                    bufz += permutedims(reshape(buf[2*Nijkl+1:3*Nijkl], Int(Nj), Int(Ni), Int(Nk), Int(Nl)), (2,1,3,4))
                end

                if x_in_A[3]
                    cint2e_ip1_sph!(buf, [k,l,i,j], BS.lib)
                    bufx += permutedims(reshape(buf[1:Nijkl],           Int(Nk), Int(Nl), Int(Ni), Int(Nj)), (3,4,1,2))
                    bufy += permutedims(reshape(buf[Nijkl+1:2*Nijkl],   Int(Nk), Int(Nl), Int(Ni), Int(Nj)), (3,4,1,2))
                    bufz += permutedims(reshape(buf[2*Nijkl+1:3*Nijkl], Int(Nk), Int(Nl), Int(Ni), Int(Nj)), (3,4,1,2))
                end

                if x_in_A[4]
                    cint2e_ip1_sph!(buf, [l,k,i,j], BS.lib)
                    bufx += permutedims(reshape(buf[1:Nijkl],           Int(Nl), Int(Nk), Int(Ni), Int(Nj)), (3,4,2,1))
                    bufy += permutedims(reshape(buf[Nijkl+1:2*Nijkl],   Int(Nl), Int(Nk), Int(Ni), Int(Nj)), (3,4,2,1))
                    bufz += permutedims(reshape(buf[2*Nijkl+1:3*Nijkl], Int(Nl), Int(Nk), Int(Ni), Int(Nj)), (3,4,2,1))
                end

                ### This block aims to retrieve unique elements within buf and map them to AO indexes
                # is, js, ks, ls are indexes within the shell e.g. for a p shell is = (1, 2, 3)
                # bl, bkl, bjkl are used to map the (i,j,k,l) index into a one-dimensional index for buf
                # That is, get the correct integrals for the AO quartet.
                for ls = 1:Nl
                    L = loff + ls
                    bl = Nijk*(ls-1)
                    for ks = 1:Nk
                        K = koff + ks
                        L < K ? break : nothing

                        # L ≥ K
                        # index2 for K,L
                        KL = ((L * (L - 1)) ÷ 2) + K - 1                            

                        bkl = Nij*(ks-1) + bl
                        for js = 1:Nj
                            J = joff + js
                            bjkl = Ni*(js-1) + bkl
                            for is = 1:Ni
                                I = ioff + is
                                J < I ? break : nothing

                                # index2 for I,J
                                IJ = ((J * (J - 1)) ÷ 2) + I - 1

                                idx = index2(IJ, KL) + 1
                                ∇x[idx] = -bufx[is + bjkl]
                                ∇y[idx] = -bufy[is + bjkl]
                                ∇z[idx] = -bufz[is + bjkl]
                                if KL ≥ IJ
                                    indexes[idx] = (I, J, K, L)
                                else
                                    indexes[idx] = (K, L, I, J)
                                end
                                mask[idx] = true
                            end
                        end
                    end
                end
            end
        #end #inbounds
        #end #spawn
    end #sync
    return indexes[mask], ∇x[mask], ∇y[mask], ∇z[mask]
end

function ∇ERI_2e3c(BS1::BasisSet, BS2::BasisSet, iA)
    # Pre allocate output
    out = zeros(BS1.nbas, BS1.nbas, BS2.nbas, 3)
    return ∇ERI_2e3c!(out, BS1, BS2, iA)
end

function ∇ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet, iA)

    atoms = unique(vcat(BS1.atoms, BS2.atoms))
    basis = vcat(BS1.basis, BS2.basis)

    Bmerged = BasisSet("$(BS1.name*BS2.name)", atoms, basis)

    if size(out) != (BS1.nbas, BS1.nbas, BS2.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (N1, N1, N2, 3)."))
    end

    A = BS1.atoms[iA]

    # Shell indexes for basis 1 in the atom A
    Ashells1 = Int[]
    notAshells1 = Int[]
    for i in 1:BS1.nshells
        b = BS1.basis[i]
        if b.atom == A
            push!(Ashells1, i)
        else
            push!(notAshells1, i)
        end
    end

    # Shell indexes for basis 2 in the atom A
    Ashells2 = Int[]
    notAshells2 = Int[]
    for i in 1:BS2.nshells
        b = BS2.basis[i]
        if b.atom == A
            push!(Ashells2, i)
        else
            push!(notAshells2, i)
        end
    end

    Nvals1 = num_basis.(BS1.basis)
    ao_offset1 = [sum(Nvals1[1:(i-1)]) for i = 1:BS1.nshells]
    Nmax1 = maximum(Nvals1)

    Nvals2 = num_basis.(BS2.basis)
    ao_offset2 = [sum(Nvals2[1:(i-1)]) for i = 1:BS2.nshells]
    Nmax2 = maximum(Nvals2)

    buf = zeros(Cdouble, 3*Nmax1^2*Nmax2)

    for i = 1:BS1.nshells
        for j = i:BS1.nshells # i <= j
            for k = 1:BS2.nshells

                x_in_A = [i in Ashells1, j in Ashells1, k in Ashells2]

                # If no basis is centered on A, skip
                # If all basis are centered on A, skip
                if !any(x_in_A) || all(x_in_A)
                    continue
                end

                Ni = Nvals1[i]
                Nj = Nvals1[j]
                Nk = Nvals2[k]
                Nijk = Ni*Nj*Nk

                ioff = ao_offset1[i]
                joff = ao_offset1[j]
                koff = ao_offset2[k]

                I = (ioff+1):(ioff+Ni)
                J = (joff+1):(joff+Nj)
                K = (koff+1):(koff+Nk)

	        # [i'j|k]
                if x_in_A[1]
                    cint3c2e_ip1_sph!(buf, [i,j,k+BS1.nshells], Bmerged.lib)
                    for q in 1:3
                        r = (1+Nijk*(q-1)):(q*Nijk)
                        ∇q = reshape(buf[r], Int(Ni), Int(Nj), Int(Nk))
                        out[I,J,K,q] += -∇q
                    end
                end

                # [ij'|k]
                if x_in_A[2]
                    cint3c2e_ip1_sph!(buf, [j,i,k+BS1.nshells], Bmerged.lib)
                    for q in 1:3
                        r = (1+Nijk*(q-1)):(q*Nijk)
                        ∇q = reshape(buf[r], Int(Nj), Int(Ni), Int(Nk))
                        out[I,J,K,q] += -permutedims(∇q, (2,1,3))
                    end
                end

                # [ij|k']
                if x_in_A[3]
                    cint3c2e_ip2_sph!(buf, [i,j,k+BS1.nshells], Bmerged.lib)
                    for q in 1:3
                        r = (1+Nijk*(q-1)):(q*Nijk)
                        ∇q = reshape(buf[r], Int(Ni), Int(Nj), Int(Nk))
                        out[I,J,K,q] += -∇q
                    end
                end

                if i != j
                    for q in 1:3
                        # i,j permutation
                        out[J, I, K, q] .= permutedims(out[I, J, K, q], (2,1,3))
                    end
                end
            end
        end
    end

    return out
end

function ∇ERI_2e2c(BS::BasisSet, iA)
    # Pre allocate output
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇ERI_2e2c!(out, BS, iA)
end

function ∇ERI_2e2c!(out, BS::BasisSet, iA)

    if size(out) != (BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (N, N, 3)."))
    end

    A = BS.atoms[iA]

    # Shell indexes for basis in the atom A
    Ashells = Int[]
    notAshells = Int[]
    for i in 1:BS.nshells
        b = BS.basis[i]
        if b.atom == A
            push!(Ashells, i)
        else
            push!(notAshells, i)
        end
    end

    Nvals = num_basis.(BS.basis)
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:BS.nshells]
    Nmax = maximum(Nvals)

    buf = zeros(Cdouble, 3*Nmax^2)

    for i = 1:BS.nshells
        for j = i:BS.nshells # i <= j

            x_in_A = [x in Ashells for x = (i,j)]

            # If no basis is centered on A, skip
            # If all basis are centered on A, skip
            if !any(x_in_A) || all(x_in_A)
                continue
            end

            Ni = Nvals[i]
            Nj = Nvals[j]
            Nij = Ni*Nj

            ioff = ao_offset[i]
            joff = ao_offset[j]

            I = (ioff+1):(ioff+Ni)
            J = (joff+1):(joff+Nj)

            # [i'|j]
            if x_in_A[1]
                cint2c2e_ip1_sph!(buf, [i,j], BS.lib)
                for q in 1:3
                    r = (1+Nij*(q-1)):(q*Nij)
                    ∇q = reshape(buf[r], Int(Ni), Int(Nj))
                    out[I,J,q] += -∇q
                end
            end

            # [i|j']
            if x_in_A[2]
                cint2c2e_ip1_sph!(buf, [j,i], BS.lib)
                for q in 1:3
                    r = (1+Nij*(q-1)):(q*Nij)
                    ∇q = reshape(buf[r], Int(Nj), Int(Ni))
                    out[I,J,q] += -permutedims(∇q, (2,1))
                end
            end

            if i != j
                for q in 1:3
                    # i,j permutation
                    out[J, I, q] .= permutedims(out[I, J, q], (2,1))
                end
            end
        end
    end

    return out
end
