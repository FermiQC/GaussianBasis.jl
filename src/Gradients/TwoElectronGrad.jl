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

    buf = zeros(Cdouble, 3*Nmax^4)

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

    for (i,j,k,l) in unique_idx

        x_in_A = [x in Ashells for x = (i,j,k,l)]

        # If no basis is centered on A, skip
        # If all basis are centered on A, skip
        if !any(x_in_A) || all(x_in_A)
            continue
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
            cint2e_ip1_sph!(buf, [i,j,k,l], BS.lib)
            for q in 1:3
                r = (1+Nijkl*(q-1)):(q*Nijkl)
                ∇q = reshape(buf[r], Int(Ni), Int(Nj), Int(Nk), Int(Nl))
                out[I,J,K,L,q] += -∇q
            end
        end

        # [ij'|kl]
        if x_in_A[2]
            cint2e_ip1_sph!(buf, [j,i,k,l], BS.lib)
            for q in 1:3
                r = (1+Nijkl*(q-1)):(q*Nijkl)
                ∇q = reshape(buf[r], Int(Nj), Int(Ni), Int(Nk), Int(Nl))
                out[I,J,K,L,q] += -permutedims(∇q, (2,1,3,4))
            end
        end

        # [ij|k'l]
        if x_in_A[3]
            cint2e_ip1_sph!(buf, [k,l,i,j], BS.lib)
            for q in 1:3
                r = (1+Nijkl*(q-1)):(q*Nijkl)
                ∇q = reshape(buf[r], Int(Nk), Int(Nl), Int(Ni), Int(Nj))
                out[I,J,K,L,q] += -permutedims(∇q, (3,4,1,2))
            end
        end

        # [ij|kl']
        if x_in_A[4]
            cint2e_ip1_sph!(buf, [l,k,i,j], BS.lib)
            for q in 1:3
                r = (1+Nijkl*(q-1)):(q*Nijkl)
                ∇q = reshape(buf[r], Int(Nl), Int(Nk), Int(Ni), Int(Nj))
                out[I,J,K,L,q] += -permutedims(∇q, (3,4,2,1))
            end
        end

        for q in 1:3
            if i != j && k != l && index2(i,j) != index2(k,l)
                # i,j permutation
                out[J, I, K, L, q] .= permutedims(out[I, J, K, L, q], (2,1,3,4))
                # k,l permutation
                out[I, J, L, K, q] .= permutedims(out[I, J, K, L, q], (1,2,4,3))

                # i,j + k,l permutatiom
                out[J, I, L, K, q] .= permutedims(out[I, J, K, L, q], (2,1,4,3))

                # ij, kl permutation
                out[K, L, I, J, q] .= permutedims(out[I, J, K, L, q], (3,4,1,2))
                # ij, kl + k,l permutation
                out[L, K, I, J, q] .= permutedims(out[I, J, K, L, q], (4,3,1,2))
                # ij, kl + i,j permutation
                out[K, L, J, I, q] .= permutedims(out[I, J, K, L, q], (3,4,2,1))
                # ij, kl + i,j + k,l permutation
                out[L, K, J, I, q] .= permutedims(out[I, J, K, L, q], (4,3,2,1))

            elseif k != l && index2(i,j) != index2(k,l)
                # k,l permutation
                out[I, J, L, K, q] .= permutedims(out[I, J, K, L, q], (1,2,4,3))
                # ij, kl permutation
                out[K, L, I, J, q] .= permutedims(out[I, J, K, L, q], (3,4,1,2))
                # ij, kl + k,l permutation
                out[L, K, I, J, q] .= permutedims(out[I, J, K, L, q], (4,3,1,2))

            elseif i != j && index2(i,j) != index2(k,l)
                # i,j permutation
                out[J, I, K, L, q] .= permutedims(out[I, J, K, L, q], (2,1,3,4))

                # ij, kl permutation
                out[K, L, I, J, q] .= permutedims(out[I, J, K, L, q], (3,4,1,2))
                # ij, kl + i,j permutation
                out[K, L, J, I, q] .= permutedims(out[I, J, K, L, q], (3,4,2,1))
        
            elseif i != j && k != l 
                # i,j permutation
                out[J, I, K, L, q] .= permutedims(out[I, J, K, L, q], (2,1,3,4))
                # k,l permutation
                out[I, J, L, K, q] .= permutedims(out[I, J, K, L, q], (1,2,4,3))

                # i,j + k,l permutatiom
                out[J, I, L, K, q] .= permutedims(out[I, J, K, L, q], (2,1,4,3))
            elseif index2(i,j) != index2(k,l) 
                # ij, kl permutation
                out[K, L, I, J, q] .= permutedims(out[I, J, K, L, q], (3,4,1,2))
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

    buf_arrays = [zeros(Cdouble, 3*Nmax^4) for _ = 1:Threads.nthreads()]
    
    # i,j,k,l => Shell indexes starting at zero
    # I, J, K, L => AO indexes starting at one
    for ij in eachindex(ij_vals)
    #@sync for ij in eachindex(ij_vals)
        #Threads.@spawn begin
        #@inbounds begin
            buf = buf_arrays[Threads.threadid()]
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