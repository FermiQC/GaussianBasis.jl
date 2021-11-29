function ∇ERI_2e4c(BS::BasisSet, iA, T::DataType = Float64)
    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas, BS.nbas, BS.nbas, 3)
    return ∇ERI_2e4c!(out, BS, iA, T)
end

function ∇ERI_2e4c!(out, BS::BasisSet, iA, T::DataType = Float64)

    if size(out) != (BS.nbas, BS.nbas, BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (N, N, N, N, 3)."))
    end

    A = BS.atoms[iA]

    # Get shell index (s0) where basis of the desired atom start
    s0 = 0
    for a in eachindex(BS.atoms)
        if BS.atoms[a] != A
            s0 += length(BS.basis[a])
        else
            break
        end
    end

    # Shell indexes for basis in the atom A
    Ashells = [(s0 + i -1) for i = 1:length(BS.basis[iA])]
    notAshells = Int[]
    for i = 1:BS.nshells
        if !((i-1) in Ashells)
            push!(notAshells, i-1)
        end
    end

    bas_per_shell = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    ao_offset = [sum(bas_per_shell[1:(i-1)]) for i = 1:BS.nshells]
    NBmax = maximum(bas_per_shell)
    buf = zeros(Cdouble, 3*NBmax^4)

    # Find unique (i,j,k,l) combinations given permutational symmetry
    unique_idx = NTuple{4,Int16}[]
    N = Int16(BS.nshells - 1)
    ZERO = zero(Int16)
    for i = ZERO:N
        for j = i:N
            for k = ZERO:N
                for l = k:N
                    if index2(i,j) < index2(k,l)
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

        Li = bas_per_shell[i+1]
        Lj = bas_per_shell[j+1]
        Lk = bas_per_shell[k+1]
        Ll = bas_per_shell[l+1]
        Lijkl = Li*Lj*Lk*Ll

        ioff = ao_offset[i+1]
        joff = ao_offset[j+1]
        koff = ao_offset[k+1]
        loff = ao_offset[l+1]

        I = (ioff+1):(ioff+Li)
        J = (joff+1):(joff+Lj)
        K = (koff+1):(koff+Lk)
        L = (loff+1):(loff+Ll)

        # [i'j|kl]
        if x_in_A[1]
            cint2e_ip1_sph!(buf, Cint.([i,j,k,l]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            for q in 1:3
                r = (1+Lijkl*(q-1)):(q*Lijkl)
                ∇q = reshape(buf[r], Int(Li), Int(Lj), Int(Lk), Int(Ll))
                out[I,J,K,L,q] += -∇q
            end
        end

        # [ij'|kl]
        if x_in_A[2]
            cint2e_ip1_sph!(buf, Cint.([j,i,k,l]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            for q in 1:3
                r = (1+Lijkl*(q-1)):(q*Lijkl)
                ∇q = reshape(buf[r], Int(Lj), Int(Li), Int(Lk), Int(Ll))
                out[I,J,K,L,q] += -permutedims(∇q, (2,1,3,4))
            end
        end

        # [ij|k'l]
        if x_in_A[3]
            cint2e_ip1_sph!(buf, Cint.([k,l,i,j]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            for q in 1:3
                r = (1+Lijkl*(q-1)):(q*Lijkl)
                ∇q = reshape(buf[r], Int(Lk), Int(Ll), Int(Li), Int(Lj))
                out[I,J,K,L,q] += -permutedims(∇q, (3,4,1,2))
            end
        end

        # [ij|kl']
        if x_in_A[4]
            cint2e_ip1_sph!(buf, Cint.([l,k,i,j]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            for q in 1:3
                r = (1+Lijkl*(q-1)):(q*Lijkl)
                ∇q = reshape(buf[r], Int(Ll), Int(Lk), Int(Li), Int(Lj))
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

function ∇sparseERI_2e4c(BS::BasisSet, iA, T::DataType = Float64)

    A = BS.atoms[iA]

    # Number of unique integral elements
    N = Int((BS.nbas^2 - BS.nbas)/2) + BS.nbas
    N = Int((N^2 - N)/2) + N

    # Pre allocate output
    ∇x = zeros(T, N)
    ∇y = zeros(T, N)
    ∇z = zeros(T, N)
    indexes = Array{NTuple{4,Int16}}(undef, N)
    mask = repeat([false], N) 

    # Get shell index (s0) where basis of the desired atom start
    s0 = 0
    for a in eachindex(BS.atoms)
        if BS.atoms[a] != A
            s0 += length(BS.basis[a])
        else
            break
        end
    end

    # Shell indexes for basis in the atom A
    # C indexing
    Ashells = [(s0 + i -1) for i = 1:length(BS.basis[iA])]
    notAshells = Int[]
    for i = 1:BS.nshells
        if !((i-1) in Ashells)
            push!(notAshells, i-1)
        end
    end

    # Pre compute a list of number of basis for each shell (2l +1)
    bas_per_shell = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    NBmax = maximum(bas_per_shell)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(bas_per_shell[1:(i-1)]) - 1 for i = 1:BS.nshells]

    # Unique shell pairs with i < j
    num_ij = Int((BS.nshells^2 - BS.nshells)/2) + BS.nshells

    # Save ij pairs
    ij_vals = Array{NTuple{2,Int32}}(undef, num_ij)
    lim = Int32(BS.nshells - 1)
    for i = UnitRange{Int32}(zero(Int32),lim)
        for j = UnitRange{Int32}(i, lim)
            idx = index2(i,j) + 1
            ij_vals[idx] = (i,j)
        end
    end

    buf_arrays = [zeros(Cdouble, 3*NBmax^4) for _ = 1:Threads.nthreads()]
    
    # i,j,k,l => Shell indexes starting at zero
    # I, J, K, L => AO indexes starting at one
    for ij in eachindex(ij_vals)
    #@sync for ij in eachindex(ij_vals)
        #Threads.@spawn begin
        #@inbounds begin
            buf = buf_arrays[Threads.threadid()]
            i,j = ij_vals[ij]
            Li, Lj = bas_per_shell[i+1], bas_per_shell[j+1]
            Lij = Li*Lj
            ioff = ao_offset[i+1]
            joff = ao_offset[j+1]
            for kl in ij:num_ij
                k,l = ij_vals[kl]

                # If all basis are centered at A, or none is, the derivative is zero
                x_in_A = [x in Ashells for x = (i,j,k,l)]
                if !any(x_in_A) || all(x_in_A)
                    continue
                end

                Lk, Ll = bas_per_shell[k+1], bas_per_shell[l+1]
                Lijk = Lij*Lk
                Lijkl = Lijk*Ll
                koff = ao_offset[k+1]
                loff = ao_offset[l+1]

                # NOTE: Using loops instead of array operations could make this more efficient
                # Compute ERI
                bufx = zeros(Cdouble, Int(Li), Int(Lj), Int(Lk), Int(Ll))
                bufy = zeros(Cdouble, Int(Li), Int(Lj), Int(Lk), Int(Ll))
                bufz = zeros(Cdouble, Int(Li), Int(Lj), Int(Lk), Int(Ll))
                if x_in_A[1]
                    cint2e_ip1_sph!(buf, [i,j,k,l], BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                    bufx += reshape(buf[1:Lijkl], Int(Li), Int(Lj), Int(Lk), Int(Ll))
                    bufy += reshape(buf[Lijkl+1:2*Lijkl], Int(Li), Int(Lj), Int(Lk), Int(Ll))
                    bufz += reshape(buf[2*Lijkl+1:3*Lijkl], Int(Li), Int(Lj), Int(Lk), Int(Ll))
                end

                if x_in_A[2]
                    cint2e_ip1_sph!(buf, [j,i,k,l], BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                    bufx += permutedims(reshape(buf[1:Lijkl],           Int(Lj), Int(Li), Int(Lk), Int(Ll)), (2,1,3,4))
                    bufy += permutedims(reshape(buf[Lijkl+1:2*Lijkl],   Int(Lj), Int(Li), Int(Lk), Int(Ll)), (2,1,3,4))
                    bufz += permutedims(reshape(buf[2*Lijkl+1:3*Lijkl], Int(Lj), Int(Li), Int(Lk), Int(Ll)), (2,1,3,4))
                end

                if x_in_A[3]
                    cint2e_ip1_sph!(buf, [k,l,i,j], BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                    bufx += permutedims(reshape(buf[1:Lijkl],           Int(Lk), Int(Ll), Int(Li), Int(Lj)), (3,4,1,2))
                    bufy += permutedims(reshape(buf[Lijkl+1:2*Lijkl],   Int(Lk), Int(Ll), Int(Li), Int(Lj)), (3,4,1,2))
                    bufz += permutedims(reshape(buf[2*Lijkl+1:3*Lijkl], Int(Lk), Int(Ll), Int(Li), Int(Lj)), (3,4,1,2))
                end

                if x_in_A[4]
                    cint2e_ip1_sph!(buf, [l,k,i,j], BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                    bufx += permutedims(reshape(buf[1:Lijkl],           Int(Ll), Int(Lk), Int(Li), Int(Lj)), (3,4,2,1))
                    bufy += permutedims(reshape(buf[Lijkl+1:2*Lijkl],   Int(Ll), Int(Lk), Int(Li), Int(Lj)), (3,4,2,1))
                    bufz += permutedims(reshape(buf[2*Lijkl+1:3*Lijkl], Int(Ll), Int(Lk), Int(Li), Int(Lj)), (3,4,2,1))
                end

                ### This block aims to retrieve unique elements within buf and map them to AO indexes
                # is, js, ks, ls are indexes within the shell e.g. for a p shell is = (1, 2, 3)
                # bl, bkl, bjkl are used to map the (i,j,k,l) index into a one-dimensional index for buf
                # That is, get the correct integrals for the AO quartet.
                for ls = 1:Ll
                    L = loff + ls
                    bl = Lijk*(ls-1)
                    for ks = 1:Lk
                        K = koff + ks
                        L < K ? break : nothing

                        # L ≥ K
                        KL = (L * (L + 1)) >> 1 + K                            
                        bkl = Lij*(ks-1) + bl
                        for js = 1:Lj
                            J = joff + js
                            bjkl = Li*(js-1) + bkl
                            for is = 1:Li
                                I = ioff + is
                                J < I ? break : nothing

                                IJ = (J * (J + 1)) >> 1 + I

                                idx = index2(IJ,KL) + 1
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