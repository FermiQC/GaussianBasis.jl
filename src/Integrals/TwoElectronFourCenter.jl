# Returned elements are in whatever order tasks happened to finish in, not
# sorted by index -- callers that need a canonical unique integral once (like
# Fock builds) don't care about order, but don't rely on it either.
function sparseERI_2e4c(BS::BasisSet, cutoff = 1e-12)
    # Pre compute a list of angular momentum numbers (l) for each shell
    Nvals = num_basis.(BS.basis)
    Nmax = maximum(Nvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(Nvals[1:(i-1)]) - 1 for i = 1:BS.nshells]

    # Unique shell pairs with i < j
    num_ij = (BS.nshells^2 - BS.nshells) ÷ 2 + BS.nshells

    # Pre allocate array to save ij pairs
    ij_vals = Vector{NTuple{2,Int32}}(undef, num_ij)

    # Pre allocate array to save σij, that is the screening parameter for Schwarz
    σvals = zeros(Cdouble, num_ij)

    ### Loop thorugh i and j such that i ≤ j. Save each pair into ij_vals and compute √σij for integral screening
    tmp = zeros(Cdouble, Nmax^4)
    lim = Int32(BS.nshells - 1)
    @inbounds for i = zero(Int32):lim
        for j = i:lim
            idx = index2(i,j) + 1
            ij_vals[idx] = (i+1,j+1)
            ERI_2e4c!(tmp, BS, i+1, i+1, j+1 ,j+1)
            σvals[idx] = √maximum(tmp)
        end
    end

    ijkl_vals = Iterators.flatten((((ij_vals[ij]..., ij_vals[kl]...) for ij in 1:kl if σvals[ij] * σvals[kl] > cutoff) for kl = 1:num_ij))

    # Upper bound on the number of surviving AO elements, from the block sizes of
    # the shell quartets that passed screening (no ERI calls needed, just Nvals
    # products). Used only to sizehint! the per-task buffers below: geometric
    # Vector growth is amortized O(1), but still copies ~2x the final size in
    # total over the doublings, so hinting close to the true size avoids that
    # overshoot (measured ~35% less peak memory on a dense system with the hint).
    size_ub = 0
    for (i,j,k,l) in ijkl_vals
        size_ub += Nvals[i]*Nvals[j]*Nvals[k]*Nvals[l]
    end

    # Results are not written into a shared nbas-sized dense buffer (whose size
    # would depend only on nbas, not on how many quartets survive screening).
    # Each task instead accumulates into its own private growable buffer, sized
    # upfront via the estimate above, and buffers are concatenated into
    # right-sized output arrays once every task is done. Peak memory then scales
    # with the number of surviving elements instead of nbas^4.
    ntasks = Threads.nthreads()
    chunksize = 10

    # ijkl_vals is a flatten-of-filtered-generators, and eltype() on that gives up
    # and infers Any even though every element is really a NTuple{4,Int32} -- so
    # the channel's element type is pinned explicitly here. Without this, (i,j,k,l)
    # inside the hot loop below become dynamically typed, forcing every (I,J,K,L)
    # tuple to be boxed before it can be pushed into a concretely-typed Vector;
    # that alone was responsible for a many-fold blowup in both time and memory.
    requests = Channel{Vector{NTuple{4,Int32}}}(Inf)
    for chunk in Iterators.partition(ijkl_vals, chunksize)
        put!(requests, chunk)
    end
    close(requests)

    per_task_hint = cld(size_ub, ntasks)
    task_vals = [sizehint!(Cdouble[], per_task_hint) for _ = 1:ntasks]
    task_idxs = [sizehint!(Vector{NTuple{4,Int16}}(), per_task_hint) for _ = 1:ntasks]

    @sync for t = 1:ntasks
        Threads.@spawn begin
            buf = zeros(Cdouble, Nmax^4)
            vals = task_vals[t]
            idxs = task_idxs[t]
            for chunk in requests
                for (i,j,k,l) in chunk
                    @inbounds begin
                        Li, Lj, Lk, Ll = map(n->Nvals[n], (i,j,k,l))
                        Lij = Li*Lj
                        Lijk = Lij*Lk

                        ioff, joff, koff, loff = map(n->ao_offset[n], (i,j,k,l))

                        # Only when the same shell-pair is reused for bra and ket does this
                        # block become self-symmetric (is,js and ks,ls range over the exact
                        # same configuration space), guaranteeing every (I,J,K,L) has a mirror
                        # (K,L,I,J) also enumerated below. Elsewhere IJ vs KL ordering carries
                        # no such redundancy and both must be kept.
                        self_paired = i == k && j == l

                        # Compute ERI
                        ERI_2e4c!(buf, BS, i, j, k ,l)

                        # is, js, ks, ls are indexes within the shell e.g. for a p shell is = (1, 2, 3)
                        # bl, bkl, bjkl are used to map the (i,j,k,l) index into a one-dimensional index for buf
                        for ls = 1:Ll
                            L = loff + ls
                            bl = Lijk*(ls-1)
                            for ks = 1:Lk
                                K = koff + ks
                                L < K && break

                                bkl = Lij*(ks-1) + bl
                                for js = 1:Lj
                                    J = joff + js
                                    bjkl = Li*(js-1) + bkl
                                    for is = 1:Li
                                        I = ioff + is
                                        J < I && break

                                        if self_paired
                                            IJ = (J * (J + 1)) >> 1 + I
                                            KL = (L * (L + 1)) >> 1 + K
                                            IJ > KL && continue
                                        end

                                        v = buf[is + bjkl]
                                        if abs(v) > cutoff
                                            push!(vals, v)
                                            push!(idxs, (I, J, K, L))
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    total = sum(length, task_vals)
    out = Vector{Cdouble}(undef, total)
    indexes = Vector{NTuple{4,Int16}}(undef, total)

    offset = 0
    for t = 1:ntasks
        n = length(task_vals[t])
        n == 0 && continue
        copyto!(out, offset+1, task_vals[t], 1, n)
        copyto!(indexes, offset+1, task_idxs[t], 1, n)
        offset += n
    end

    return indexes, out
end

function ERI_2e4c(BS::BasisSet, i, j, k, l)
    out = zeros(eltype(BS.atoms[1].xyz), num_basis(BS.basis[i]), num_basis(BS.basis[j]),
                num_basis(BS.basis[k]), num_basis(BS.basis[l]))
    ERI_2e4c!(out, BS, i, j, k, l)
    return out
end

function ERI_2e4c!(out, BS::BasisSet{LCint}, i, j, k, l)
    cint2e_sph!(out, @SVector([i,j,k,l]), BS.lib)
end

function ERI_2e4c!(out, BS::BasisSet, i, j, k, l)
    generate_ERI_quartet!(out, BS, i, j, k, l)
end

function ERI_2e4c(BS::BasisSet)
    N = BS.nbas
    out = zeros(N, N, N, N)
    ERI_2e4c!(out, BS)
end

function ERI_2e4c!(out, BS::BasisSet)
    # Save a list containing the number of basis for each shell
    Nvals = num_basis.(BS.basis)
    Nmax = maximum(Nvals)

    # Get slice corresponding to the address in S where the compute chunk goes
    ranges = UnitRange{Int64}[]
    iaccum = 1
    for i = 1:BS.nshells
        push!(ranges, iaccum:(iaccum+ Nvals[i] -1))
        iaccum += Nvals[i]
    end

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

    # Initialize array for results
    allocate(body) = body(zeros(Cdouble, Nmax^4))
    workerpool(allocate, unique_idx; chunksize=10) do (i,j,k,l), buf
        # Shift indexes (C starts with 0, Julia 1)
        id, jd, kd, ld = i+1, j+1, k+1, l+1
        Ni, Nj, Nk, Nl = Nvals[id], Nvals[jd], Nvals[kd], Nvals[ld]

        # Compute ERI
        ERI_2e4c!(buf, BS, id, jd, kd, ld)

        # Move results to output array
        ri, rj, rk, rl = ranges[id], ranges[jd], ranges[kd], ranges[ld]
        out[ri, rj, rk, rl] .= reshape(@view(buf[1:Ni*Nj*Nk*Nl]), (Ni, Nj, Nk, Nl))

        if i != j && k != l && index2(i,j) != index2(k,l)
            @inbounds for ni = ri, nj = rj, nk = rk, nl = rl
                out[nj, ni, nk, nl] = out[ni, nj, nk, nl]
                out[ni, nj, nl, nk] = out[ni, nj, nk, nl]
                out[nj, ni, nl, nk] = out[ni, nj, nk, nl]
                out[nk, nl, ni, nj] = out[ni, nj, nk, nl]
                out[nl, nk, ni, nj] = out[ni, nj, nk, nl]
                out[nk, nl, nj, ni] = out[ni, nj, nk, nl]
                out[nl, nk, nj, ni] = out[ni, nj, nk, nl]
            end
        elseif k != l && index2(i,j) != index2(k,l)
            @inbounds for ni = ri, nj = rj, nk = rk, nl = rl
                out[ni, nj, nl, nk] = out[ni, nj, nk, nl]
                out[nk, nl, ni, nj] = out[ni, nj, nk, nl]
                out[nl, nk, ni, nj] = out[ni, nj, nk, nl]
            end
        elseif i != j && index2(i,j) != index2(k,l)
            @inbounds for ni = ri, nj = rj, nk = rk, nl = rl
                out[nj, ni, nk, nl] = out[ni, nj, nk, nl]
                out[nk, nl, ni, nj] = out[ni, nj, nk, nl]
                out[nk, nl, nj, ni] = out[ni, nj, nk, nl]
            end
        elseif i != j && k != l
            @inbounds for ni = ri, nj = rj, nk = rk, nl = rl
                out[nj, ni, nk, nl] = out[ni, nj, nk, nl]
                out[ni, nj, nl, nk] = out[ni, nj, nk, nl]
                out[nj, ni, nl, nk] = out[ni, nj, nk, nl]
            end
        elseif index2(i,j) != index2(k,l)
            @inbounds for ni = ri, nj = rj, nk = rk, nl = rl
                out[nk, nl, ni, nj] = out[ni, nj, nk, nl]
            end
        end
    end #sync

    return out
end
