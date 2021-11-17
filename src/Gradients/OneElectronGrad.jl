function ∇1e(BS::BasisSet, compute::String, iA, T::DataType = Float64)
    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas, 3)
    return ∇1e!(out, BS, compute, iA, T)
end

function ∇1e!(out, BS::BasisSet, compute::String, iA, T::DataType = Float64)

    if size(out) != (BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3)"))
    end

    if compute == "overlap"
        libcint_1e! =  cint1e_ipovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_ipkin_sph!
    elseif compute == "nuclear"
        return ∇nuclear(BS, iA)
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

    # Shell indexes for basis in the atom A (C notation: Starts from 0)
    Ashells = [(s0 + i -1) for i = 1:length(BS.basis[iA])]
    notAshells = Int[]
    for i = 1:BS.nshells
        if !((i-1) in Ashells)
            push!(notAshells, i-1)
        end
    end

    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]
    Lmax = maximum(lvals)
    buf_arrays = [zeros(Cdouble, 3*Lmax^2) for _ = 1:Threads.nthreads()]
    @sync for i in Ashells
        Threads.@spawn begin
        @inbounds begin
            Li = lvals[i+1]
            ioff = ao_offset[i+1]
            buf = buf_arrays[Threads.threadid()]
            for j in notAshells
                Lj = lvals[j+1]
                joff = ao_offset[j+1]
                Lij = Li*Lj
                # Call libcint
                libcint_1e!(buf, Cint.([i,j]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                I = (ioff+1):(ioff+Li)
                J = (joff+1):(joff+Lj)

                # Get strides for each cartesian
                for k in 1:3
                    r = (1+Lij*(k-1)):(k*Lij)
                    @views ∇k = buf[r]
                    out[I,J,k] .-= reshape(∇k, Int(Li), Int(Lj))
                end

                # Copy over the transpose
                out[J,I,:] .+= permutedims(out[I,J,:], (2,1,3))
            end
        end #inbounds
        end #spawn 
    end #sync
    return out
end

function ∇overlap(BS::BasisSet, iA, T::DataType = Float64)
    return ∇1e(BS, "overlap", iA, T)
end

function ∇overlap!(out, BS::BasisSet, iA, T::DataType = Float64)
    return ∇1e!(out, BS, "overlap", iA, T)
end

function ∇kinetic(BS::BasisSet, iA, T::DataType = Float64)
    return ∇1e(BS, "kinetic", iA, T)
end

function ∇kinetic!(out, BS::BasisSet, iA, T::DataType = Float64)
    return ∇1e!(out, BS, "kinetic", iA, T)
end

function ∇nuclear(BS::BasisSet, iA, T::DataType = Float64)
    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas, 3)
    return ∇nuclear!(out, BS, iA, T)
end

function ∇nuclear!(out, BS::BasisSet, iA, T::DataType = Float64)

    if size(out) != (BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3)"))
    end

    A = BS.atoms[iA]

    # Fudge lc_atoms
    lc_atoms_A = deepcopy(BS.lc_atoms)
    lc_atoms_woA = deepcopy(BS.lc_atoms)
    for i = eachindex(BS.atoms)
        if i == iA
            lc_atoms_woA[1 + 6*(i-1)] = 0
            continue 
        end
        lc_atoms_A[1 + 6*(i-1)] = 0
    end

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
        # i-1 is used to convert to C notation (starting with 0)
        if !((i-1) in Ashells)
            push!(notAshells, i-1)
        end
    end
    
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]
    Lmax = maximum(lvals)
    buf_arrays = [zeros(Cdouble, 3*Lmax^2) for _ = 1:Threads.nthreads()]

    # i ∉ A & j ∉ A
    @sync for i in notAshells
        Threads.@spawn begin
        @inbounds begin
            Li = lvals[i+1]
            ioff = ao_offset[i+1]
            I = (ioff+1):(ioff+Li)
            buf = buf_arrays[Threads.threadid()]
            for j in notAshells
                Lj = lvals[j+1]
                Lij = Li*Lj
                joff = ao_offset[j+1]
                J = (joff+1):(joff+Lj)

                # + ⟨i'|Va|j⟩ + ⟨i|Va|j'⟩   (Note that Va is the potential of the nuclei A alone!!)
                cint1e_ipnuc_sph!(buf, Cint.([i,j]), lc_atoms_A, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                # Get strides for each cartesian
                for k in 1:3
                    r = (1+Lij*(k-1)):(k*Lij)
                    ∇k = reshape(buf[r], Int(Li), Int(Lj))
                    out[I,J,k] .+= ∇k  # ⟨i'|Va|j ⟩
                end
            end
        end # inbounds
        end # #spawn 
    end # sync

    # i ∈ A & j ∈ A
    @sync for i in Ashells
        Threads.@spawn begin
        @inbounds begin
            Li = lvals[i+1]
            ioff = ao_offset[i+1]
            I = (ioff+1):(ioff+Li)
            buf = buf_arrays[Threads.threadid()]
            for j in Ashells
                Lj = lvals[j+1]
                Lij = Li*Lj
                joff = ao_offset[j+1]
                J = (joff+1):(joff+Lj)

                # - ⟨i'|∑Vc|j⟩ - ⟨i|∑Vc|j'⟩ c != a
                cint1e_ipnuc_sph!(buf, Cint.([i,j]), lc_atoms_woA, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                for k in 1:3
                    r = (1+Lij*(k-1)):(k*Lij)
                    ∇k = reshape(buf[r], Int(Li), Int(Lj))
                    out[I,J,k] .-= ∇k  # ⟨i'|∑Vc|j ⟩ c != a
                end
            end
        end #inbounds
        end #spawn
    end #sync

    # i ∈ A & j ∉ A
    @sync for i in Ashells
        Threads.@spawn begin
        @inbounds begin
            Li = lvals[i+1]
            ioff = ao_offset[i+1]
            I = (ioff+1):(ioff+Li)
            buf = buf_arrays[Threads.threadid()]
            for j in notAshells
                Lj = lvals[j+1]
                Lij = Li*Lj
                joff = ao_offset[j+1]
                J = (joff+1):(joff+Lj)

                # - ⟨i'|∑Vc|j⟩ + ⟨i|Va|j'⟩ c != a
                cint1e_ipnuc_sph!(buf, Cint.([i,j]), lc_atoms_woA, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                for k in 1:3
                    r = (1+Lij*(k-1)):(k*Lij)
                    ∇k = buf[r]
                    out[I,J,k] .-= reshape(∇k, Int(Li), Int(Lj))
                end

                cint1e_ipnuc_sph!(buf, Cint.([j,i]), lc_atoms_A, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                for k in 1:3
                    r = (1+Lij*(k-1)):(k*Lij)
                    ∇k = buf[r]
                    out[I,J,k] .+= transpose(reshape(∇k, Int(Lj), Int(Li)))
                end
            end
        end #inbounds
        end #spawn
    end #sync

    # Add transpose values
    # This must be done outside the threaded loops
    # to avoid race conditions. 
    for k in 1:3
        out[:,:,k] .+= out[:,:,k]'
    end

    return out
end