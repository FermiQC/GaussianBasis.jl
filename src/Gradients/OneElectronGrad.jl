function ∇1e(BS::BasisSet, compute::String, iA)
    # Pre allocate output
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇1e!(out, BS, compute, iA)
end

function ∇1e!(out, BS::BasisSet, compute::String, iA)

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
    buf_arrays = [zeros(Cdouble, 3*Nmax^2) for _ = 1:Threads.nthreads()]
    @sync for i in Ashells
        Threads.@spawn begin
        @inbounds begin
            Ni = Nvals[i]
            ioff = ao_offset[i]
            buf = buf_arrays[Threads.threadid()]
            for j in notAshells
                Nj = Nvals[j]
                joff = ao_offset[j]
                Nij = Ni*Nj
                # Call libcint
                libcint_1e!(buf, [i,j], BS.lib)
                I = (ioff+1):(ioff+Ni)
                J = (joff+1):(joff+Nj)

                # Get strides for each cartesian
                for k in 1:3
                    r = (1+Nij*(k-1)):(k*Nij)
                    @views ∇k = buf[r]
                    #∇k = buf[r]
                    out[I,J,k] .-= reshape(∇k, Int(Ni), Int(Nj))
                end

                # Copy over the transpose
                out[J,I,:] .+= permutedims(out[I,J,:], (2,1,3))
            end
        end #inbounds
        end #spawn 
    end #sync
    return out
end

∇overlap(BS::BasisSet, iA) = ∇1e(BS, "overlap", iA)
∇overlap!(out, BS::BasisSet, iA) = ∇1e!(out, BS, "overlap", iA)
∇kinetic(BS::BasisSet, iA) = ∇1e(BS, "kinetic", iA)
∇kinetic!(out, BS::BasisSet, iA) = ∇1e!(out, BS, "kinetic", iA)

function ∇nuclear(BS::BasisSet, iA)
    # Pre allocate output
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇nuclear!(out, BS, iA)
end

function ∇nuclear!(out, BS::BasisSet, iA)

    if size(out) != (BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3)"))
    end

    A = BS.atoms[iA]

    # Fudge lc_atoms
    lc_atoms_A = deepcopy(BS.lib.atm)
    lc_atoms_woA = deepcopy(BS.lib.atm)
    for i = eachindex(BS.atoms)
        if i == iA
            lc_atoms_woA[1 + 6*(i-1)] = 0
            continue 
        end
        lc_atoms_A[1 + 6*(i-1)] = 0
    end

    # Shell indexes for basis in the atom A (C notation: Starts from 0)
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
    buf_arrays = [zeros(Cdouble, 3*Nmax^2) for _ = 1:Threads.nthreads()]
    # i ∉ A & j ∉ A
    @sync for i in notAshells
        Threads.@spawn begin
        @inbounds begin
            Ni = Nvals[i]
            ioff = ao_offset[i]
            I = (ioff+1):(ioff+Ni)
            buf = buf_arrays[Threads.threadid()]
            for j in notAshells
                Nj = Nvals[j]
                Nij = Ni*Nj
                joff = ao_offset[j]
                J = (joff+1):(joff+Nj)

                # + ⟨i'|Va|j⟩ + ⟨i|Va|j'⟩   (Note that Va is the potential of the nuclei A alone!!)
                cint1e_ipnuc_sph!(buf, Cint.([i-1,j-1]), lc_atoms_A, BS.lib.natm, BS.lib.bas, BS.lib.nbas, BS.lib.env)

                # Get strides for each cartesian
                for k in 1:3
                    r = (1+Nij*(k-1)):(k*Nij)
                    ∇k = reshape(buf[r], Int(Ni), Int(Nj))
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
            Ni = Nvals[i]
            ioff = ao_offset[i]
            I = (ioff+1):(ioff+Ni)
            buf = buf_arrays[Threads.threadid()]
            for j in Ashells
                Nj = Nvals[j]
                Nij = Ni*Nj
                joff = ao_offset[j]
                J = (joff+1):(joff+Nj)

                # - ⟨i'|∑Vc|j⟩ - ⟨i|∑Vc|j'⟩ c != a
                cint1e_ipnuc_sph!(buf, Cint.([i-1,j-1]), lc_atoms_woA, BS.lib.natm, BS.lib.bas, BS.lib.nbas, BS.lib.env)

                for k in 1:3
                    r = (1+Nij*(k-1)):(k*Nij)
                    ∇k = reshape(buf[r], Int(Ni), Int(Nj))
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
            Ni = Nvals[i]
            ioff = ao_offset[i]
            I = (ioff+1):(ioff+Ni)
            buf = buf_arrays[Threads.threadid()]
            for j in notAshells
                Nj = Nvals[j]
                Nij = Ni*Nj
                joff = ao_offset[j]
                J = (joff+1):(joff+Nj)

                # - ⟨i'|∑Vc|j⟩ + ⟨i|Va|j'⟩ c != a
                cint1e_ipnuc_sph!(buf, Cint.([i-1,j-1]), lc_atoms_woA, BS.lib.natm, BS.lib.bas, BS.lib.nbas, BS.lib.env)
                for k in 1:3
                    r = (1+Nij*(k-1)):(k*Nij)
                    ∇k = buf[r]
                    out[I,J,k] .-= reshape(∇k, Int(Ni), Int(Nj))
                end

                cint1e_ipnuc_sph!(buf, Cint.([j-1,i-1]), lc_atoms_A, BS.lib.natm, BS.lib.bas, BS.lib.nbas, BS.lib.env)
                for k in 1:3
                    r = (1+Nij*(k-1)):(k*Nij)
                    ∇k = buf[r]
                    out[I,J,k] .+= transpose(reshape(∇k, Int(Nj), Int(Ni)))
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