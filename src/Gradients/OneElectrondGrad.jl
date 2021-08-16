function grad_ao_1e(BS::BasisSet, iA, compute::String, T::DataType = Float64)

    if compute == "overlap"
        libcint_1e! =  cint1e_ipovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_ipkin_sph!
    # NUCLEAR DOES NOT WORK
    elseif compute == "nuclear"
        libcint_1e! =  cint1e_ipnuc_sph!
    end

    A = BS.molecule.atoms[iA]

    # Pre allocate output
    out = zeros(T, 3, BS.nbas, BS.nbas)

    # Get shell index (s0) where basis of the desired atom start
    s0 = 0
    for a in BS.molecule.atoms
        if a != A
            s0 += length(BS.basis[a])
        else
            break
        end
    end

    # Shell indexes for basis in the atom A
    Ashells = [(s0 + i -1) for i = 1:length(BS.basis[A])]
    notAshells = Int[]
    for i = 1:BS.nshells
        if !((i-1) in Ashells)
            push!(notAshells, i-1)
        end
    end

    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]
    Lmax = maximum(lvals)
    buf = zeros(Cdouble, 3*Lmax^2)
    for i in Ashells
        Li = lvals[i+1]
        ioff = ao_offset[i+1]
        for j in notAshells
            Lj = lvals[j+1]
            joff = ao_offset[j+1]
            # Call libcint
            libcint_1e!(buf, Cint.([j,i]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            I = (ioff+1):(ioff+Li)
            J = (joff+1):(joff+Lj)
            out[:,I,J] .= reshape(buf[1:(3*Li*Lj)], Int.([3,Li,Lj])...)
            out[:,J,I] .= permutedims(out[:,I,J], (1,3,2))
        end
    end
    return out

    # Pre compute a list of angular momentum numbers (l) for each shell
    #lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    #Lmax = maximum(lvals)

    ## Offset list for each shell, used to map shell index to AO index
    #ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    #buf_arrays = [zeros(Cdouble, Lmax^4) for _ = 1:Threads.nthreads()]

    #@sync for i in 1:BS.nshells
    #    Threads.@spawn begin
    #        @inbounds begin
    #            Li = lvals[i]
    #            buf = buf_arrays[Threads.threadid()]
    #            ioff = ao_offset[i]
    #            for j in i:BS.nshells
    #                Lj = lvals[j]
    #                joff = ao_offset[j]

    #                # Call libcint
    #                libcint_1e!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

    #                # Loop through shell block and save unique elements
    #                for js = 1:Lj
    #                    J = joff + js
    #                    for is = 1:Li
    #                        I = ioff + is
    #                        J < I ? break : nothing
    #                        out[I,J] = buf[is + Li*(js-1)]
    #                        out[J,I] = out[I,J]
    #                    end
    #                end
    #            end
    #        end #inbounds
    #    end #spawn
    #end #sync
    #return out
end

function grad_ao_nuc(BS::BasisSet, iA, T::DataType = Float64)

    A = BS.molecule.atoms[iA]

    # Pre allocate output
    out = zeros(T, 3, BS.nbas, BS.nbas)

    # Get shell index (s0) where basis of the desired atom start
    s0 = 0
    for a in BS.molecule.atoms
        if a != A
            s0 += length(BS.basis[a])
        else
            break
        end
    end

    # Shell indexes for basis in the atom A
    Ashells = [(s0 + i -1) for i = 1:length(BS.basis[A])]
    notAshells = Int[]
    for i = 1:BS.nshells
        if !((i-1) in Ashells)
            push!(notAshells, i-1)
        end
    end

    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]
    Lmax = maximum(lvals)
    buf = zeros(Cdouble, 3*Lmax^2)
    for i in Ashells
        Li = lvals[i+1]
        ioff = ao_offset[i+1]
        for j in notAshells
            Lj = lvals[j+1]
            joff = ao_offset[j+1]
            # Call libcint
            cint1e_ipnuc_sph!(buf, Cint.([j,i]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            I = (ioff+1):(ioff+Li)
            J = (joff+1):(joff+Lj)
            out[:,I,J] .= reshape(buf[1:(3*Li*Lj)], 3,Int(Li),Int(Lj))
            out[:,J,I] .= permutedims(out[:,I,J], (1,3,2))
        end
    end
    return out

    # Pre compute a list of angular momentum numbers (l) for each shell
    #lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    #Lmax = maximum(lvals)

    ## Offset list for each shell, used to map shell index to AO index
    #ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    #buf_arrays = [zeros(Cdouble, Lmax^4) for _ = 1:Threads.nthreads()]

    #@sync for i in 1:BS.nshells
    #    Threads.@spawn begin
    #        @inbounds begin
    #            Li = lvals[i]
    #            buf = buf_arrays[Threads.threadid()]
    #            ioff = ao_offset[i]
    #            for j in i:BS.nshells
    #                Lj = lvals[j]
    #                joff = ao_offset[j]

    #                # Call libcint
    #                libcint_1e!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

    #                # Loop through shell block and save unique elements
    #                for js = 1:Lj
    #                    J = joff + js
    #                    for is = 1:Li
    #                        I = ioff + is
    #                        J < I ? break : nothing
    #                        out[I,J] = buf[is + Li*(js-1)]
    #                        out[J,I] = out[I,J]
    #                    end
    #                end
    #            end
    #        end #inbounds
    #    end #spawn
    #end #sync
    #return out
end