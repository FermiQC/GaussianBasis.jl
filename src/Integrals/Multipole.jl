function dipole(BS::BasisSet, T::DataType = Float64)

    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas, 3)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    buf_arrays = [zeros(Cdouble, 3*Lmax^2) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = lvals[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset[i]
                for j in i:BS.nshells
                    Lj = lvals[j]
                    Lij = Li*Lj
                    joff = ao_offset[j]

                    # Call libcint
                    cint1e_r_sph!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                    I = (ioff+1):(ioff+Li)
                    J = (joff+1):(joff+Lj)

                    # Get strides for each cartesian
                    for k in 1:3
                        r = (1+Lij*(k-1)):(k*Lij)
                        @views dk = buf[r]
                        kvals = reshape(dk, Int(Li), Int(Lj))
                        out[I,J,k] .+= kvals
                        if i != j
                            out[J,I,k] .+= kvals'
                        end
                    end
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end

function quadrupole(BS::BasisSet, T::DataType = Float64)

    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas, 3, 3)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    buf_arrays = [zeros(Cdouble, 9*Lmax^2) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = lvals[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset[i]
                for j in i:BS.nshells
                    Lj = lvals[j]
                    Lij = Li*Lj
                    joff = ao_offset[j]

                    # Call libcint
                    cint1e_rr_sph!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                    I = (ioff+1):(ioff+Li)
                    J = (joff+1):(joff+Lj)

                    # Get strides for each cartesian product k1 * k2
                    for k1 in 1:3, k2 in 1:3
                        k = k2 + 3*(k1-1)
                        r = (1+Lij*(k-1)):(k*Lij)
                        @views dk = buf[r]
                        kvals = reshape(dk, Int(Li), Int(Lj))
                        out[I,J,k1,k2] .+= kvals
                        if i != j
                            out[J,I,k1,k2] .+= kvals'
                        end
                    end
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end

function octupole(BS::BasisSet, T::DataType = Float64)

    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas, 3, 3, 3)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    buf_arrays = [zeros(Cdouble, 27*Lmax^2) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = lvals[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset[i]
                for j in i:BS.nshells
                    Lj = lvals[j]
                    Lij = Li*Lj
                    joff = ao_offset[j]

                    # Call libcint
                    cint1e_rrr_sph!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                    I = (ioff+1):(ioff+Li)
                    J = (joff+1):(joff+Lj)

                    # Get strides for each cartesian product k1 * k2 * k3
                    for k1 in 1:3, k2 in 1:3, k3 in 1:3
                        k = k3 + 3*(k2-1) + 9*(k1-1)
                        r = (1+Lij*(k-1)):(k*Lij)
                        @views dk = buf[r]
                        kvals = reshape(dk, Int(Li), Int(Lj))
                        out[I,J,k1,k2,k3] .+= kvals
                        if i != j
                            out[J,I,k1,k2,k3] .+= kvals'
                        end
                    end
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end

function hexadecapole(BS::BasisSet, T::DataType = Float64)

    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas, 3, 3, 3, 3)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    buf_arrays = [zeros(Cdouble, 81*Lmax^2) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = lvals[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset[i]
                for j in i:BS.nshells
                    Lj = lvals[j]
                    Lij = Li*Lj
                    joff = ao_offset[j]

                    # Call libcint
                    cint1e_rrrr_sph!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                    I = (ioff+1):(ioff+Li)
                    J = (joff+1):(joff+Lj)

                    # Get strides for each cartesian product k1 * k2 * k3
                    for k1 in 1:3, k2 in 1:3, k3 in 1:3, k4 in 1:3
                        k = k4 + 3*(k3-1) + 9*(k2-1) + 27*(k1-1)
                        r = (1+Lij*(k-1)):(k*Lij)
                        @views dk = buf[r]
                        kvals = reshape(dk, Int(Li), Int(Lj))
                        out[I,J,k1,k2,k3,k4] .+= kvals
                        if i != j
                            out[J,I,k1,k2,k3,k4] .+= kvals'
                        end
                    end
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end