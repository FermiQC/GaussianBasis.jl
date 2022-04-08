function ERI_2e2c!(out::Array{Float64}, BS::BasisSet, A,i,B,j)

    # i represents the ith basis centered on A
    # j represents the jth basis centered on B

    # Get absolute indexes I and J
    I = i-1
    for b in 1:(A-1)
        I += length(BS[b])
    end
    J = j-1
    for b in 1:(B-1)
        J += length(BS[b])
    end

    Ni = 2*BS[A,i].l+1
    Nj = 2*BS[B,j].l+1

    # Check if theres enough space in 'out'
    if length(out) < Ni*Nj
        throw(BoundsError(out, Ni*Nj))
    end

    cint2c2e_sph!(out, Cint.([I,J]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
end

function ERI_2e2c(BS::BasisSet, T::DataType = Float64)

    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    buf_arrays = [zeros(Cdouble, Lmax^4) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = lvals[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset[i]
                for j in i:BS.nshells
                    Lj = lvals[j]
                    joff = ao_offset[j]

                    # Call libcint
                    cint2c2e_sph!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                    # Loop through shell block and save unique elements
                    for js = 1:Lj
                        J = joff + js
                        for is = 1:Li
                            I = ioff + is
                            J < I ? break : nothing
                            out[I,J] = buf[is + Li*(js-1)]
                            out[J,I] = out[I,J]
                        end
                    end
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end