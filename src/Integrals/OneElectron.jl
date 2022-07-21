# Backend: Libcint
function overlap!(out, BS::BasisSet{LCint}, i, j)
    cint1e_ovlp_sph!(out, [i,j], BS.lib)
end

function kinetic!(out, BS::BasisSet{LCint}, i, j)
    cint1e_kin_sph!(out, [i,j], BS.lib)
end

function nuclear!(out, BS::BasisSet{LCint}, i, j)
    cint1e_nuc_sph!(out, [i,j], BS.lib)
end

# Backend: ACSint -- fall back 
function overlap!(out, BS::BasisSet, i, j)
    generate_S_pair!(out, BS, i ,j)
end

function kinetic!(out, BS::BasisSet, i, j)
    generate_T_pair!(out, BS, i ,j)
end

function nuclear!(out, BS::BasisSet, i, j)
    generate_V_pair!(out, BS, i ,j)
end

# General
function overlap(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.basis[i]), num_basis(BS.basis[j]))
    overlap!(out, BS, i, j)
    return out
end

function kinetic(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.basis[i]), num_basis(BS.basis[j]))
    kinetic!(out, BS, i, j)
    return out
end

function nuclear(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.basis[i]), num_basis(BS.basis[j]))
    nuclear!(out, BS, i, j)
    return out
end

overlap(BS::BasisSet) = get_1e_matrix(overlap!, BS)
kinetic(BS::BasisSet) = get_1e_matrix(kinetic!, BS)
nuclear(BS::BasisSet) = get_1e_matrix(nuclear!, BS)
overlap!(out, BS::BasisSet) = get_1e_matrix!(overlap!, out, BS)
kinetic!(out, BS::BasisSet) = get_1e_matrix!(kinetic!, out, BS)
nuclear!(out, BS::BasisSet) = get_1e_matrix!(nuclear!, out, BS)

function get_1e_matrix(callback, BS::BasisSet)
    out = zeros(eltype(BS.atoms[1].xyz), BS.nbas, BS.nbas)
    get_1e_matrix!(callback, out, BS)
end

function get_1e_matrix!(callback, out, BS::BasisSet)

    # Pre compute number of basis per shell
    Nvals = num_basis.(BS.basis)
    Nmax = maximum(Nvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:BS.nshells]

    buf_arrays = [zeros(eltype(BS.atoms[1].xyz), Nmax^2) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = Nvals[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset[i]
                for j in i:BS.nshells
                    Lj = Nvals[j]
                    joff = ao_offset[j]

                    # Call libcint
                    callback(buf, BS, i, j)

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

########################################################### 
########################################################### 
#                       Mixed basis                       #
########################################################### 
########################################################### 

# Backend: ACSint -- fall back
function overlap!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    generate_S_pair!(out, BS1.basis[i], BS2.basis[j])
end

function kinetic!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    generate_T_pair!(out, BS1.basis[i], BS2.basis[j])
end

function nuclear!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    generate_V_pair!(out, BS1.basis[i], BS2.basis[j], BS1.atoms)
end

# General
function overlap(BS1::BasisSet, BS2::BasisSet, i, j)
    out = zeros(num_basis(BS1.basis[i]), num_basis(BS2.basis[j]))
    overlap!(out, BS1, BS2, i, j)
    return out
end

function kinetic(BS1::BasisSet, BS2::BasisSet, i, j)
    out = zeros(num_basis(BS1.basis[i]), num_basis(BS2.basis[j]))
    kinetic!(out, BS1, BS2, i, j)
    return out
end

function nuclear(BS1::BasisSet, BS2::BasisSet, i, j)
    out = zeros(num_basis(BS1.basis[i]), num_basis(BS2.basis[j]))
    nuclear!(out, BS1, BS2, i, j)
    return out
end

overlap(BS1::BasisSet, BS2::BasisSet) = get_1e_matrix(overlap!, BS1, BS2)
kinetic(BS1::BasisSet, BS2::BasisSet) = get_1e_matrix(kinetic!, BS1, BS2)
nuclear(BS1::BasisSet, BS2::BasisSet) = get_1e_matrix(nuclear!, BS1, BS2)
overlap!(out, BS1::BasisSet, BS2::BasisSet) = get_1e_matrix!(overlap!, out, BS1, BS2)
kinetic!(out, BS1::BasisSet, BS2::BasisSet) = get_1e_matrix!(kinetic!, out, BS1, BS2)
nuclear!(out, BS1::BasisSet, BS2::BasisSet) = get_1e_matrix!(nuclear!, out, BS1, BS2)

function get_1e_matrix(callback, BS1::BasisSet, BS2::BasisSet)
    out = zeros(BS1.nbas, BS2.nbas)
    get_1e_matrix!(callback, out, BS1, BS2)
end

function get_1e_matrix!(callback, out, BS1::BasisSet, BS2::BasisSet)

    # Pre compute number of basis per shell
    Nvals1 = num_basis.(BS1.basis)
    Nvals2 = num_basis.(BS2.basis)
    Nmax1 = maximum(Nvals1)
    Nmax2 = maximum(Nvals2)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset1 = [sum(Nvals1[1:(i-1)]) for i = 1:BS1.nshells]
    ao_offset2 = [sum(Nvals2[1:(i-1)]) for i = 1:BS2.nshells]

    buf_arrays = [zeros(Cdouble, Nmax1*Nmax2) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS1.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = Nvals1[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset1[i]
                for j in 1:BS2.nshells
                    Lj = Nvals2[j]
                    joff = ao_offset2[j]

                    # Call libcint
                    callback(buf, BS1, BS2, i, j)

                    # Loop through shell block and save unique elements
                    for js = 1:Lj
                        J = joff + js
                        for is = 1:Li
                            I = ioff + is
                            out[I,J] = buf[is + Li*(js-1)]
                        end
                    end
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end

function get_1e_matrix!(callback, out, BS1::BasisSet{LCint}, BS2::BasisSet{LCint})

    atoms = unique(vcat(BS1.atoms, BS2.atoms))
    basis = vcat(BS1.basis, BS2.basis)

    Bmerged = BasisSet("$(BS1.name*BS2.name)", atoms, basis)

    # Pre compute number of basis per shell
    Nvals1 = num_basis.(BS1.basis)
    Nvals2 = num_basis.(BS2.basis)
    Nmax1 = maximum(Nvals1)
    Nmax2 = maximum(Nvals2)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset1 = [sum(Nvals1[1:(i-1)]) for i = 1:BS1.nshells]
    ao_offset2 = [sum(Nvals2[1:(i-1)]) for i = 1:BS2.nshells]

    buf_arrays = [zeros(Cdouble, Nmax1*Nmax2) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS1.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = Nvals1[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset1[i]
                for j in 1:BS2.nshells
                    Lj = Nvals2[j]
                    joff = ao_offset2[j]

                    # Call libcint
                    callback(buf, Bmerged, i, j+BS1.nshells)

                    # Loop through shell block and save unique elements
                    for js = 1:Lj
                        J = joff + js
                        for is = 1:Li
                            I = ioff + is
                            out[I,J] = buf[is + Li*(js-1)]
                        end
                    end
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end