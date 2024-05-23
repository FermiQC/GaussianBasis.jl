function gen_multipole(BS::BasisSet, callback, rank::Integer)
    out = zeros(BS.nbas, BS.nbas, Iterators.repeated(3, rank)...)

    Nvals = num_basis.(BS.basis)
    Nmax = maximum(Nvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:BS.nshells]

    buf_arrays = [zeros(Cdouble, 3^rank * Nmax^2) for _ = 1:Threads.nthreads()]

    Threads.@threads :static for i in 1:BS.nshells
        Ni = Nvals[i]
        buf = buf_arrays[Threads.threadid()]
        ioff = ao_offset[i]
        for j in i:BS.nshells
            Nj = Nvals[j]
            Nij = Ni*Nj
            joff = ao_offset[j]

            callback(buf, BS, i, j)
            I = (ioff+1):(ioff+Ni)
            J = (joff+1):(joff+Nj)

            # Get strides for each cartesian product
            for (n, ks) in enumerate(Iterators.product(Iterators.repeated(1:3, rank)...))
                r = Nij*(n-1)+1:n*Nij
                kvals = reshape(view(buf, r), Ni, Nj)
                out[I,J,ks...] .+= kvals
                if i != j
                    out[J,I,ks...] .+= kvals'
                end
            end
        end
    end
    return out
end

dipole!(out, BS::BasisSet{LCint}, i, j) = cint1e_r_sph!(out, [i,j], BS.lib)
quadrupole!(out, BS::BasisSet{LCint}, i, j) = cint1e_rr_sph!(out, [i,j], BS.lib)
octupole!(out, BS::BasisSet{LCint}, i, j) = cint1e_rrr_sph!(out, [i,j], BS.lib)
hexadecapole!(out, BS::BasisSet{LCint}, i, j) = cint1e_rrrr_sph!(out, [i,j], BS.lib)

function dipole(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.basis[i]), num_basis(BS.basis[j]), 3)
    dipole!(out, BS, i, j)
    return out
end
function quadrupole(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.basis[i]), num_basis(BS.basis[j]), 3, 3)
    quadrupole!(out, BS, i, j)
    return out
end
function octupole(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.basis[i]), num_basis(BS.basis[j]), 3, 3, 3)
    octupole!(out, BS, i, j)
    return out
end
function hexadecapole(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.basis[i]), num_basis(BS.basis[j]), 3, 3, 3, 3)
    hexadecapole!(out, BS, i, j)
    return out
end

dipole(BS::BasisSet) = gen_multipole(BS, dipole!, 1)
quadrupole(BS::BasisSet) = gen_multipole(BS, quadrupole!, 2)
octupole(BS::BasisSet) = gen_multipole(BS, octupole!, 3)
hexadecapole(BS::BasisSet) = gen_multipole(BS, hexadecapole!, 4)
