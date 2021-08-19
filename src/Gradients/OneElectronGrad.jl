function grad_ao_1e(BS::BasisSet, compute::String, iA, T::DataType = Float64)

    if compute == "overlap"
        libcint_1e! =  cint1e_ipovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_ipkin_sph!
    elseif compute == "nuclear"
        return grad_ao_nuc(BS, iA)
    end

    A = BS.atoms[iA]

    # Pre allocate output
    out = zeros(T, 3, BS.nbas, BS.nbas)

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
            Lij = Li*Lj
            # Call libcint
            libcint_1e!(buf, Cint.([i,j]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            I = (ioff+1):(ioff+Li)
            J = (joff+1):(joff+Lj)

            # Get strides for each cartesian
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] .= -reshape(∇k, Int(Li), Int(Lj))
            end
            out[:,J,I] .= permutedims(out[:,I,J], (1,3,2))
        end
    end
    return out
end

# Needs improvement!!!!
function grad_ao_nuc(BS::BasisSet, iA, T::DataType = Float64)

    A = BS.atoms[iA]

    # Fudge lc_atoms
    lc_atoms_mod = deepcopy(BS.lc_atoms)
    for i = eachindex(BS.atoms)
        if i == iA
            continue 
        end
        lc_atoms_mod[1 + 6*(i-1)] = 0
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
    
    # Pre allocate output
    out = zeros(T, 3, BS.nbas, BS.nbas)

    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]
    Lmax = maximum(lvals)
    buf = zeros(Cdouble, 3*Lmax^2)

    #Case 3
    for i in notAshells
        Li = lvals[i+1]
        ioff = ao_offset[i+1]
        I = (ioff+1):(ioff+Li)
        for j in notAshells
            Lj = lvals[j+1]
            Lij = Li*Lj
            joff = ao_offset[j+1]
            J = (joff+1):(joff+Lj)

            # + ⟨i'|Va|j⟩ + ⟨i|Va|j'⟩   (Note that Va is the potential of the nuclei A alone!!)
            cint1e_ipnuc_sph!(buf, Cint.([i,j]), lc_atoms_mod, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            # Get strides for each cartesian
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] .= reshape(∇k, Int(Li), Int(Lj))
            end

            cint1e_ipnuc_sph!(buf, Cint.([j,i]), lc_atoms_mod, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] += transpose(reshape(∇k, Int(Lj), Int(Li)))
            end
        end
    end

    #Case 2
    for i in Ashells
        Li = lvals[i+1]
        ioff = ao_offset[i+1]
        I = (ioff+1):(ioff+Li)
        for j in Ashells
            Lj = lvals[j+1]
            Lij = Li*Lj
            joff = ao_offset[j+1]
            J = (joff+1):(joff+Lj)

            # - ⟨i'|∑Vc|j⟩ - ⟨i|∑Vc|j'⟩
            # + ⟨i'|Va|j⟩ + ⟨i|Va|j'⟩   (Note that Va is the potential of the nuclei A alone!!)
            cint1e_ipnuc_sph!(buf, Cint.([i,j]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            #out[:,I,J] .= -reshape(buf[1:(3*Li*Lj)], 3, Int(Li), Int(Lj))
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] -= reshape(∇k, Int(Li), Int(Lj))
            end

            cint1e_ipnuc_sph!(buf, Cint.([i,j]), lc_atoms_mod, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            #out[:,I,J] += reshape(buf[1:(3*Li*Lj)], 3, Int(Li), Int(Lj))
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] += reshape(∇k, Int(Li), Int(Lj))
            end

            cint1e_ipnuc_sph!(buf, Cint.([j,i]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            #out[:,I,J] -= permutedims(reshape(buf[1:(3*Li*Lj)], 3, Int(Lj), Int(Li)), (1,3,2))
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] -= transpose(reshape(∇k, Int(Lj), Int(Li)))
            end

            cint1e_ipnuc_sph!(buf, Cint.([j,i]), lc_atoms_mod, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            #out[:,I,J] += permutedims(reshape(buf[1:(3*Li*Lj)], 3, Int(Lj), Int(Li)), (1,3,2))
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] += transpose(reshape(∇k, Int(Lj), Int(Li)))
            end
        end
    end

    #Case 1
    for i in Ashells
        Li = lvals[i+1]
        ioff = ao_offset[i+1]
        I = (ioff+1):(ioff+Li)
        for j in notAshells
            Lj = lvals[j+1]
            Lij = Li*Lj
            joff = ao_offset[j+1]
            J = (joff+1):(joff+Lj)

            # - ⟨i'|∑Vc|j⟩ + ⟨i'|Va|j⟩ + ⟨i|Va|j'⟩
            cint1e_ipnuc_sph!(buf, Cint.([i,j]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            #out[:,I,J] .= -reshape(buf[1:(3*Li*Lj)], 3, Int(Li), Int(Lj))
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] -= reshape(∇k, Int(Li), Int(Lj))
            end

            cint1e_ipnuc_sph!(buf, Cint.([i,j]), lc_atoms_mod, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            #out[:,I,J] += reshape(buf[1:(3*Li*Lj)], 3, Int(Li), Int(Lj))
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] += reshape(∇k, Int(Li), Int(Lj))
            end

            cint1e_ipnuc_sph!(buf, Cint.([j,i]), lc_atoms_mod, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            #out[:,I,J] += permutedims(reshape(buf[1:(3*Li*Lj)], 3, Int(Lj), Int(Li)), (1,3,2))
            for k in 1:3
                r = (1+Lij*(k-1)):(k*Lij)
                ∇k = buf[r]
                out[k,I,J] += transpose(reshape(∇k, Int(Lj), Int(Li)))
            end

            out[:,J,I] .= permutedims(out[:,I,J], (1,3,2))
        end
    end
    return out
end