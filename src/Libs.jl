# struct with backend specific data for integral computation
abstract type IntLib end
# struct flag for computations using Acsint.jl
struct ACSint <: IntLib end
# struct flag for computations using Libcint.jl
struct LCint <: IntLib 
    atm::Vector{Cint}
    natm::Cint
    bas::Vector{Cint}
    nbas::Cint
    env::Vector{Cdouble}
end

function LCint(atoms::Vector{<:Atom}, basis::Vector{<:BasisFunction})
    ATM_SLOTS = 6
    BAS_SLOTS = 8

    natm = length(atoms)

    nshells = length(basis)
    nbas    = 0
    nexps   = 0
    nprims  = 0

    for i in eachindex(atoms)
        for b in basis
            if b.atom == atoms[i]
                nbas    += num_basis(b)
                nexps   += length(b.exp)
                nprims  += length(b.coef)
            end
        end
    end

    lc_atm = zeros(Cint, natm*ATM_SLOTS)
    lc_bas = zeros(Cint, nshells*BAS_SLOTS)
    env = zeros(Cdouble, 20+4*natm+nexps+nprims)

    # Prepare the lc_atom input 
    off = 20
    ib = 0 
    for i = eachindex(atoms)
        A = atoms[i]
        # lc_atom has ATM_SLOTS (6) "spaces" for each atom
        # The first one (Z_INDEX) is the atomic number
        lc_atm[1 + ATM_SLOTS*(i-1)] = A.Z
        # The second one is the env index address for xyz
        lc_atm[2 + ATM_SLOTS*(i-1)] = off
        env[off+1:off+3] .= A.xyz ./ Molecules.bohr_to_angstrom
        off += 4 # Skip an extra slot for the kappa (nuclear model parameter)
        # The remaining 4 slots are zero.
    end

    for j = eachindex(basis)
        B = basis[j]
        # Find atom index
        a = findfirst(a->a==B.atom, atoms)
        Ne = length(B.exp)
        Nc = length(B.coef)
        # lc_bas has BAS_SLOTS for each basis set
        # The first one is the index of the atom starting from 0
        lc_bas[1 + BAS_SLOTS*ib] = a-1#i-1
        # The second one is the angular momentum
        lc_bas[2 + BAS_SLOTS*ib] = B.l
        # The third is the number of primitive functions
        lc_bas[3 + BAS_SLOTS*ib] = Nc
        # The fourth is the number of contracted functions
        lc_bas[4 + BAS_SLOTS*ib] = 1
        # The fifth is a κ parameter
        lc_bas[5 + BAS_SLOTS*ib] = 0
        # Sixth is the env index address for exponents
        lc_bas[6 + BAS_SLOTS*ib] = off
        env[off+1:off+Ne] .= B.exp
        off += Ne
        # Seventh is the env index address for contraction coeff
        lc_bas[7 + BAS_SLOTS*ib] = off
        env[off+1:off+Nc] .= B.coef
        off += Nc
        # Eigth, nothing
        ib += 1
    end
    return LCint(lc_atm, Cint(natm), lc_bas, Cint(nbas), env)
end

#function LCint(BS1::BasisSet, BS2::BasisSet)
#
#    atoms = vcat(BS1.atoms, BS2.atoms)
#    
#
#    ATM_SLOTS = 6
#    BAS_SLOTS = 8
#
#    natm = BS1.natom + BS2.natm
#    #nbas = BS1.nbas  + BS2.nbas
#    nshells = BS1.nshells + BS2.nshells
#
#    lc_atm = zeros(Cint, natm*ATM_SLOTS)
#    lc_bas = zeros(Cint, nshells*BAS_SLOTS)
#    env = zeros(Cdouble, length(BS1.lib.env) + length(BS2.lib.env) -20)
#
#    # Prepare the lc_atom input 
#    off = 20
#    ib = 0 
#    atoms = vcat(BS1.atoms, BS2.atoms)
#    for i = eachindex(atoms)
#        A = atoms[i]
#        # lc_atom has ATM_SLOTS (6) "spaces" for each atom
#        # The first one (Z_INDEX) is the atomic number
#        lc_atm[1 + ATM_SLOTS*(i-1)] = Cint(A.Z)
#        # The second one is the env index address for xyz
#        lc_atm[2 + ATM_SLOTS*(i-1)] = off
#        env[off+1:off+3] .= A.xyz ./ Molecules.bohr_to_angstrom
#        off += 4 # Skip an extra slot reserved for nuclear model
#        # The remaining 4 slots are zero.
#    end
#
#    for i = eachindex(BS1.atoms)
#        # Prepare the lc_bas input
#        for j = eachindex(BS1.basis[i])
#            B = BS1[i][j] 
#            Ne = length(B.exp)
#            Nc = length(B.coef)
#            # lc_bas has BAS_SLOTS for each basis set
#            # The first one is the index of the atom starting from 0
#            lc_bas[1 + BAS_SLOTS*ib] = i-1
#            # The second one is the angular momentum
#            lc_bas[2 + BAS_SLOTS*ib] = B.l
#            # The third is the number of primitive functions
#            lc_bas[3 + BAS_SLOTS*ib] = Nc
#            # The fourth is the number of contracted functions
#            lc_bas[4 + BAS_SLOTS*ib] = 1
#            # The fifth is a κ parameter
#            lc_bas[5 + BAS_SLOTS*ib] = 0
#            # Sixth is the env index address for exponents
#            lc_bas[6 + BAS_SLOTS*ib] = off
#            env[off+1:off+Ne] .= B.exp
#            off += Ne
#            # Seventh is the env index address for contraction coeff
#            lc_bas[7 + BAS_SLOTS*ib] = off
#            env[off+1:off+Nc] .= B.coef
#            off += Nc
#            # Eigth, nothing
#            ib += 1
#        end
#    end
#
#    for i = eachindex(BS2.atoms)
#        # Prepare the lc_bas input
#        for j = eachindex(BS2.basis[i])
#            B = BS2[i][j] 
#            Ne = length(B.exp)
#            Nc = length(B.coef)
#            # lc_bas has BAS_SLOTS for each basis set
#            # The first one is the index of the atom starting from 0
#            lc_bas[1 + BAS_SLOTS*ib] = i-1
#            # The second one is the angular momentum
#            lc_bas[2 + BAS_SLOTS*ib] = B.l
#            # The third is the number of primitive functions
#            lc_bas[3 + BAS_SLOTS*ib] = Nc
#            # The fourth is the number of contracted functions
#            lc_bas[4 + BAS_SLOTS*ib] = 1
#            # The fifth is a κ parameter
#            lc_bas[5 + BAS_SLOTS*ib] = 0
#            # Sixth is the env index address for exponents
#            lc_bas[6 + BAS_SLOTS*ib] = off
#            env[off+1:off+Ne] .= B.exp
#            off += Ne
#            # Seventh is the env index address for contraction coeff
#            lc_bas[7 + BAS_SLOTS*ib] = off
#            env[off+1:off+Nc] .= B.coef
#            off += Nc
#            # Eigth, nothing
#            ib += 1
#        end
#    end
#
#    return LCint(natom, )
#end