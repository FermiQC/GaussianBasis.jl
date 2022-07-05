# Given a list of atoms, create a cartesian displamentent in one of them (A)
function create_displacement(BS::BasisSet, A::Int, i::Int, h)

    disp = zeros(3)
    disp[i] += h
    Z, m, xyz = BS.atoms[A].Z, BS.atoms[A].mass, BS.atoms[A].xyz

    Aplus = deepcopy(BS.atoms)
    Aplus[A] = Atom(Z, m, xyz + disp)

    Aminus = deepcopy(BS.atoms)
    Aminus[A] = Atom(Z, m, xyz - disp)

    bs_plus = BasisSet(BS.name, Aplus)
    bs_minus = BasisSet(BS.name, Aminus)

    return bs_plus, bs_minus
end

function ∇FD_overlap(BS::BasisSet, A, i, h = 1e-5)
    return ∇FD_1e(BS, "overlap", A, i, h)
end

function ∇FD_kinetic(BS::BasisSet, A, i, h = 1e-5)
    return ∇FD_1e(BS, "kinetic", A, i, h)
end

function ∇FD_nuclear(BS::BasisSet, A, i, h = 1e-5)
    return ∇FD_1e(BS, "nuclear", A, i, h)
end

function ∇FD_1e(BS::BasisSet, compute::String, A, i, h)

    bs_plus, bs_minus = create_displacement(BS, A, i, h)

    Xplus = ao1e(bs_plus, compute)
    Xminus = ao1e(bs_minus, compute)

    return (Xplus - Xminus) ./ (2*h/Molecules.bohr_to_angstrom)
end

function ∇FD_ERI_2e4c(BS::BasisSet, A, i, h=1e-5)

    bs_plus, bs_minus = create_displacement(BS, A, i, h)

    Xplus = ERI_2e4c(bs_plus)
    Xminus = ERI_2e4c(bs_minus)

    return (Xplus - Xminus) ./ (2*h/Molecules.bohr_to_angstrom)

end

function ∇FD_sparseERI_2e4c(BS::BasisSet, A, i, h=1e-5)

    bs_plus, bs_minus = create_displacement(BS, A, i, h)

    Iplus, Xplus = sparseERI_2e4c(bs_plus)
    Iminus, Xminus = sparseERI_2e4c(bs_minus)

    @assert Iplus == Iminus

    Xout = (Xplus - Xminus) ./ (2*h/Molecules.bohr_to_angstrom)

    return Iplus, Xout
end

function ∇FD_ERI_2e3c(BS::BasisSet, auxBS::BasisSet, A, i, h=1e-5)

    bs_plus, bs_minus = create_displacement(BS, A, i, h)
    Abs_plus, Abs_minus = create_displacement(auxBS, A, i, h)

    Xplus  = ERI_2e3c(bs_plus, Abs_plus)
    Xminus = ERI_2e3c(bs_minus, Abs_minus)

    return (Xplus - Xminus) ./ (2*h/Molecules.bohr_to_angstrom)

end

function ∇FD_ERI_2e2c(BS::BasisSet, A, i, h=1e-5)

    bs_plus, bs_minus = create_displacement(BS, A, i, h)

    Xplus = ERI_2e2c(bs_plus)
    Xminus = ERI_2e2c(bs_minus)

    return (Xplus - Xminus) ./ (2*h/Molecules.bohr_to_angstrom)

end