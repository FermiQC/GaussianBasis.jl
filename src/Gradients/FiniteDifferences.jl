# Given a list of atoms, create a cartesian displamentent in one of them (A)
function create_displacement(bs::BasisSet, A::Int, i::Int, h)

    Aplus = deepcopy(bs.atoms)
    Aplus[A].xyz[i] += h

    Aminus = deepcopy(bs.atoms)
    Aminus[A].xyz[i] -= h

    bs_plus = BasisSet(bs.name, Aplus)
    bs_minus = BasisSet(bs.name, Aminus)

    return bs_plus, bs_minus
end

function findif_ao_1e(bs::BasisSet, compute::String, A, i, h=0.005)

    bs_plus, bs_minus = create_displacement(bs, A, i, h)

    Xplus = ao_1e(bs_plus, compute)
    Xminus = ao_1e(bs_minus, compute)

    return (Xplus - Xminus) ./ (2*h)
end

function findif_ao_2e4c(bs::BasisSet, A, i, h=0.005)

    bs_plus, bs_minus = create_displacement(bs, A, i, h)

    Xplus = ao_2e4c(bs_plus)
    Xminus = ao_2e4c(bs_minus)

    return (Xplus - Xminus) ./ (2*h)

end

function findif_sparse_ao_2e4c(bs::BasisSet, A, i, h=0.005)

    bs_plus, bs_minus = create_displacement(bs, A, i, h)

    Iplus, Xplus = sparse_ao_2e4c(bs_plus)
    Iminus, Xminus = sparse_ao_2e4c(bs_minus)

    @assert Iplus == Iminus

    Xout = (Xplus - Xminus) ./ (2*h)

    return Iplus, Xout
end