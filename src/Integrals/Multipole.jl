
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

dipole(BS::BasisSet) = get_1e_matrix(dipole!, BS, 1)
quadrupole(BS::BasisSet) = get_1e_matrix(quadrupole!, BS, 2)
octupole(BS::BasisSet) = get_1e_matrix(octupole!, BS, 3)
hexadecapole(BS::BasisSet) = get_1e_matrix(hexadecapole!, BS, 4)
