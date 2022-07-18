function ERI_2e2c!(out, BS::BasisSet, i, j)
    generate_ERI_quartet!(out, BS.basis[i], _ghostBF, BS.basis[j], _ghostBF)
end

function ERI_2e2c!(out, BS::BasisSet{LCint}, i, j)
    cint2c2e_sph!(out, [i,j], BS.lib)
end

function ERI_2e2c(BS::BasisSet)
    out = zeros(BS.nbas, BS.nbas)
    ERI_2e2c!(out, BS)
end

function ERI_2e2c!(out, BS::BasisSet)
    get_1e_matrix!(ERI_2e2c!, out, BS)
end