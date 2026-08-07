# Finite-difference reference for the analytic Hessian-integral functions,
# mirroring Gradients/FiniteDifferences.jl one derivative order up: central
# difference of the already-trusted, already-tested first-derivative
# functions (∇overlap, ∇kinetic) rather than differencing raw energies/
# integrals twice, so this only ever validates the *new* Hessian code against
# machinery that's already independently verified.

function ∇2FD_1e(BS::BasisSet, gradcallback, iA, iB, k, h=1e-5)
    bs_plus, bs_minus = create_displacement(BS, iB, k, h)
    Xplus = gradcallback(bs_plus, iA)
    Xminus = gradcallback(bs_minus, iA)
    return (Xplus .- Xminus) ./ (2*h/Molecules.bohr_to_angstrom)
end

function ∇2FD_overlap(BS::BasisSet, iA, iB, k, h=1e-5)
    return ∇2FD_1e(BS, ∇overlap, iA, iB, k, h)
end

function ∇2FD_kinetic(BS::BasisSet, iA, iB, k, h=1e-5)
    return ∇2FD_1e(BS, ∇kinetic, iA, iB, k, h)
end

function ∇2FD_nuclear(BS::BasisSet, iA, iB, k, h=1e-5)
    return ∇2FD_1e(BS, ∇nuclear, iA, iB, k, h)
end
