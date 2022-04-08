using GaussianBasis.Libcint

export nuclear, overlap, kinetic
export ERI_2e2c, ERI_2e3c, ERI_2e4c, sparseERI_2e4c
export ERI_2e2c!

function index2(i, j)
    if i < j
        return (j * (j + 1)) >> 1 + i
    else
        return (i * (i + 1)) >> 1 + j
    end
end

function index4(i,j,k,l)
    return index2(index2(i,j), index2(k,l))
end

include("Integrals/OneElectron.jl")
include("Integrals/TwoElectronTwoCenter.jl")
include("Integrals/TwoElectronThreeCenter.jl")
include("Integrals/TwoElectronFourCenter.jl")
include("Integrals/Multipole.jl")