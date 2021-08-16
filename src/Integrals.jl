using GaussianBasis.Libcint

function index2(i, j)
    if i < j
        return (j * (j + 1)) >> 1 + i
    else
        return (i * (i + 1)) >> 1 + j
    end
end

include("Integrals/OneElectron.jl")
include("Integrals/TwoElectronTwoCenter.jl")
include("Integrals/TwoElectronThreeCenter.jl")
include("Integrals/TwoElectronFourCenter.jl")