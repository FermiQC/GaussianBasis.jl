# Test if elements in the sparse array (S) are present
# in the dense array (D)
function SinD(S, D, verbose = false)
    Si = S[1]
    Sx = S[2]
    Sy = S[3]
    Sz = S[4]

    Dxs = D[:,:,:,:,1]
    Dys = D[:,:,:,:,2]
    Dzs = D[:,:,:,:,3]
    for i in eachindex(Si)
        idx = Si[i]
        x = Sx[i]
        y = Sy[i]
        z = Sz[i]
        Dx = Dxs[idx...]
        Dy = Dys[idx...]
        Dz = Dzs[idx...]

        for (k,Dk) in [(x,Dx), (y,Dy), (z,Dz)]
            if abs(k - Dk) > 1e-10
                println("$idx diff found: $k against $Dk")
                return false
            end
        end
        verbose ? println("$idx OK") : nothing
    end
    return true
end

atoms = Molecules.parse_string("""
C   -2.131551124300    2.286168823700    0.000000000000
H   -1.061551124300    2.286168823700    0.000000000000
H   -2.488213906200    1.408104616400    0.496683911300
H   -2.488218762100    2.295059432700   -1.008766153900
H   -2.488220057000    3.155340844300    0.512081313000""")

bs = BasisSet("cc-pvdz", atoms)

@testset "∂S/∂X" begin
    for iA = 1:5
        dS = ∇overlap(bs, iA)
        for k = 1:3
            @test dS[:,:,k] ≈ ∇FD_overlap(bs, iA, k)
        end
    end
end

@testset "∂T/∂X" begin
    for iA = 1:5
        dT = ∇kinetic(bs, iA)
        for k = 1:3
            @test dT[:,:,k] ≈ ∇FD_kinetic(bs, iA, k)
        end
    end
end

@testset "∂V/∂X" begin
    for iA = 1:5
        dV = ∇nuclear(bs, iA)
        for k = 1:3
            @test dV[:,:,k] ≈ ∇FD_nuclear(bs, iA, k)
        end
    end
end

@testset "∂[ij|kl]/∂X" begin
    @testset "Dense" begin
        for iA = 1:5
            dERI = ∇ERI_2e4c(bs, iA)
            for k = 1:3
                @test dERI[:,:,:,:,k] ≈ ∇FD_ERI_2e4c(bs, iA, k)
            end
        end
    end

    @testset "Sparse" begin
        for iA = 1:5
            sparse = ∇sparseERI_2e4c(bs, iA)
            dense = ∇ERI_2e4c(bs, iA)
            @test SinD(sparse, dense) 
        end
    end
end