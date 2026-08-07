atoms = Molecules.parse_string("""
C   -2.131551124300    2.286168823700    0.000000000000
H   -1.061551124300    2.286168823700    0.000000000000
H   -2.488213906200    1.408104616400    0.496683911300
H   -2.488218762100    2.295059432700   -1.008766153900
H   -2.488220057000    3.155340844300    0.512081313000""")

bs = BasisSet("cc-pvdz", atoms)

@testset "∂²S/∂X²" begin
    for iA = 1:length(atoms), iB = 1:length(atoms)
        d2S = ∇2overlap(bs, iA, iB)
        for k = 1:3
            @test d2S[:,:,:,k] ≈ ∇2FD_overlap(bs, iA, iB, k) atol=1e-6
        end
    end
end

@testset "∂²T/∂X²" begin
    for iA = 1:length(atoms), iB = 1:length(atoms)
        d2T = ∇2kinetic(bs, iA, iB)
        for k = 1:3
            @test d2T[:,:,:,k] ≈ ∇2FD_kinetic(bs, iA, iB, k) atol=1e-6
        end
    end
end

@testset "∂²V/∂X²" begin
    for iA = 1:length(atoms), iB = 1:length(atoms)
        d2V = ∇2nuclear(bs, iA, iB)
        for k = 1:3
            @test d2V[:,:,:,k] ≈ ∇2FD_nuclear(bs, iA, iB, k) atol=1e-6
        end
    end
end

@testset "∂²V/∂X² translational invariance" begin
    # sum_B H(A,B) == 0 exactly (to machine precision), for any fixed A --
    # the whole molecule (shells + all nuclear charges) is translation
    # invariant, at every derivative order. Independent of finite
    # differences, so it catches a different class of bug (global sign/axis
    # errors that a per-entry FD comparison can also miss if consistently
    # applied on both sides).
    for iA = 1:length(atoms)
        total = zeros(bs.nbas, bs.nbas, 3, 3)
        for iB = 1:length(atoms)
            total .+= ∇2nuclear(bs, iA, iB)
        end
        @test maximum(abs.(total)) < 1e-8
    end
end
