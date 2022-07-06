using StaticArrays

atoms = Molecules.parse_string("""
C   -2.131551124300    2.286168823700    0.000000000000
H   -1.061551124300    2.286168823700    0.000000000000
H   -2.488213906200    1.408104616400    0.496683911300
H   -2.488218762100    2.295059432700   -1.008766153900
H   -2.488220057000    3.155340844300    0.512081313000""")

shells = [
    SphericalShell(0, SVector(1/√2, 1/√2), SVector(5.0, 1.2), atoms[1])
    SphericalShell(1, SVector(1/√2, 1/√2), SVector(5.0, 1.2), atoms[1])
    SphericalShell(0, SVector(1/√2, 1/√2), SVector(5.0, 1.2), atoms[2])
    SphericalShell(0, SVector(1/√2, 1/√2), SVector(5.0, 1.2), atoms[3])
    SphericalShell(0, SVector(1/√2, 1/√2), SVector(5.0, 1.2), atoms[4])
    SphericalShell(0, SVector(1/√2, 1/√2), SVector(5.0, 1.2), atoms[5])
]

bset = BasisSet("Test", atoms, shells)

@testset "Misc" begin
    @test GaussianBasis.string_repr(shells[1]) == "S shell with 1 basis built from 2 primitive gaussians\n\nχ₀₀  =    0.7071067812⋅Y₀₀⋅exp(-5.0⋅r²)\n     +    0.7071067812⋅Y₀₀⋅exp(-1.2⋅r²)" 
    @test occursin(r"Test Basis\s+?Set\nType:\s+?Spherical", GaussianBasis.string_repr(bset))
    @test occursin(r"Number of shells:\s?+6\nNumber of basis:\s+?8", GaussianBasis.string_repr(bset))
    @test occursin(r"C: 1s 1p\s+?\nH: 1s\s+?\nH: 1s\s+?\nH: 1s\s+?\nH: 1s\s*", GaussianBasis.string_repr(bset))
    @test bset[1] == bset.basis[1]
end
