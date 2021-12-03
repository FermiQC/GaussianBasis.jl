atoms = Molecules.parse_string("""
C   -2.131551124300    2.286168823700    0.000000000000
H   -1.061551124300    2.286168823700    0.000000000000
H   -2.488213906200    1.408104616400    0.496683911300
H   -2.488218762100    2.295059432700   -1.008766153900
H   -2.488220057000    3.155340844300    0.512081313000""")

s = BasisFunction(0, [1/√2, 1/√2], [5.0, 1.2])
p = BasisFunction(1, [1/√2, 1/√2], [5.0, 1.2])
shells = [ [s,p], [s], [s], [s], [s]]
bset = BasisSet("Test", atoms, shells)

@testset "Misc" begin
    @test GaussianBasis.string_repr(s) == "S shell with 1 basis built from 2 primitive gaussians\n\nχ₀₀  =    0.7071067812⋅Y₀₀⋅exp(-5.0⋅r²)\n     +    0.7071067812⋅Y₀₀⋅exp(-1.2⋅r²)" 
    @test GaussianBasis.string_repr(bset) == "Test Basis Set\nNumber of shells: 6\nNumber of basis:  8\n\nC: 1s 1p \nH: 1s \nH: 1s \nH: 1s \nH: 1s"
    @test bset[1] == [s,p]
    @test bset[1,2] == p
end
