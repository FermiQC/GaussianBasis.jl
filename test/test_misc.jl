using StaticArrays

atoms = Molecules.parse_string("""
C   -2.131551124300    2.286168823700    0.000000000000
H   -1.061551124300    2.286168823700    0.000000000000
H   -2.488213906200    1.408104616400    0.496683911300
H   -2.488218762100    2.295059432700   -1.008766153900
H   -2.488220057000    3.155340844300    0.512081313000""")

shells = [
    SphericalShell(0, [1/√2, 1/√2], [5.0, 1.2], atoms[1])
    SphericalShell(1, [1/√2, 1/√2], [5.0, 1.2], atoms[1])
    SphericalShell(0, [1/√2, 1/√2], [5.0, 1.2], atoms[2])
    SphericalShell(0, [1/√2, 1/√2], [5.0, 1.2], atoms[3])
    SphericalShell(0, [1/√2, 1/√2], [5.0, 1.2], atoms[4])
    SphericalShell(0, [1/√2, 1/√2], [5.0, 1.2], atoms[5])
]

bset = BasisSet("Test", atoms, shells)

@testset "Misc" begin
    @testset "string_repr" begin
        @test GaussianBasis.string_repr(shells[1]) == "S shell with 1 basis built from 2 primitive gaussians\n\nχ₀₀  =    0.7071067812⋅Y₀₀⋅exp(-5.0⋅r²)\n     +    0.7071067812⋅Y₀₀⋅exp(-1.2⋅r²)"
        @test occursin(r"Test Basis\s+?Set\nType:\s+?Spherical", GaussianBasis.string_repr(bset))
        @test occursin(r"Number of shells:\s?+6\nNumber of basis:\s+?8", GaussianBasis.string_repr(bset))
        @test occursin(r"C: 1s 1p\s+?\nH: 1s\s+?\nH: 1s\s+?\nH: 1s\s+?\nH: 1s\s*", GaussianBasis.string_repr(bset))
    end

    @test bset[1] == bset.basis[1]

    @testset "legendre_polynomial" begin
        examples = [
            (0, 0, -0.42, +1.000000000000)
            (1, 0, +0.33, +0.330000000000)
            (1, 1, +0.33, -0.943980932011)
            (2, 0, -0.80, +0.460000000000)
            (2, 1, -0.80, +1.440000000000)
            (2, 2, -0.80, +1.080000000000)
            (3, 0, +0.70, -0.192500000000)
            (3, 1, +0.70, -1.553260683208)
            (3, 2, +0.70, +5.355000000000)
            (3, 3, +0.70, -5.463192747835)
            (4, 0, -0.30, +0.072937500000)
            (4, 1, -0.30, -1.695626930519)
            (4, 2, -0.30, -2.525250000000)
            (4, 3, -0.30, +27.344667208617)
            (4, 4, -0.30, +86.950500000000)
            (5, 0, +0.50, +0.089843750000)
            (5, 1, +0.50, +1.928259688114)
            (5, 2, +0.50, -4.921875000000)
            (5, 3, +0.50, -42.624687842515)
            (5, 4, +0.50, +265.781250000000)
            (5, 5, +0.50, -460.346628699166)
        ]

        for (l,m,x,expected) in examples
            @test GaussianBasis.legendre_polynomial(m,l,x) ≈ expected
        end
    end

    @testset "atomic_orbital_amplitude" begin
        dx = 0.07

        bset = BasisSet("cc-pvdz", "C 0.1 -0.5 -0.13\nC -0.1 0.05 0.13")

        # test s, p, and d shells
        idxs = [1, 4, 5, 25, 28]
        coef = sin.(idxs)

        b = range(-3, 3; step=dx)

        xs = [(x, y, z)[i] for i = 1:3, x in b, y in b, z in b]

        r = 0.0
        refidx = 1
        ref = atomic_orbital_amplitude(bset, refidx, xs)
        for (i, c) in zip(idxs, coef)
            r += c * sum(ref .* atomic_orbital_amplitude(bset, i, xs))
        end
        r *= dx^3

        ovlp = overlap(bset)
        @test isapprox(r, sum(c*ovlp[refidx, i] for (i,c) in zip(idxs, coef)), rtol=0.05)

        @test atomic_orbital_amplitude(bset, 1, [0,0,1]) isa Real
        @test_throws BoundsError atomic_orbital_amplitude(bset, 100, [0,0,1])
    end
end
