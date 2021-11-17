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

    #@testset "Sparse" begin
    #    for iA = 1:5
    #        idx, ∇q... = ∇sparseERI_2e4c(bs, iA)
            ## Make sure indexes are sorted by index4
    #        for k = 1:3
    #            ∇k = ∇q[k]
    #            fd_data = ∇FD_ERI_2e4c(bs, iA, k)
    #           for i in idx

                #@assert issorted(fd_idx, by = x -> GaussianBasis.index4(x[1], x[2], x[3], x[4]))
                #@test dIDX == fd_idx
            #end
        #end
    #end
end

#function compare_sparse_arrays(id1, d1, id2, d2)
#
#    l1 = length(id1)    
#    l2 = length(id2)    
#
#    i = 1
#    j = 1
#
#    while l1 > i && l2 > j
#        i4 = GaussianBasis.index4(id1[i]...)
#        j4 = GaussianBasis.index4(id2[j]...)
#        if i4 == j4
#            if d1[i] ≈ d2[j]
#                i += 1
#                j += 1
#            else
#                return false
#            end
#        elseif i4 < j4
#            if d1[i] ≈ 0
#                i += 1
#            else
#                return false    
#            end
#        elseif i4 > j4
#            if d2[i] ≈ 0
#                j += 1
#            else
#                return false    
#            end
#        end
#    end
#
#    return true
#end