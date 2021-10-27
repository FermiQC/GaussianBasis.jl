function compare(E, idx, ∇x, ∇y, ∇z)
    N = size(E, 5)
    println("Checking X")
    for i = 1:N
        f = true
        for j = i:N
            for k = 1:N
                for l = k:N
                    if GaussianBasis.index2(i-1,j-1) > GaussianBasis.index2(k-1,l-1)
                        continue
                    end
                    A = 0.0
                    try
                        Aidx = findfirst(x->x==(i-1, j-1, k-1, l-1), idx)
                        A = ∇x[Aidx]
                    catch ArgumentError
                        A = 0.0
                    end
                    B = E[1,i,j,k,l]
                    test = A ≈ B
                    if !test
                        println("$i $j $k $l")
                        f = false
                    end
                end
                !f ? break : nothing
            end
            !f ? break : nothing
        end
        println(f)
        !f ? break : nothing
    end

    println("Checking Y")
    for i = 1:N
        f = true
        for j = i:N
            for k = 1:N
                for l = k:N
                    if GaussianBasis.index2(i-1,j-1) > GaussianBasis.index2(k-1,l-1)
                        continue
                    end
                    A = 0.0
                    try
                        Aidx = findfirst(x->x==(i-1, j-1, k-1, l-1), idx)
                        A = ∇y[Aidx]
                    catch ArgumentError
                        A = 0.0
                    end
                    B = E[2,i,j,k,l]
                    test = A ≈ B
                    if !test
                        println("$i $j $k $l")
                        f = false
                    end
                end
                !f ? break : nothing
            end
            !f ? break : nothing
        end
        println(f)
        !f ? break : nothing
    end

    println("Checking Z")
    for i = 1:N
        f = true
        for j = i:N
            for k = 1:N
                for l = k:N
                    if GaussianBasis.index2(i-1,j-1) > GaussianBasis.index2(k-1,l-1)
                        continue
                    end
                    A = 0.0
                    try
                        Aidx = findfirst(x->x==(i-1, j-1, k-1, l-1), idx)
                        A = ∇z[Aidx]
                    catch ArgumentError
                        A = 0.0
                    end
                    B = E[3,i,j,k,l]
                    test = A ≈ B
                    if !test
                        println("$i $j $k $l")
                        f = false
                    end
                end
                !f ? break : nothing
            end
            !f ? break : nothing
        end
        println(f)
        !f ? break : nothing
    end
end