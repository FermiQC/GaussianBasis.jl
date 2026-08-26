const _ghostBF = CartesianShell(0, [1.0], [0.0], Atom(1, 1.0, [0.0, 0.0, 0.0]))

# Mutating, shell-triple-level backend for ERI_2e3c(BS1, BS2): writes into a
# caller-supplied `out` instead of allocating. Dispatches on the integral
# backend: LCint uses libcint's native 3-center kernel directly; the ACSint
# fallback below instead evaluates it as a 4-center integral against a ghost
# (zero-charge, s-type) basis function standing in for the missing 4th center.
"""
    ERI_2e3c!(out, BS::BasisSet, i, j, k)
    ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet)

Mutating counterpart of [`ERI_2e3c`](@ref): writes into the caller-supplied
`out` instead of allocating.

# Methods

  - `ERI_2e3c!(out, BS, i, j, k)`: `out` must be `(Ni,Nj,Nk)`, the `(ij|k)`
    block for shells `i,j,k` of `BS` (shell indices, not AO indices). This
    is the shell-triple primitive the full-tensor form builds on -- but it
    takes a single `BasisSet` with the regular and auxiliary shells already
    merged together, since libcint's 3-center kernel resolves shell indices
    against one basis. See [Three Centers](@ref) for how to build one and
    map shell indices into it.
  - `ERI_2e3c!(out, BS1, BS2)`: `out` must be a dense
    `BS1.nbas × BS1.nbas × BS2.nbas` array.
"""
function ERI_2e3c!(out, BS::BasisSet{LCint}, i, j, k)
    cint3c2e_sph!(out, @SVector([i,j,k]), BS.lib)
end

function ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet, i, j, k)
    generate_ERI_quartet!(out, BS1.shells[i], BS1.shells[j], BS2.shells[k], _ghostBF)
end

"""
    ERI_2e3c(BS1::BasisSet, BS2::BasisSet) -> Array{Float64,3}

Compute the full two-electron three-center integral tensor `(μν|P)`, with
`μ,ν` running over `BS1`'s AOs (the "regular" orbital basis) and `P` over
`BS2`'s (the auxiliary/fitting basis) -- the building block for density
fitting / resolution-of-the-identity approximations. Returns a dense
`BS1.nbas × BS1.nbas × BS2.nbas` array, symmetric under `μ↔ν` swap. For
repeated calls, see `ERI_2e3c!`, which writes into a preallocated array
instead of allocating.
"""
function ERI_2e3c(BS1::BasisSet, BS2::BasisSet)
    out = zeros(BS1.nbas, BS1.nbas, BS2.nbas)
    ERI_2e3c!(out, BS1, BS2)
end

function ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet)

    # NOTE: `out` is deliberately not zeroed -- the loops below cover every
    # (i<=j, k) shell triple and mirror each block over μ↔ν, so every element
    # is written and prior contents fully overwritten. Adding any skipping
    # (e.g. shell-pair screening) requires a `fill!(out, 0.0)` first. Same
    # invariant as ERI_2e4c!.

    # Pre compute number of basis per shell
    Nvals1 = num_basis.(BS1.shells)
    Nvals2 = num_basis.(BS2.shells)
    Nmax1 = maximum(Nvals1)
    Nmax2 = maximum(Nvals2)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset1 = cumsum(Nvals1) .- Nvals1
    ao_offset2 = cumsum(Nvals2) .- Nvals2

    allocate(body) = body(zeros(Cdouble, Nmax1^2*Nmax2))
    workerpool(allocate, 1:BS2.nshells; chunksize=1) do k, buf
        @inbounds begin
            Nk = Nvals2[k]
            koff = ao_offset2[k]
            for i in 1:BS1.nshells
                Ni = Nvals1[i]
                ioff = ao_offset1[i]
                for j in i:BS1.nshells
                    Nj = Nvals1[j]
                    joff = ao_offset1[j]

                    # Call libcint
                    ERI_2e3c!(buf, BS1, BS2, i, j, k)

                    # Loop through shell block and save unique elements
                    for ks = 1:Nk
                        K = koff + ks
                        for js = 1:Nj
                            J = joff + js
                            for is = 1:Ni
                                I = ioff + is
                                J < I ? break : nothing
                                out[I,J,K] = buf[is + Ni*(js-1) + Ni*Nj*(ks-1)]
                                out[J,I,K] = out[I,J,K]
                            end
                        end
                    end
                end
            end
        end #inbounds
    end
    return out
end

function ERI_2e3c!(out, BS1::BasisSet{LCint}, BS2::BasisSet{LCint}; Bmerged::Union{Nothing,BasisSet}=nothing)

    # Bmerged depends only on BS1/BS2 -- callers making several calls against
    # the same basis pair (e.g. this integral plus ∇ERI_2e3c!/∇2ERI_2e3c!
    # over every atom, which take the same keyword) can build it once and
    # pass it in rather than reconstructing it here each time. It's a small
    # fraction of the runtime but over half of this function's allocation.
    if Bmerged === nothing
        Bmerged = merge_basis(BS1, BS2)
    end

    # Pre compute number of basis per shell
    Nvals1 = num_basis.(BS1.shells)
    Nvals2 = num_basis.(BS2.shells)
    Nmax1 = maximum(Nvals1)
    Nmax2 = maximum(Nvals2)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset1 = cumsum(Nvals1) .- Nvals1
    ao_offset2 = cumsum(Nvals2) .- Nvals2


    allocate(body) = body(zeros(Cdouble, Nmax1^2*Nmax2))
    workerpool(allocate, 1:BS2.nshells; chunksize = 1) do k, buf
        @inbounds begin
            Nk = Nvals2[k]
            koff = ao_offset2[k]
            for i in 1:BS1.nshells
                Ni = Nvals1[i]
                ioff = ao_offset1[i]
                for j in i:BS1.nshells
                    Nj = Nvals1[j]
                    joff = ao_offset1[j]

                    # Call libcint
                    ERI_2e3c!(buf, Bmerged, i, j, k+BS1.nshells)

                    # Loop through shell block and save unique elements
                    for ks = 1:Nk
                        K = koff + ks
                        for js = 1:Nj
                            J = joff + js
                            for is = 1:Ni
                                I = ioff + is
                                J < I ? break : nothing
                                out[I,J,K] = buf[is + Ni*(js-1) + Ni*Nj*(ks-1)]
                                out[J,I,K] = out[I,J,K]
                            end
                        end
                    end
                end
            end
        end #inbounds
    end
    return out
end
