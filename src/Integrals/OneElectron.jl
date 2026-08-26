# Mutating, shell-pair-level forms of overlap/kinetic/nuclear(BS, i, j):
# write into a caller-supplied `out` (sized (Ni,Nj)) instead of allocating.
# Dispatches on the integral backend (LCint vs. the ACSint fallback below).
# Backend: Libcint
"""
    overlap!(out, BS::BasisSet, i, j)
    overlap!(out, BS::BasisSet)
    overlap!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    overlap!(out, BS1::BasisSet, BS2::BasisSet)

Mutating counterpart of [`overlap`](@ref): writes into the caller-supplied
`out` instead of allocating. This shell-pair, single-basis form is the
primitive every other method builds on.

# Methods

  - `overlap!(out, BS, i, j)`: `out` must be `(Ni,Nj)`, the block for shells
    `i`/`j` of `BS` (shell indices, not AO indices).
  - `overlap!(out, BS)`: `out` must be a dense `nbas × nbas` matrix.
  - `overlap!(out, BS1, BS2, i, j)`: `out` must be `(Ni,Nj)`, the block for
    shell `i` of `BS1` and shell `j` of `BS2`.
  - `overlap!(out, BS1, BS2)`: `out` must be `BS1.nbas × BS2.nbas`.
"""
function overlap!(out, BS::BasisSet{LCint}, i, j)
    cint1e_ovlp_sph!(out, @SVector([i,j]), BS.lib)
end

"""
    kinetic!(out, BS::BasisSet, i, j)
    kinetic!(out, BS::BasisSet)
    kinetic!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    kinetic!(out, BS1::BasisSet, BS2::BasisSet)

Mutating counterpart of [`kinetic`](@ref): writes into the caller-supplied
`out` instead of allocating. This shell-pair, single-basis form is the
primitive every other method builds on.

# Methods

  - `kinetic!(out, BS, i, j)`: `out` must be `(Ni,Nj)`, the block for shells
    `i`/`j` of `BS` (shell indices, not AO indices).
  - `kinetic!(out, BS)`: `out` must be a dense `nbas × nbas` matrix.
  - `kinetic!(out, BS1, BS2, i, j)`: `out` must be `(Ni,Nj)`, the block for
    shell `i` of `BS1` and shell `j` of `BS2`.
  - `kinetic!(out, BS1, BS2)`: `out` must be `BS1.nbas × BS2.nbas`.
"""
function kinetic!(out, BS::BasisSet{LCint}, i, j)
    cint1e_kin_sph!(out, @SVector([i,j]), BS.lib)
end

"""
    nuclear!(out, BS::BasisSet, i, j)
    nuclear!(out, BS::BasisSet)
    nuclear!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    nuclear!(out, BS1::BasisSet, BS2::BasisSet)

Mutating counterpart of [`nuclear`](@ref): writes into the caller-supplied
`out` instead of allocating. This shell-pair, single-basis form is the
primitive every other method builds on.

# Methods

  - `nuclear!(out, BS, i, j)`: `out` must be `(Ni,Nj)`, the block for shells
    `i`/`j` of `BS` (shell indices, not AO indices).
  - `nuclear!(out, BS)`: `out` must be a dense `nbas × nbas` matrix.
  - `nuclear!(out, BS1, BS2, i, j)`: `out` must be `(Ni,Nj)`, the block for
    shell `i` of `BS1` and shell `j` of `BS2`.
  - `nuclear!(out, BS1, BS2)`: `out` must be `BS1.nbas × BS2.nbas`.
"""
function nuclear!(out, BS::BasisSet{LCint}, i, j)
    cint1e_nuc_sph!(out, @SVector([i,j]), BS.lib)
end

# Backend: ACSint -- fall back
function overlap!(out, BS::BasisSet, i, j)
    generate_S_pair!(out, BS, i ,j)
end

function kinetic!(out, BS::BasisSet, i, j)
    generate_T_pair!(out, BS, i ,j)
end

function nuclear!(out, BS::BasisSet, i, j)
    generate_V_pair!(out, BS, i ,j)
end

# General

function overlap(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.shells[i]), num_basis(BS.shells[j]))
    overlap!(out, BS, i, j)
    return out
end

function kinetic(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.shells[i]), num_basis(BS.shells[j]))
    kinetic!(out, BS, i, j)
    return out
end

function nuclear(BS::BasisSet, i, j)
    out = zeros(num_basis(BS.shells[i]), num_basis(BS.shells[j]))
    nuclear!(out, BS, i, j)
    return out
end

"""
    overlap(BS::BasisSet) -> Matrix{Float64}
    overlap(BS::BasisSet, i, j) -> Matrix{Float64}
    overlap(BS1::BasisSet, BS2::BasisSet) -> Matrix{Float64}
    overlap(BS1::BasisSet, BS2::BasisSet, i, j) -> Matrix{Float64}

Compute the AO overlap matrix `S_{μν} = ⟨χ_μ|χ_ν⟩`.

# Methods

  - `overlap(BS)`: full, dense, symmetric `nbas × nbas` matrix for `BS`.
  - `overlap(BS, i, j)`: just the `(Ni,Nj)` block for shells `i` and `j` of
    `BS` (shell indices, not AO indices).
  - `overlap(BS1, BS2)`: mixed-basis matrix, `μ` running over `BS1` and `ν`
    over `BS2` -- a dense `BS1.nbas × BS2.nbas` matrix, not generally
    symmetric. Useful e.g. for projecting a density or orbitals from one
    basis onto another.
  - `overlap(BS1, BS2, i, j)`: mixed-basis block between shell `i` of `BS1`
    and shell `j` of `BS2`.

For repeated calls (e.g. in a hot loop or a geometry scan), see `overlap!`,
which writes into a preallocated array instead of allocating.
"""
overlap(BS::BasisSet) = get_1e_matrix(overlap!, BS)

"""
    kinetic(BS::BasisSet) -> Matrix{Float64}
    kinetic(BS::BasisSet, i, j) -> Matrix{Float64}
    kinetic(BS1::BasisSet, BS2::BasisSet) -> Matrix{Float64}
    kinetic(BS1::BasisSet, BS2::BasisSet, i, j) -> Matrix{Float64}

Compute the AO kinetic energy matrix `T`.

# Methods

  - `kinetic(BS)`: full, dense, symmetric `nbas × nbas` matrix for `BS`.
  - `kinetic(BS, i, j)`: just the `(Ni,Nj)` block for shells `i` and `j` of
    `BS` (shell indices, not AO indices).
  - `kinetic(BS1, BS2)`: mixed-basis matrix, `BS1.nbas × BS2.nbas`. See
    [`overlap`](@ref) for the general mixed-basis convention.
  - `kinetic(BS1, BS2, i, j)`: mixed-basis block between shell `i` of `BS1`
    and shell `j` of `BS2`.

For repeated calls, see `kinetic!`, which writes into a preallocated array
instead of allocating.
"""
kinetic(BS::BasisSet) = get_1e_matrix(kinetic!, BS)

"""
    nuclear(BS::BasisSet) -> Matrix{Float64}
    nuclear(BS::BasisSet, i, j) -> Matrix{Float64}
    nuclear(BS1::BasisSet, BS2::BasisSet) -> Matrix{Float64}
    nuclear(BS1::BasisSet, BS2::BasisSet, i, j) -> Matrix{Float64}

Compute the AO nuclear attraction matrix `V` (potential from every nucleus
in `BS.atoms`, summed).

# Methods

  - `nuclear(BS)`: full, dense, symmetric `nbas × nbas` matrix for `BS`.
  - `nuclear(BS, i, j)`: just the `(Ni,Nj)` block for shells `i` and `j` of
    `BS` (shell indices, not AO indices).
  - `nuclear(BS1, BS2)`: mixed-basis matrix (nuclei taken from `BS1.atoms`),
    `BS1.nbas × BS2.nbas`. See [`overlap`](@ref) for the general mixed-basis
    convention.
  - `nuclear(BS1, BS2, i, j)`: mixed-basis block between shell `i` of `BS1`
    and shell `j` of `BS2`, using the nuclei of `BS1`.

For repeated calls, see `nuclear!`, which writes into a preallocated array
instead of allocating.
"""
nuclear(BS::BasisSet) = get_1e_matrix(nuclear!, BS)

overlap!(out, BS::BasisSet) = get_1e_matrix!(overlap!, out, BS)
kinetic!(out, BS::BasisSet) = get_1e_matrix!(kinetic!, out, BS)
nuclear!(out, BS::BasisSet) = get_1e_matrix!(nuclear!, out, BS)

function get_1e_matrix(callback, BS::BasisSet)
    out = zeros(eltype(BS.atoms[1].xyz), BS.nbas, BS.nbas)
    return get_1e_matrix!(callback, out, BS)
end

# Specialized for rank-0 (overlap/kinetic/nuclear) integrals only -- each
# shell pair's (Ni,Nj) block is copied into `out` (and its (j,i) mirror,
# exploiting S_ji = S_ij^T) via a plain scalar double loop instead of
# range-indexed dot-broadcast (`out[I,J] .+= kvals`): the latter allocates
# a SubArray/Broadcasted wrapper per shell pair that Julia's escape
# analysis doesn't elide, ~90 bytes/pair, growing linearly with the number
# of shell pairs; the scalar loop is allocation-free. See
# `get_multipole_matrix!` in Multipole.jl for the rank>0 (dipole etc.)
# analogue.
function get_1e_matrix!(callback, out, BS::BasisSet)

    # Zero out the output array
    fill!(out, 0.0)

    Nvals = num_basis.(BS.shells)
    Nmax = maximum(Nvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = cumsum(Nvals) .- Nvals

    ijs = unique_ij(BS.nshells)

    if Threads.nthreads() == 1 || BS.nbas^2 < THREADING_THRESHOLD_1E
        # Below the threshold the whole matrix takes less time to build than
        # spawning the task pool costs -- see THREADING_THRESHOLD_1E.
        buf = zeros(Cdouble, Nmax^2)
        for (i,j) in ijs
            _scatter_1e!(callback, out, BS, i, j, buf, Nvals, ao_offset)
        end
    else
        allocate(body) = body(zeros(Cdouble, Nmax^2))
        workerpool(allocate, ijs; chunksize=10) do (i,j), buf
            _scatter_1e!(callback, out, BS, i, j, buf, Nvals, ao_offset)
        end
    end
    return out
end

# One shell pair: evaluate it into `buf` and scatter the (Ni,Nj) block into
# `out`, mirroring over i<->j. Shared by the serial and threaded paths above so
# the two can't drift apart.
@inline function _scatter_1e!(callback, out, BS, i, j, buf, Nvals, ao_offset)
    @inbounds begin
        Ni = Nvals[i]
        Nj = Nvals[j]
        ioff = ao_offset[i]
        joff = ao_offset[j]

        callback(buf, BS, i, j)

        for js = 1:Nj, is = 1:Ni
            v = buf[is + Ni*(js-1)]
            out[ioff+is, joff+js] = v
            if i != j
                out[joff+js, ioff+is] = v
            end
        end
    end
end

###########################################################
###########################################################
#                       Mixed basis                       #
###########################################################
###########################################################

# Mutating, shell-pair-level forms of the mixed-basis overlap/kinetic/nuclear
# above: write into a caller-supplied `out` instead of allocating.
# Backend: ACSint -- fall back
function overlap!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    generate_S_pair!(out, BS1.shells[i], BS2.shells[j])
end

function kinetic!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    generate_T_pair!(out, BS1.shells[i], BS2.shells[j])
end

function nuclear!(out, BS1::BasisSet, BS2::BasisSet, i, j)
    generate_V_pair!(out, BS1.shells[i], BS2.shells[j], BS1.atoms)
end

# General

function overlap(BS1::BasisSet, BS2::BasisSet, i, j)
    out = zeros(num_basis(BS1.shells[i]), num_basis(BS2.shells[j]))
    overlap!(out, BS1, BS2, i, j)
    return out
end

function kinetic(BS1::BasisSet, BS2::BasisSet, i, j)
    out = zeros(num_basis(BS1.shells[i]), num_basis(BS2.shells[j]))
    kinetic!(out, BS1, BS2, i, j)
    return out
end

function nuclear(BS1::BasisSet, BS2::BasisSet, i, j)
    out = zeros(num_basis(BS1.shells[i]), num_basis(BS2.shells[j]))
    nuclear!(out, BS1, BS2, i, j)
    return out
end

overlap(BS1::BasisSet, BS2::BasisSet) = get_1e_matrix(overlap!, BS1, BS2)
kinetic(BS1::BasisSet, BS2::BasisSet) = get_1e_matrix(kinetic!, BS1, BS2)
nuclear(BS1::BasisSet, BS2::BasisSet) = get_1e_matrix(nuclear!, BS1, BS2)

# `overlap!`/`kinetic!`/`nuclear!(out, BS1, BS2)` write into a caller-supplied
# `BS1.nbas × BS2.nbas` `out` instead of allocating.
overlap!(out, BS1::BasisSet, BS2::BasisSet) = get_1e_matrix!(overlap!, out, BS1, BS2)
kinetic!(out, BS1::BasisSet, BS2::BasisSet) = get_1e_matrix!(kinetic!, out, BS1, BS2)
nuclear!(out, BS1::BasisSet, BS2::BasisSet) = get_1e_matrix!(nuclear!, out, BS1, BS2)

function get_1e_matrix(callback, BS1::BasisSet, BS2::BasisSet)
    out = zeros(BS1.nbas, BS2.nbas)
    get_1e_matrix!(callback, out, BS1, BS2)
end

function get_1e_matrix!(callback, out, BS1::BasisSet, BS2::BasisSet)

    # Pre compute number of basis per shell
    Nvals1 = num_basis.(BS1.shells)
    Nvals2 = num_basis.(BS2.shells)
    Nmax1 = maximum(Nvals1)
    Nmax2 = maximum(Nvals2)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset1 = cumsum(Nvals1) .- Nvals1
    ao_offset2 = cumsum(Nvals2) .- Nvals2

    allocate(body) = body(zeros(Cdouble, Nmax1*Nmax2))
    workerpool(allocate, 1:BS1.nshells; chunksize = 1) do i, buf
        @inbounds begin
            Li = Nvals1[i]
            ioff = ao_offset1[i]
            for j in 1:BS2.nshells
                Lj = Nvals2[j]
                joff = ao_offset2[j]

                # Call libcint
                callback(buf, BS1, BS2, i, j)

                # Loop through shell block and save unique elements
                for js = 1:Lj
                    J = joff + js
                    for is = 1:Li
                        I = ioff + is
                        out[I,J] = buf[is + Li*(js-1)]
                    end
                end
            end
        end #inbounds
    end
    return out
end

function get_1e_matrix!(callback, out, BS1::BasisSet{LCint}, BS2::BasisSet{LCint})

    Bmerged = merge_basis(BS1, BS2)

    # Pre compute number of basis per shell
    Nvals1 = num_basis.(BS1.shells)
    Nvals2 = num_basis.(BS2.shells)
    Nmax1 = maximum(Nvals1)
    Nmax2 = maximum(Nvals2)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset1 = cumsum(Nvals1) .- Nvals1
    ao_offset2 = cumsum(Nvals2) .- Nvals2

    allocate(body) = body(zeros(Cdouble, Nmax1*Nmax2))
    workerpool(allocate, 1:BS1.nshells; chunksize = 1) do i, buf
        @inbounds begin
            Li = Nvals1[i]
            ioff = ao_offset1[i]
            for j in 1:BS2.nshells
                Lj = Nvals2[j]
                joff = ao_offset2[j]

                # Call libcint
                callback(buf, Bmerged, i, j+BS1.nshells)

                # Loop through shell block and save unique elements
                for js = 1:Lj
                    J = joff + js
                    for is = 1:Li
                        I = ioff + is
                        out[I,J] = buf[is + Li*(js-1)]
                    end
                end
            end
        end #inbounds
    end
    return out
end
