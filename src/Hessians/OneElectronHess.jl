# --- Bare libcint primitives -------------------------------------------
#
# A Hessian needs TWO kernels per integral, not one: `ipip` puts both
# derivatives on the same shell center, `ipXip` puts one on each. Unlike the
# gradients -- where every case was the single `_μ!` primitive with its
# arguments permuted -- these are genuinely different libcint kernels and
# neither can be expressed as an argument permutation of the other. Same
# situation as ∇ERI_2e3c_μ!/∇ERI_2e3c_P!.
#
# Note there is no sign flip here, unlike the gradient primitives. libcint
# differentiates w.r.t. the electron coordinate and a GTO depends on r and
# its center only through (r-R), so each derivative contributes one factor of
# -1 -- and a second derivative applies two, which cancel.

"""
    ∇2overlap_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    ∇2kinetic_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    ∇2nuclear_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    ∇2ERI_2e2c_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int)

Libcint call placing BOTH derivatives on shell `i` (the "μ" shell), for
shells `i,j` of `BS`. `out` comes back as a raw `(Ni,Nj,3,3)` block and may
be a contiguous view.

`∇2nuclear_μμ!` covers only the piece of the nuclear-attraction Hessian
where both derivatives land on shell centers, summed over every nucleus;
the terms where a derivative lands on a nuclear *charge* position are added
separately by [`∇2nuclear!`](@ref).

> No bounds checking, no zero-block skipping, no output-size validation.
"""
function ∇2overlap_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint1e_ipipovlp_sph!(out, @SVector(Cint[i-1, j-1]),
                         lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

"""
    ∇2overlap_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    ∇2kinetic_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    ∇2nuclear_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    ∇2ERI_2e2c_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int)

Libcint call placing one derivative on shell `i` and one on shell `j`, for
shells `i,j` of `BS`. `out` comes back as a raw `(Ni,Nj,3,3)` block and may
be a contiguous view.

!!! warning "Derivative-axis order"
    The two derivative-component axes come back in **(ket, bra)** order:
    `out[:,:,q,p]` is `d²X / d(shell i)_p d(shell j)_q`. For overlap and
    kinetic this is invisible -- those depend only on the shell-shell
    separation, so the cross-Hessian is symmetric under `p<->q` -- but it is
    *not* invisible in general (nuclear attraction's shell-only piece has a
    fixed background charge distribution that breaks the symmetry), so
    callers must orient it deliberately.

> No bounds checking, no zero-block skipping, no output-size validation.
"""
function ∇2overlap_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint1e_ipovlpip_sph!(out, @SVector(Cint[i-1, j-1]),
                         lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2overlap_μμ!)
function ∇2kinetic_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint1e_ipipkin_sph!(out, @SVector(Cint[i-1, j-1]),
                        lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2overlap_μν!)
function ∇2kinetic_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint1e_ipkinip_sph!(out, @SVector(Cint[i-1, j-1]),
                        lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2overlap_μμ!)
function ∇2nuclear_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint1e_ipipnuc_sph!(out, @SVector(Cint[i-1, j-1]),
                        lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2overlap_μν!)
function ∇2nuclear_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint1e_ipnucip_sph!(out, @SVector(Cint[i-1, j-1]),
                        lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2overlap_μμ!)
function ∇2ERI_2e2c_μμ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint2c2e_ipip1_sph!(out, @SVector(Cint[i-1, j-1]),
                        lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

@doc (@doc ∇2overlap_μν!)
function ∇2ERI_2e2c_μν!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint2c2e_ip1ip2_sph!(out, @SVector(Cint[i-1, j-1]),
                         lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    return out
end

# Shell-pair Hessian block for shells (i,j) w.r.t. atoms (A,B), accumulated
# into `out` (which must be (Ni,Nj,3,3) and is NOT zeroed here -- the dense
# loop and the public shell-pair wrappers each handle that). `μμ!`/`μν!` are
# the two primitives for the integral being differentiated.
#
#   X_i,X_j = shell i,j sits on atom A;  Y_i,Y_j = ... on atom B
#
#   X_i && Y_i : both derivatives on shell i          -> μμ(i,j)
#   X_i && Y_j : d/dA on i, d/dB on j                 -> μν(i,j), axes swapped
#   X_j && Y_i : d/dA on j, d/dB on i                 -> μν(i,j), as-is
#   X_j && Y_j : both derivatives on shell j          -> μμ(j,i), AO axes swapped
#
# Every permutation is folded into the index expression of a scalar loop, so
# no permutedims temporaries are needed (permutedims! allocates ~384 B per
# call even though it is the in-place form).
@inline function _∇21e_pair!(μμ!, μν!, out, BS, i, j, Ni, Nj,
                             X_i, X_j, Y_i, Y_j, buf)
    Nij = Ni*Nj
    bufv = view(buf, 1:9*Nij)

    if X_i && Y_i
        μμ!(bufv, BS, i, j)                       # raw (Ni,Nj,3,3)
        @inbounds for q = 1:3, p = 1:3
            o = Nij*(p-1) + 3*Nij*(q-1)
            for b = 1:Nj, a = 1:Ni
                out[a,b,p,q] += bufv[o + a + Ni*(b-1)]
            end
        end
    end

    if (X_i && Y_j) || (X_j && Y_i)
        μν!(bufv, BS, i, j)                       # raw (Ni,Nj,3,3), axes (ket,bra)
        if X_i && Y_j
            # want d2/d(i=A)_p d(j=B)_q = raw[:,:,q,p]
            @inbounds for q = 1:3, p = 1:3
                o = Nij*(q-1) + 3*Nij*(p-1)
                for b = 1:Nj, a = 1:Ni
                    out[a,b,p,q] += bufv[o + a + Ni*(b-1)]
                end
            end
        end
        if X_j && Y_i
            # want d2/d(j=A)_p d(i=B)_q = raw[:,:,p,q]
            @inbounds for q = 1:3, p = 1:3
                o = Nij*(p-1) + 3*Nij*(q-1)
                for b = 1:Nj, a = 1:Ni
                    out[a,b,p,q] += bufv[o + a + Ni*(b-1)]
                end
            end
        end
    end

    if X_j && Y_j
        μμ!(bufv, BS, j, i)                       # raw (Nj,Ni,3,3)
        @inbounds for q = 1:3, p = 1:3
            o = Nij*(p-1) + 3*Nij*(q-1)
            for b = 1:Nj, a = 1:Ni
                out[a,b,p,q] += bufv[o + b + Nj*(a-1)]
            end
        end
    end
    return out
end

function ∇21e(BS::BasisSet, compute::String, iA, iB)
    out = zeros(BS.nbas, BS.nbas, 3, 3)
    return ∇21e!(out, BS, compute, iA, iB)
end

# Second derivative of overlap/kinetic w.r.t. two nuclear positions (atoms iA,
# iB). Unlike nuclear attraction, these integrals only depend on two
# positions total (the two shell centers), so -- for a given shell pair (i,j)
# with atom(i)=Ai, atom(j)=Aj -- the only combinations of (dR_iA, dR_iB) that
# can be nonzero are built from whichever of shell i's/shell j's centers
# happen to coincide with atom iA and/or atom iB:
#   - both derivatives land on shell i (Ai==iA and Ai==iB, i.e. iA==iB==Ai):
#     same-center kernel (ipip) on shell i, i.e. cint1e_ipip*_sph!([i,j],...)
#   - one derivative on shell i, one on shell j (cross terms): cross-center
#     kernel (ipXip) on [i,j]
#   - both derivatives land on shell j: same-center kernel on shell j,
#     obtained by calling it with shells reversed ([j,i]) and transposing
#     back into (i,j) AO ordering
# When iA==iB and both shells sit on that atom, all four combinations apply
# simultaneously and add together (chain rule for a doubled variable), which
# fixes the size of the diagonal Hessian block.
# Maps the legacy `compute::String` selector onto the two primitives. Kept so
# existing callers (and ∇2nuclear!'s "nuclear_shellonly" request) keep working;
# new code should call the ∇2overlap!/∇2kinetic!/∇2ERI_2e2c! wrappers, which
# name the primitives directly.
function _∇21e_kernels(compute::String)
    compute == "overlap"           && return (∇2overlap_μμ!,   ∇2overlap_μν!)
    compute == "kinetic"           && return (∇2kinetic_μμ!,   ∇2kinetic_μν!)
    compute == "nuclear_shellonly" && return (∇2nuclear_μμ!,   ∇2nuclear_μν!)
    compute == "metric"            && return (∇2ERI_2e2c_μμ!,  ∇2ERI_2e2c_μν!)
    throw(ArgumentError("compute must be \"overlap\", \"kinetic\", \"nuclear_shellonly\", or \"metric\" (use ∇2nuclear! for the full nuclear attraction Hessian)"))
end

∇21e!(out, BS::BasisSet, compute::String, iA, iB) =
    ∇21e!(_∇21e_kernels(compute)..., out, BS, iA, iB)

# Dense whole-matrix Hessian, parametrized by the two primitives rather than
# by a string -- mirroring how ∇1e! in Gradients/OneElectronGrad.jl takes a
# callback. Serial, like the gradient drivers: callers wanting parallelism
# should thread over atom pairs, which are independent.
function ∇21e!(μμ!, μν!, out, BS::BasisSet, iA, iB)

    if size(out) != (BS.nbas, BS.nbas, 3, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3, 3)"))
    end
    fill!(out, 0.0)

    Aat = BS.atoms[iA]
    Bat = BS.atoms[iB]

    Nvals = num_basis.(BS.shells)
    ao_offset = cumsum(Nvals) .- Nvals
    Nmax = maximum(Nvals)
    buf = Vector{Cdouble}(undef, 9*Nmax^2)
    block = Vector{Cdouble}(undef, 9*Nmax^2)

    @inbounds for i in 1:BS.nshells
        atom_i = BS.shells[i].atom
        X_i = (atom_i == Aat)
        Y_i = (atom_i == Bat)
        # No early-exit on X_i||Y_i alone: shell i not touching A or B is
        # still relevant if shell j does.
        Ni = Nvals[i]
        ioff = ao_offset[i]

        for j in 1:BS.nshells
            atom_j = BS.shells[j].atom
            X_j = (atom_j == Aat)
            Y_j = (atom_j == Bat)

            (X_i || X_j) || continue
            (Y_i || Y_j) || continue

            Nj = Nvals[j]
            joff = ao_offset[j]

            blk = reshape(view(block, 1:9*Ni*Nj), Ni, Nj, 3, 3)
            fill!(blk, 0.0)
            _∇21e_pair!(μμ!, μν!, blk, BS, i, j, Ni, Nj, X_i, X_j, Y_i, Y_j, buf)

            for q = 1:3, p = 1:3, b = 1:Nj, a = 1:Ni
                out[ioff+a, joff+b, p, q] += blk[a,b,p,q]
            end
        end
    end

    return out
end

"""
    ∇2overlap(BS::BasisSet, iA, iB) -> Array{Float64,4}

Second derivative (Hessian) of the AO overlap matrix `S` w.r.t. atoms
`iA`,`iB`'s three Cartesian coordinates each, `∂²S/∂R_iA∂R_iB`, with
`R_iA`/`R_iB` in bohr (see [Hessians](@ref) for units). Returns a dense
`nbas × nbas × 3 × 3` array. For repeated calls, see `∇2overlap!`.
"""
∇2overlap(BS::BasisSet, iA, iB) = ∇21e(BS, "overlap", iA, iB)
∇2overlap!(out, BS::BasisSet, iA, iB) = ∇21e!(out, BS, "overlap", iA, iB)

"""
    ∇2kinetic(BS::BasisSet, iA, iB) -> Array{Float64,4}

Second derivative (Hessian) of the AO kinetic energy matrix `T` w.r.t.
atoms `iA`,`iB`'s three Cartesian coordinates each, `∂²T/∂R_iA∂R_iB`, with
`R_iA`/`R_iB` in bohr (see [Hessians](@ref) for units). Returns a dense
`nbas × nbas × 3 × 3` array. For repeated calls, see `∇2kinetic!`.
"""
∇2kinetic(BS::BasisSet, iA, iB) = ∇21e(BS, "kinetic", iA, iB)
∇2kinetic!(out, BS::BasisSet, iA, iB) = ∇21e!(out, BS, "kinetic", iA, iB)

"""
    ∇2ERI_2e2c(auxbset::BasisSet, iA, iB)

Second derivative (Hessian, atoms `iA`,`iB`, with `R_iA`/`R_iB` in bohr --
see [Hessians](@ref) for units) of the 2-center two-electron auxiliary
metric `(P|Q)` (density fitting's `J_PQ`). Output `(naux,naux,3,3)`.
Reuses `∇21e!`'s same-shell/cross-shell structure via the `"metric"`
compute-case -- see this file's header comment.
"""
∇2ERI_2e2c(BS::BasisSet, iA, iB) = ∇21e(BS, "metric", iA, iB)
∇2ERI_2e2c!(out, BS::BasisSet, iA, iB) = ∇21e!(out, BS, "metric", iA, iB)

# --- Shell-pair level -------------------------------------------------
#
# The level a Hessian caller usually wants: for the 2-electron integrals a
# dense Hessian is a (nbas^4,3,3) tensor and simply not materializable, so
# integral-direct code has to work block by block anyway. These give the same
# per-block entry point for the 1-electron integrals, so a CPHF or
# integral-direct Hessian loop has one uniform shape across the whole family.

"""
    ∇2overlap!(out, BS::BasisSet{LCint}, iA, iB, i::Int, j::Int; scratch=nothing)
    ∇2kinetic!(out, BS::BasisSet{LCint}, iA, iB, i::Int, j::Int; scratch=nothing)
    ∇2ERI_2e2c!(out, BS::BasisSet{LCint}, iA, iB, i::Int, j::Int; scratch=nothing)

Shell-pair Hessian block: `∂²X_ij/∂R_iA∂R_iB` for shells `i,j` of `BS`
(shell indices, not AO indices). `out` must be `(Ni,Nj,3,3)` and is
overwritten, so a reused buffer is safe.

Returns the free zero -- no libcint call -- whenever neither shell sits on
`iA`, or neither sits on `iB`: the block cannot depend on both coordinates
then.

Pass `scratch` (a vector of at least `9*Ni*Nj` elements, or `9*Nmax^2` to
cover any pair) to make the call allocation-free in a loop.

See [`∇2overlap`](@ref) for the dense whole-matrix form, and
[`∇2overlap_μμ!`](@ref)/[`∇2overlap_μν!`](@ref) for the bare primitives
underneath.
"""
function ∇2overlap!(out, BS::BasisSet{LCint}, iA, iB, i::Int, j::Int; scratch=nothing)
    _∇21e_pair_checked!(∇2overlap_μμ!, ∇2overlap_μν!, out, BS, iA, iB, i, j, scratch)
end

@doc (@doc ∇2overlap!(::Any, ::BasisSet{LCint}, ::Any, ::Any, ::Int, ::Int))
function ∇2kinetic!(out, BS::BasisSet{LCint}, iA, iB, i::Int, j::Int; scratch=nothing)
    _∇21e_pair_checked!(∇2kinetic_μμ!, ∇2kinetic_μν!, out, BS, iA, iB, i, j, scratch)
end

@doc (@doc ∇2overlap!(::Any, ::BasisSet{LCint}, ::Any, ::Any, ::Int, ::Int))
function ∇2ERI_2e2c!(out, BS::BasisSet{LCint}, iA, iB, i::Int, j::Int; scratch=nothing)
    _∇21e_pair_checked!(∇2ERI_2e2c_μμ!, ∇2ERI_2e2c_μν!, out, BS, iA, iB, i, j, scratch)
end

function _∇21e_pair_checked!(μμ!, μν!, out, BS, iA, iB, i, j, scratch)
    Ni = num_basis(BS[i]); Nj = num_basis(BS[j])
    if size(out) != (Ni, Nj, 3, 3)
        throw(DimensionMismatch("Size of the output array needs to be ($Ni, $Nj, 3, 3)"))
    end
    fill!(out, 0.0)

    X_i, X_j = GaussianBasis.on_atom_flags(BS, iA, i, j)
    Y_i, Y_j = GaussianBasis.on_atom_flags(BS, iB, i, j)
    ((X_i || X_j) && (Y_i || Y_j)) || return out   # free zero

    buf = scratch === nothing ? Vector{Cdouble}(undef, 9*Ni*Nj) : scratch
    length(buf) >= 9*Ni*Nj ||
        throw(DimensionMismatch("scratch must hold at least $(9*Ni*Nj) elements"))
    return _∇21e_pair!(μμ!, μν!, out, BS, i, j, Ni, Nj, X_i, X_j, Y_i, Y_j, buf)
end

# Allocating shell-pair forms.
function ∇2overlap(BS::BasisSet, iA, iB, i::Int, j::Int)
    out = zeros(num_basis(BS[i]), num_basis(BS[j]), 3, 3)
    return ∇2overlap!(out, BS, iA, iB, i, j)
end
function ∇2kinetic(BS::BasisSet, iA, iB, i::Int, j::Int)
    out = zeros(num_basis(BS[i]), num_basis(BS[j]), 3, 3)
    return ∇2kinetic!(out, BS, iA, iB, i, j)
end
function ∇2ERI_2e2c(BS::BasisSet, iA, iB, i::Int, j::Int)
    out = zeros(num_basis(BS[i]), num_basis(BS[j]), 3, 3)
    return ∇2ERI_2e2c!(out, BS, iA, iB, i, j)
end
