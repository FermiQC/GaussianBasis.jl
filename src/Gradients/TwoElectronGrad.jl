"""
    ∇ERI_2e4c(BS::BasisSet, iA) -> Array{Float64,5}

Gradient of the full two-electron four-center integral tensor `(ij|kl)`
w.r.t. atom `iA`'s three Cartesian coordinates, with `R_iA` in bohr (see
[Gradients](@ref) for units). Returns a dense `nbas × nbas × nbas × nbas ×
3` array respecting the same 8-fold permutational symmetry as `ERI_2e4c`.
This is the full, uncompressed tensor -- for large basis sets prefer
`∇sparseERI_2e4c`, which screens and stores only the unique elements. For a
single shell quartet, see [`∇ERI_2e4c(BS,iA,i,j,k,l)`](@ref
∇ERI_2e4c(::BasisSet, ::Int, ::Int, ::Int, ::Int, ::Int)). For repeated
calls, see `∇ERI_2e4c!`.
"""
function ∇ERI_2e4c(BS::BasisSet, iA)
    # Pre allocate output
    out = zeros(BS.nbas, BS.nbas, BS.nbas, BS.nbas, 3)
    return ∇ERI_2e4c!(out, BS, iA)
end

"""
    ∇ERI_2e4c!(out, BS::BasisSet, iA)

Mutating counterpart of [`∇ERI_2e4c`](@ref): writes the dense
`nbas × nbas × nbas × nbas × 3` gradient into the caller-supplied `out`
instead of allocating it.

For a single shell quartet see
[`∇ERI_2e4c(BS,iA,i,j,k,l)`](@ref ∇ERI_2e4c(::BasisSet, ::Int, ::Int, ::Int, ::Int, ::Int));
for large basis sets prefer [`∇sparseERI_2e4c`](@ref).
"""
function ∇ERI_2e4c!(out, BS::BasisSet, iA)

    if size(out) != (BS.nbas, BS.nbas, BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (N, N, N, N, 3)."))
    end

    A = BS.atoms[iA]

    # Shell indexes for basis in the atom A
    Ashells = Int[]
    notAshells = Int[]
    for i in 1:BS.nshells
        b = BS.shells[i]
        if b.atom == A
            push!(Ashells, i)
        else
            push!(notAshells, i)
        end
    end

    Nvals = num_basis.(BS.shells)
    ao_offset = cumsum(Nvals) .- Nvals
    Nmax = maximum(Nvals)
    #buf_arrays = [zeros(Cdouble, 3*Nmax^4) for _ = 1:Threads.nthreads()]


    # Find unique (i,j,k,l) combinations given permutational symmetry
    unique_idx = Tuple{Int, Int, Int, Int}[]
    for i = 1:BS.nshells
        for j = i:BS.nshells # i <= j
            for k = 1:BS.nshells
                for l = k:BS.nshells # k <= l
                    if index2(i-1,j-1) < index2(k-1,l-1)
                        continue
                    end
                    push!(unique_idx, (i,j,k,l))
                end
            end
        end
    end

    allocate(body) = body(zeros(Cdouble, 3*Nmax^4))
    workerpool(allocate, unique_idx; chunksize=10) do (i,j,k,l), buf
        x_in_A = map(in(Ashells), (i,j,k,l))

        # If no basis is centered on A, skip
        # If all basis are centered on A, skip
        if !any(x_in_A) || all(x_in_A)
            return
        end

        Ni = Nvals[i]
        Nj = Nvals[j]
        Nk = Nvals[k]
        Nl = Nvals[l]
        Nijkl = Ni*Nj*Nk*Nl

        ioff = ao_offset[i]
        joff = ao_offset[j]
        koff = ao_offset[k]
        loff = ao_offset[l]

        I = (ioff+1):(ioff+Ni)
        J = (joff+1):(joff+Nj)
        K = (koff+1):(koff+Nk)
        L = (loff+1):(loff+Nl)

        # [i'j|kl]
        if x_in_A[1]
            cint2e_ip1_sph!(buf, @SVector([i,j,k,l]), BS.lib)
            ∇q = reshape(view(buf, 1:3*Nijkl), Ni, Nj, Nk, Nl, 3)
            for q in 1:3
                @inbounds for (rl, nl) in enumerate(L), (rk, nk) in enumerate(K), (rj, nj) in enumerate(J), (ri, ni) in enumerate(I)
                    out[ni, nj, nk, nl, q] -= ∇q[ri, rj, rk, rl, q]
                end
            end
        end

        # [ij'|kl]
        if x_in_A[2]
            cint2e_ip1_sph!(buf, @SVector([j,i,k,l]), BS.lib)
            ∇q = reshape(view(buf, 1:3*Nijkl), Nj, Ni, Nk, Nl, 3)
            for q in 1:3
                @inbounds for (rl, nl) in enumerate(L), (rk, nk) in enumerate(K), (rj, nj) in enumerate(J), (ri, ni) in enumerate(I)
                    out[ni, nj, nk, nl, q] -= ∇q[rj, ri, rk, rl, q]
                end
            end
        end

        # [ij|k'l]
        if x_in_A[3]
            cint2e_ip1_sph!(buf, @SVector([k,l,i,j]), BS.lib)
            ∇q = reshape(view(buf, 1:3*Nijkl), Nk, Nl, Ni, Nj, 3)
            for q in 1:3
                @inbounds for (rl, nl) in enumerate(L), (rk, nk) in enumerate(K), (rj, nj) in enumerate(J), (ri, ni) in enumerate(I)
                    out[ni, nj, nk, nl, q] -= ∇q[rk, rl, ri, rj, q]
                end
            end
        end

        # [ij|kl']
        if x_in_A[4]
            cint2e_ip1_sph!(buf, @SVector([l,k,i,j]), BS.lib)
            ∇q = reshape(view(buf, 1:3*Nijkl), Nl, Nk, Ni, Nj, 3)
            for q in 1:3
                @inbounds for (rl, nl) in enumerate(L), (rk, nk) in enumerate(K), (rj, nj) in enumerate(J), (ri, ni) in enumerate(I)
                    out[ni, nj, nk, nl, q] -= ∇q[rl, rk, ri, rj, q]
                end
            end
        end

        for q in 1:3
            if i != j && k != l && index2(i,j) != index2(k,l)
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[nj, ni, nk, nl, q] = out[ni, nj, nk, nl, q]
                    out[ni, nj, nl, nk, q] = out[ni, nj, nk, nl, q]
                    out[nj, ni, nl, nk, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, ni, nj, q] = out[ni, nj, nk, nl, q]
                    out[nl, nk, ni, nj, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, nj, ni, q] = out[ni, nj, nk, nl, q]
                    out[nl, nk, nj, ni, q] = out[ni, nj, nk, nl, q]
                end
            elseif k != l && index2(i,j) != index2(k,l)
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[ni, nj, nl, nk, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, ni, nj, q] = out[ni, nj, nk, nl, q]
                    out[nl, nk, ni, nj, q] = out[ni, nj, nk, nl, q]
                end
            elseif i != j && index2(i,j) != index2(k,l)
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[nj, ni, nk, nl, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, ni, nj, q] = out[ni, nj, nk, nl, q]
                    out[nk, nl, nj, ni, q] = out[ni, nj, nk, nl, q]
                end
            elseif i != j && k != l 
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[nj, ni, nk, nl, q] = out[ni, nj, nk, nl, q]
                    out[ni, nj, nl, nk, q] = out[ni, nj, nk, nl, q]
                    out[nj, ni, nl, nk, q] = out[ni, nj, nk, nl, q]
                end
            elseif index2(i,j) != index2(k,l) 
                @inbounds for nl = L, nk = K, nj = J, ni = I
                    out[nk, nl, ni, nj, q] = out[ni, nj, nk, nl, q]
                end
            end
        end
    end

    return out
end

# --- Shell-quartet-level gradient ---
#
# GaussianBasis-idiomatic single-quartet primitive, mirroring the shell-pair
# primitives added for the 1-electron gradients (see OneElectronGrad.jl's
# header comment for the general rationale) -- motivated here by two things
# `∇ERI_2e4c!`'s own whole-array loop can't offer: (1) a building block for
# genuinely integral-direct gradient/CPHF code (compute one quartet,
# contract it immediately, discard -- unlike the plain energy ERI, which is
# reused across every SCF/CPHF iteration, a derivative quartet is only ever
# needed once to help form Jq/Kq or a Hessian block, so there's no caching
# benefit being given up by not materializing anything), and (2) avoiding
# GaussianBasis.Libcint's raw kernel calls at the call site the way Fermi.jl's
# `TwoElectronHess.jl` currently has to.
#
# Unlike `∇nuclear`, the bare Coulomb operator has no third "operator
# center" to differentiate -- so the free-zero case is exactly the same
# shape as overlap/kinetic's: a quartet's derivative is zero, no libcint
# call needed, whenever ALL FOUR shells share atom-membership status (all on
# atom `iA`, or all off it), by translational invariance. `∇ERI_2e4c!`'s
# own loop already exploits exactly this (`!any(x_in_A) || all(x_in_A)`,
# line 54 above).
#
# For a surviving quartet, this reuses `∇ERI_2e4c!`'s own per-shell
# derivative terms (`[i'j|kl]`, `[ij'|kl]`, `[ij|k'l]`, `[ij|kl']`, one
# `cint2e_ip1_sph!` call and index-permutation per shell that's actually on
# atom `iA`) directly, for an ARBITRARY shell ordering -- unlike the
# whole-array loop, this function does not need `i≤j`, `k≤l`, or any
# canonical-quartet reordering first: that reordering in `∇ERI_2e4c!` is
# purely a performance optimization (compute an 8-fold-symmetric block once,
# propagate it to every symmetric position instead of recomputing), not a
# correctness requirement -- the underlying `[i'j|kl]+[ij'|kl]+[ij|k'l]+
# [ij|kl']` formula is already well-defined for any `(i,j,k,l)` on its own.

"""
    ∇ERI_2e4c_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int, l::Int)

Libcint call for the 4-center ERI gradient block differentiated with respect
to shell `i` (the "μ" shell), for shells `i,j,k,l` of `BS`, already
sign-flipped to the nuclear-coordinate convention -- the direct analogue of
[`∇overlap_μ!`](@ref), and the only place this file touches libcint.

Always differentiates its FIRST shell argument. To differentiate any other
center, pass that shell first and permute the result, exactly as the 1e
gradients do with `∇overlap_μ!(out, BS, j, i)`. `out` comes back in
libcint's raw `(Ni,Nj,Nk,Nl,3)` layout for the shell order **as passed**, so
a permuted call returns a permuted block -- see `∇ERI_2e4c!` for the four
index maps that fold those permutations into the scatter.

Handles the 1-based to 0-based shell index conversion and passes the indices
as an `SVector`, so no caller-owned `shls` buffer is needed and the call is
allocation-free. `out` may be a contiguous view (e.g. `view(buf, 1:3*Nijkl)`).

> No bounds checking, no zero-block skipping, no output-size validation --
> same segfault/heap-corruption risks as [`∇overlap_μ!`](@ref).
"""
function ∇ERI_2e4c_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int, l::Int)
    lib = BS.lib
    cint2e_ip1_sph!(out, @SVector(Cint[i-1, j-1, k-1, l-1]),
                    lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    out .*= -1.0
    return out
end

"""
    ∇ERI_2e4c(BS::BasisSet, iA::Int, i::Int, j::Int, k::Int, l::Int)
    ∇ERI_2e4c(BS::BasisSet, on_A::NTuple{4,Bool}, i::Int, j::Int, k::Int, l::Int)

Two forms. The `iA::Int` form is the convenient, standalone-safe one:
computes `on_A` (which of shells `i,j,k,l` sit on atom `iA`, via `===` --
see `on_atom_flags`) and returns the free zero (no libcint call) when
`i,j,k,l` are ALL on atom `iA` or ALL off it. The `on_A::NTuple{4,Bool}`
form is the unchecked core: it does whatever `on_A` says unconditionally,
no membership computation and no zero-skip, for callers (like Fermi.jl's
gradient loop) that have already screened at a higher level and precomputed
`on_A` once outside a hot loop -- calling this form on an all-on/all-off
`on_A` does NOT raise an error, it just wastefully computes an answer of
exactly zero the long way, so getting that screening right is entirely the
caller's responsibility for this form.
"""
function ∇ERI_2e4c(BS::BasisSet, iA::Int, i::Int, j::Int, k::Int, l::Int)
    on_A = on_atom_flags(BS, iA, i, j, k, l)
    Ni, Nj, Nk, Nl = num_basis(BS.shells[i]), num_basis(BS.shells[j]), num_basis(BS.shells[k]), num_basis(BS.shells[l])
    out = zeros(Ni, Nj, Nk, Nl, 3)
    (!any(on_A) || all(on_A)) && return out
    return ∇ERI_2e4c!(out, BS, on_A, i, j, k, l)
end

function ∇ERI_2e4c!(out, BS::BasisSet, iA::Int, i::Int, j::Int, k::Int, l::Int)
    on_A = on_atom_flags(BS, iA, i, j, k, l)
    if !any(on_A) || all(on_A)
        out .= 0.0
        return out
    end
    return ∇ERI_2e4c!(out, BS, on_A, i, j, k, l)
end

function ∇ERI_2e4c(BS::BasisSet, on_A::NTuple{4,Bool}, i::Int, j::Int, k::Int, l::Int)
    Ni, Nj, Nk, Nl = num_basis(BS.shells[i]), num_basis(BS.shells[j]), num_basis(BS.shells[k]), num_basis(BS.shells[l])
    out = zeros(Ni, Nj, Nk, Nl, 3)
    return ∇ERI_2e4c!(out, BS, on_A, i, j, k, l)
end

function ∇ERI_2e4c!(out, BS::BasisSet, on_A::NTuple{4,Bool}, i::Int, j::Int, k::Int, l::Int)
    Ni, Nj, Nk, Nl = num_basis(BS.shells[i]), num_basis(BS.shells[j]), num_basis(BS.shells[k]), num_basis(BS.shells[l])
    buf = Vector{Cdouble}(undef, 3*Ni*Nj*Nk*Nl)
    return ∇ERI_2e4c!(out, BS, on_A, i, j, k, l, buf)
end

"""
    ∇ERI_2e4c!(out, BS::BasisSet, on_A::NTuple{4,Bool}, i, j, k, l, buf)

Scratch-buffer-accepting core: identical math to the 7-argument form above,
but takes a caller-owned `buf` (sized `>= 3*Nmax^4`, `Nmax` = the largest
`num_basis` over any shell the caller will ever pass) instead of allocating
one fresh every call -- zero-allocation, for callers in a hot per-quartet
loop. Not thread-safe to share: each concurrent caller (e.g. each worker
task) needs its own `buf`.
"""
function ∇ERI_2e4c!(out, BS::BasisSet, on_A::NTuple{4,Bool}, i::Int, j::Int, k::Int, l::Int,
                     buf::Vector{Cdouble})
    # Exists because this sits in the innermost loop of Fermi.jl's
    # integral-direct gradient (called once per (atom, canonical-quartet)
    # visit, millions of times for a real molecule) -- profiling found the
    # allocating form costing ~1.3 KB/call, several GB of GC churn over a
    # full gradient. Note that permutedims! is NOT a fix for that: despite
    # being the in-place form it still allocates ~384 B per call, so a
    # quartet hitting two transposing branches leaked 768 B even with
    # caller-owned buffers. Folding the permutations into the index
    # expressions below is what actually makes this allocation-free.

    Ni, Nj, Nk, Nl = num_basis(BS.shells[i]), num_basis(BS.shells[j]), num_basis(BS.shells[k]), num_basis(BS.shells[l])

    if size(out) != (Ni, Nj, Nk, Nl, 3)
        throw(DimensionMismatch("Size of the output array needs to be ($Ni, $Nj, $Nk, $Nl, 3)."))
    end

    fill!(out, 0.0)

    Nijkl = Ni*Nj*Nk*Nl
    bufv = view(buf, 1:3*Nijkl)

    # Every branch is the same primitive, ∇ERI_2e4c_μ! -- which always
    # differentiates its first shell argument -- called with the shell to be
    # differentiated moved to the front. The resulting block therefore comes
    # back permuted, and each branch folds that permutation into the index
    # expression rather than materializing a transposed copy.
    #
    # Index maps follow permutedims' convention, dest[j...] =
    # src[j[invperm(perm)]...], with the raw block laid out in the shell order
    # that branch passed to libcint.

    # [i'j|kl] -- raw (Ni,Nj,Nk,Nl,3), no permutation
    if on_A[1]
        ∇ERI_2e4c_μ!(bufv, BS, i, j, k, l)
        @inbounds for q = 1:3
            oq = Nijkl*(q-1)
            for d = 1:Nl, c = 1:Nk, b = 1:Nj, a = 1:Ni
                out[a,b,c,d,q] += bufv[oq + a + Ni*(b-1) + Ni*Nj*(c-1) + Ni*Nj*Nk*(d-1)]
            end
        end
    end

    # [ij'|kl] -- raw (Nj,Ni,Nk,Nl,3), perm (2,1,3,4,5)
    if on_A[2]
        ∇ERI_2e4c_μ!(bufv, BS, j, i, k, l)
        @inbounds for q = 1:3
            oq = Nijkl*(q-1)
            for d = 1:Nl, c = 1:Nk, b = 1:Nj, a = 1:Ni
                out[a,b,c,d,q] += bufv[oq + b + Nj*(a-1) + Nj*Ni*(c-1) + Nj*Ni*Nk*(d-1)]
            end
        end
    end

    # [ij|k'l] -- raw (Nk,Nl,Ni,Nj,3), perm (3,4,1,2,5)
    if on_A[3]
        ∇ERI_2e4c_μ!(bufv, BS, k, l, i, j)
        @inbounds for q = 1:3
            oq = Nijkl*(q-1)
            for d = 1:Nl, c = 1:Nk, b = 1:Nj, a = 1:Ni
                out[a,b,c,d,q] += bufv[oq + c + Nk*(d-1) + Nk*Nl*(a-1) + Nk*Nl*Ni*(b-1)]
            end
        end
    end

    # [ij|kl'] -- raw (Nl,Nk,Ni,Nj,3), perm (3,4,2,1,5)
    if on_A[4]
        ∇ERI_2e4c_μ!(bufv, BS, l, k, i, j)
        @inbounds for q = 1:3
            oq = Nijkl*(q-1)
            for d = 1:Nl, c = 1:Nk, b = 1:Nj, a = 1:Ni
                out[a,b,c,d,q] += bufv[oq + d + Nl*(c-1) + Nl*Nk*(a-1) + Nl*Nk*Ni*(b-1)]
            end
        end
    end

    return out
end

"""
    schwarz_bounds(BS::BasisSet)

Per-shell-pair Cauchy-Schwarz screening bound `σ_ij := sqrt(max|(ij|ij)|)`,
used by both `sparseERI_2e4c` and `∇sparseERI_2e4c` to skip shell quartets
whose bound falls below a cutoff. Independent of which atom is being
differentiated -- callers making several `∇sparseERI_2e4c` calls (e.g.
looping over atoms) should compute it once and pass it through via the
`ij_vals`/`σvals` kwargs rather than recomputing it each time.
"""
function schwarz_bounds(BS::BasisSet)
    Nvals = num_basis.(BS.shells)
    Nmax = maximum(Nvals)
    num_ij = Int((BS.nshells^2 - BS.nshells)/2) + BS.nshells
    ij_vals = Array{NTuple{2,Int32}}(undef, num_ij)
    σvals = zeros(Cdouble, num_ij)
    tmp = zeros(Cdouble, Nmax^4)
    for i = 1:BS.nshells
        for j = i:BS.nshells
            idx = index2(i-1,j-1) + 1
            ij_vals[idx] = (i,j)
            ERI_2e4c!(tmp, BS, i, i, j, j)
            σvals[idx] = √maximum(tmp)
        end
    end
    return ij_vals, σvals
end

"""
    ∇sparseERI_2e4c(BS::BasisSet, iA, cutoff=1e-12; ij_vals=nothing, σvals=nothing)

Derivative (w.r.t. atom `iA`'s three Cartesian directions, in bohr -- see
[Gradients](@ref) for units) of the unique (permutation-compressed)
two-electron four-center integrals, Schwarz-screened the same way
`sparseERI_2e4c` screens the energy integrals (see `schwarz_bounds`).
`ij_vals`/`σvals` default to a fresh `schwarz_bounds(BS)` call if not
supplied -- pass a precomputed pair in when calling this repeatedly across
atoms to avoid recomputing it every time.
"""
function ∇sparseERI_2e4c(BS::BasisSet, iA, cutoff = 1e-12; ij_vals = nothing, σvals = nothing)
    # The energy-integral Schwarz bound is a valid screening proxy for the
    # derivative too -- differentiating a Gaussian-product ERI w.r.t. a
    # nuclear center only introduces a bounded polynomial prefactor via the
    # chain rule, so the exponential shell-pair falloff with distance that
    # makes the bound work for the energy integral is unchanged for its
    # derivative. Quartets entirely on atom iA or entirely off it are also
    # skipped (exactly zero by translational invariance, not a cutoff).
    #
    # Results accumulate into a growable buffer (sizehint!ed from a cheap
    # upper-bound estimate, no extra integral evaluations needed for it)
    # rather than a fixed array pre-sized at the full theoretical unique-
    # element count -- peak memory scales with surviving elements, not
    # nbas^4/8. Deliberately NOT threaded, unlike sparseERI_2e4c (called
    # once per SCF run, amortizing threading overhead over every Fock build
    # that follows): this is typically called once per atom per gradient/
    # Hessian evaluation, too little work per call to amortize a thread
    # pool's synchronization cost -- measured ~2x slower threaded than
    # serial for a small molecule.
    A = BS.atoms[iA]

    # Shell indexes for basis in the atom A
    in_A = falses(BS.nshells)
    for i in 1:BS.nshells
        BS.shells[i].atom == A && (in_A[i] = true)
    end

    # Pre compute a list of number of basis for each shell (2l +1)
    Nvals = num_basis.(BS.shells)
    ao_offset = cumsum(Nvals) .- Nvals
    Nmax = maximum(Nvals)

    # Unique shell pairs with i ≤ j
    num_ij = Int((BS.nshells^2 - BS.nshells)/2) + BS.nshells

    if ij_vals === nothing || σvals === nothing
        ij_vals, σvals = schwarz_bounds(BS)
    end

    # Candidate shell quartets (ij,kl), ij ≤ kl: Schwarz-screened AND
    # touching atom A in some but not all of the four slots. A plain nested
    # loop rather than a filtered/flattened generator -- eltype() on the
    # latter gives up and infers Any (see sparseERI_2e4c's comment on the
    # same issue), which would box every (i,j,k,l) tuple in the hot loop
    # below; measured as the dominant cost here, more than any of the actual
    # per-quartet integral work.
    #
    # Upper bound on the number of surviving AO elements (Nvals products
    # only, no ERI calls), used to sizehint! the buffers below.
    size_ub = 0
    @inbounds for ij in 1:num_ij
        i, j = ij_vals[ij]
        for kl in ij:num_ij
            σvals[ij]*σvals[kl] <= cutoff && continue
            k, l = ij_vals[kl]
            nA = in_A[i] + in_A[j] + in_A[k] + in_A[l]
            (nA == 0 || nA == 4) && continue
            size_ub += Nvals[i]*Nvals[j]*Nvals[k]*Nvals[l]
        end
    end

    indexes = sizehint!(Vector{NTuple{4,Int16}}(), size_ub)
    ∇x = sizehint!(Cdouble[], size_ub)
    ∇y = sizehint!(Cdouble[], size_ub)
    ∇z = sizehint!(Cdouble[], size_ub)

    buf = Vector{Cdouble}(undef, 3*Nmax^4)
    # Block staging area for ∇ERI_2e4c!; contiguous so the emit loop can index
    # it linearly. `tmp`/`shls` are ignored by that core (see its docstring)
    # but are still positional, so pass empty vectors.
    blkbuf = Vector{Cdouble}(undef, 3*Nmax^4)

    # i,j,k,l => Shell indexes starting at one
    # I, J, K, L => AO indexes starting at one
    @inbounds for ij in 1:num_ij
        i, j = ij_vals[ij]
        for kl in ij:num_ij
            σvals[ij]*σvals[kl] <= cutoff && continue
            k, l = ij_vals[kl]
            nA = in_A[i] + in_A[j] + in_A[k] + in_A[l]
            (nA == 0 || nA == 4) && continue
        begin
            Ni, Nj, Nk, Nl = Nvals[i], Nvals[j], Nvals[k], Nvals[l]
            Nij = Ni*Nj
            Nijk = Nij*Nk
            Nijkl = Nijk*Nl
            ioff, joff, koff, loff = ao_offset[i], ao_offset[j], ao_offset[k], ao_offset[l]

            # One call to the shared shell-quartet core instead of
            # re-implementing its four branches here. It writes a
            # (Ni,Nj,Nk,Nl,3) block into `blkbuf`; because that view is
            # contiguous, component q of the block occupies the linear range
            # Nijkl*(q-1)+1 : q*Nijkl, so the emit loop below can index
            # `blkbuf` directly with the same `is + bjkl` offsets the old
            # bufx/bufy/bufz used.
            #
            # The core already folds libcint's sign flip into ∇ERI_2e4c_μ!,
            # where the old code here accumulated raw kernel output and
            # negated at push! time -- hence the pushes below no longer carry
            # a minus sign.
            blk = reshape(view(blkbuf, 1:3*Nijkl), Ni, Nj, Nk, Nl, 3)
            ∇ERI_2e4c!(blk, BS, (in_A[i], in_A[j], in_A[k], in_A[l]),
                       Int(i), Int(j), Int(k), Int(l), buf)

            ### This block aims to retrieve unique elements within buf and map them to AO indexes
            # is, js, ks, ls are indexes within the shell e.g. for a p shell is = (1, 2, 3)
            # bl, bkl, bjkl are used to map the (i,j,k,l) index into a one-dimensional index for buf
            # That is, get the correct integrals for the AO quartet.
            #
            # Only when the same shell-pair is reused for bra and ket (i==k && j==l)
            # does this block become self-symmetric (is,js and ks,ls range over the
            # exact same configuration space), so every (I,J,K,L) generated below has
            # a mirror (K,L,I,J) also generated in this same shell quartet -- guarded
            # by IJ>KL&&continue to push each AO quartet once, not twice.
            self_paired = i == k && j == l
            for ls = 1:Nl
                L = loff + ls
                bl = Nijk*(ls-1)
                for ks = 1:Nk
                    K = koff + ks
                    L < K && break

                    # L ≥ K
                    # index2 for K,L
                    KL = ((L * (L - 1)) ÷ 2) + K - 1

                    bkl = Nij*(ks-1) + bl
                    for js = 1:Nj
                        J = joff + js
                        bjkl = Ni*(js-1) + bkl
                        for is = 1:Ni
                            I = ioff + is
                            J < I && break

                            # index2 for I,J
                            IJ = ((J * (J - 1)) ÷ 2) + I - 1

                            self_paired && IJ > KL && continue

                            n = is + bjkl
                            push!(∇x, blkbuf[n])
                            push!(∇y, blkbuf[Nijkl + n])
                            push!(∇z, blkbuf[2*Nijkl + n])
                            push!(indexes, KL ≥ IJ ? (I, J, K, L) : (K, L, I, J))
                        end
                    end
                end
            end
        end
        end
    end

    return indexes, ∇x, ∇y, ∇z
end

"""
    ∇ERI_2e3c(BS1::BasisSet, BS2::BasisSet, iA) -> Array{Float64,4}

Gradient of the full two-electron three-center integral tensor `(μν|P)`
(`BS1`=regular basis, `BS2`=auxiliary/fitting basis) w.r.t. atom `iA`'s
three Cartesian coordinates, with `R_iA` in bohr (see [Gradients](@ref) for
units). `iA` indexes into `BS1.atoms`. Returns a dense `BS1.nbas ×
BS1.nbas × BS2.nbas × 3` array, symmetric under `μ↔ν` swap. For repeated
calls, see `∇ERI_2e3c!`.
"""
function ∇ERI_2e3c(BS1::BasisSet, BS2::BasisSet, iA)
    # Pre allocate output
    out = zeros(BS1.nbas, BS1.nbas, BS2.nbas, 3)
    return ∇ERI_2e3c!(out, BS1, BS2, iA)
end

"""
    ∇ERI_2e3c_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int)
    ∇ERI_2e3c_P!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int)

Libcint calls for the 3-center `(μν|P)` gradient block, sign-flipped to the
nuclear-coordinate convention -- the analogues of [`∇overlap_μ!`](@ref) for
this integral.

`BS` is the **merged** basis (regular shells followed by auxiliary ones, as
`ERI_2e3c!(out, BS, i, j, k)` also expects), so an auxiliary shell `k` is
addressed as `k + BS1.nshells`.

Two primitives are needed rather than one, because the three centers do not
live in interchangeable slots: `μ` and `ν` are the two shells of the bra,
while `P` sits alone in the ket. `∇ERI_2e3c_μ!` (libcint's `ip1`)
differentiates the FIRST shell argument, so `ν` is reached by swapping the
first two arguments -- the same trick used throughout this package.
`∇ERI_2e3c_P!` (libcint's `ip2`) differentiates the third center; no
argument permutation can express it in terms of `ip1`.

`out` comes back in libcint's raw layout for the shell order as passed, and
may be a contiguous view.

> No bounds checking or output-size validation -- same risks as
> [`∇overlap_μ!`](@ref).
"""
function ∇ERI_2e3c_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int)
    lib = BS.lib
    cint3c2e_ip1_sph!(out, @SVector(Cint[i-1, j-1, k-1]),
                      lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    out .*= -1.0
    return out
end

function ∇ERI_2e3c_P!(out, BS::BasisSet{LCint}, i::Int, j::Int, k::Int)
    lib = BS.lib
    cint3c2e_ip2_sph!(out, @SVector(Cint[i-1, j-1, k-1]),
                      lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    out .*= -1.0
    return out
end

"""
    ∇ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet, iA; Bmerged=nothing)

Mutating counterpart of [`∇ERI_2e3c`](@ref): writes the dense
`BS1.nbas × BS1.nbas × BS2.nbas × 3` gradient into `out` instead of
allocating it. `out` is zeroed first, so a reused buffer is safe.

`Bmerged` depends only on `BS1`/`BS2`, never on `iA`. Callers looping over
atoms should build it once and pass it in -- it is a small fraction of the
runtime but over half of this routine's allocation.
"""
function ∇ERI_2e3c!(out, BS1::BasisSet, BS2::BasisSet, iA; Bmerged::Union{Nothing,BasisSet}=nothing)

    # Bmerged depends only on BS1/BS2, never on iA -- callers looping over
    # atoms (e.g. Fermi.jl's DF-gradient atom loops) can build it once and
    # pass it in via this keyword instead of paying for a fresh BasisSet
    # construction on every call.
    if Bmerged === nothing
        Bmerged = merge_basis(BS1, BS2)
    end

    if size(out) != (BS1.nbas, BS1.nbas, BS2.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (N1, N1, N2, 3)."))
    end

    A = BS1.atoms[iA]

    # Per-shell membership, looked up in O(1) below. Previously this was a
    # pair of index lists searched with `in`, and the per-triple flags were
    # built as a Vector literal -- one heap allocation and two linear scans
    # for every shell triple.
    onA1 = [BS1.shells[i].atom == A for i in 1:BS1.nshells]
    onA2 = [BS2.shells[i].atom == A for i in 1:BS2.nshells]

    Nvals1 = num_basis.(BS1.shells)
    ao_offset1 = cumsum(Nvals1) .- Nvals1
    Nmax1 = maximum(Nvals1)

    Nvals2 = num_basis.(BS2.shells)
    ao_offset2 = cumsum(Nvals2) .- Nvals2
    Nmax2 = maximum(Nvals2)

    buf = Vector{Cdouble}(undef, 3*Nmax1^2*Nmax2)

    # Blocks where all three shells share atom-membership status are zero and
    # are skipped below, so `out` must start clean for a reused buffer.
    fill!(out, 0.0)

    for i = 1:BS1.nshells
        for j = i:BS1.nshells # i <= j
            for k = 1:BS2.nshells

                x_in_A = (onA1[i], onA1[j], onA2[k])

                # If no basis is centered on A, skip
                # If all basis are centered on A, skip
                if !any(x_in_A) || all(x_in_A)
                    continue
                end

                Ni = Nvals1[i]
                Nj = Nvals1[j]
                Nk = Nvals2[k]
                Nijk = Ni*Nj*Nk
                bufv = view(buf, 1:3*Nijk)

                ioff = ao_offset1[i]
                joff = ao_offset1[j]
                koff = ao_offset2[k]

                I = (ioff+1):(ioff+Ni)
                J = (joff+1):(joff+Nj)
                K = (koff+1):(koff+Nk)

                # Each branch is a primitive call plus a scalar scatter that
                # folds any index permutation into the write, so there are no
                # `buf[r]` copies, no permutedims temporaries and no
                # range-indexed broadcasts.
                kk = k + BS1.nshells

	        # [i'j|k] -- raw (Ni,Nj,Nk,3), no permutation
                if x_in_A[1]
                    ∇ERI_2e3c_μ!(bufv, Bmerged, i, j, kk)
                    @inbounds for q = 1:3
                        oq = Nijk*(q-1)
                        for c = 1:Nk, b = 1:Nj, a = 1:Ni
                            out[ioff+a, joff+b, koff+c, q] += bufv[oq + a + Ni*(b-1) + Ni*Nj*(c-1)]
                        end
                    end
                end

                # [ij'|k] -- swap the bra shells; raw (Nj,Ni,Nk,3)
                if x_in_A[2]
                    ∇ERI_2e3c_μ!(bufv, Bmerged, j, i, kk)
                    @inbounds for q = 1:3
                        oq = Nijk*(q-1)
                        for c = 1:Nk, b = 1:Nj, a = 1:Ni
                            out[ioff+a, joff+b, koff+c, q] += bufv[oq + b + Nj*(a-1) + Nj*Ni*(c-1)]
                        end
                    end
                end

                # [ij|k'] -- the ket center, its own kernel; raw (Ni,Nj,Nk,3)
                if x_in_A[3]
                    ∇ERI_2e3c_P!(bufv, Bmerged, i, j, kk)
                    @inbounds for q = 1:3
                        oq = Nijk*(q-1)
                        for c = 1:Nk, b = 1:Nj, a = 1:Ni
                            out[ioff+a, joff+b, koff+c, q] += bufv[oq + a + Ni*(b-1) + Ni*Nj*(c-1)]
                        end
                    end
                end

                # (μν|P) is symmetric under μ<->ν, so mirror the block.
                if i != j
                    @inbounds for q = 1:3, c = 1:Nk, b = 1:Nj, a = 1:Ni
                        out[joff+b, ioff+a, koff+c, q] = out[ioff+a, joff+b, koff+c, q]
                    end
                end
            end
        end
    end

    return out
end

"""
    ∇ERI_2e2c(BS::BasisSet, iA) -> Array{Float64,3}

Gradient of the full two-electron two-center integral matrix `(P|Q)` (the
density-fitting Coulomb metric `J_PQ`) w.r.t. atom `iA`'s three Cartesian
coordinates, with `R_iA` in bohr (see [Gradients](@ref) for units). Returns
a dense `nbas × nbas × 3` array, symmetric under `P↔Q` swap. For repeated
calls, see `∇ERI_2e2c!`.
"""
function ∇ERI_2e2c(BS::BasisSet, iA)
    # Pre allocate output
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇ERI_2e2c!(out, BS, iA)
end

"""
    ∇ERI_2e2c_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int)

Libcint call for the 2-center `(P|Q)` gradient block differentiated with
respect to shell `i`, sign-flipped to the nuclear-coordinate convention --
the analogue of [`∇overlap_μ!`](@ref) for the DF Coulomb metric.

Differentiates its FIRST shell argument, so `Q` is reached by swapping the
two arguments, which returns a `(Nj,Ni,3)` block to be transposed in. `out`
may be a contiguous view.

> No bounds checking or output-size validation -- same risks as
> [`∇overlap_μ!`](@ref).
"""
function ∇ERI_2e2c_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    lib = BS.lib
    cint2c2e_ip1_sph!(out, @SVector(Cint[i-1, j-1]),
                      lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
    out .*= -1.0
    return out
end

"""
    ∇ERI_2e2c!(out, BS::BasisSet, iA)

Mutating counterpart of [`∇ERI_2e2c`](@ref): writes the dense
`nbas × nbas × 3` gradient of the DF Coulomb metric `(P|Q)` into `out`
instead of allocating it. `out` is zeroed first, so a reused buffer is safe.
"""
function ∇ERI_2e2c!(out, BS::BasisSet, iA)

    if size(out) != (BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (N, N, 3)."))
    end

    A = BS.atoms[iA]

    # Shell indexes for basis in the atom A
    Ashells = Int[]
    notAshells = Int[]
    for i in 1:BS.nshells
        b = BS.shells[i]
        if b.atom == A
            push!(Ashells, i)
        else
            push!(notAshells, i)
        end
    end

    Nvals = num_basis.(BS.shells)
    ao_offset = cumsum(Nvals) .- Nvals
    Nmax = maximum(Nvals)

    buf = Vector{Cdouble}(undef, 3*Nmax^2)

    # Blocks with both shells on A (or both off it) are zero and are skipped
    # below, so `out` must start clean for a reused buffer.
    fill!(out, 0.0)

    for i = 1:BS.nshells
        for j = i:BS.nshells # i <= j

            x_in_A = on_atom_flags(BS, iA, i, j)

            # If no basis is centered on A, skip
            # If all basis are centered on A, skip
            if !any(x_in_A) || all(x_in_A)
                continue
            end

            Ni = Nvals[i]
            Nj = Nvals[j]
            Nij = Ni*Nj

            ioff = ao_offset[i]
            joff = ao_offset[j]

            I = (ioff+1):(ioff+Ni)
            J = (joff+1):(joff+Nj)

            bufv = view(buf, 1:3*Nij)

            # [i'|j] -- raw (Ni,Nj,3), no permutation
            if x_in_A[1]
                ∇ERI_2e2c_μ!(bufv, BS, i, j)
                @inbounds for q = 1:3
                    oq = Nij*(q-1)
                    for b = 1:Nj, a = 1:Ni
                        out[ioff+a, joff+b, q] += bufv[oq + a + Ni*(b-1)]
                    end
                end
            end

            # [i|j'] -- swap the arguments; raw (Nj,Ni,3), transposed in
            if x_in_A[2]
                ∇ERI_2e2c_μ!(bufv, BS, j, i)
                @inbounds for q = 1:3
                    oq = Nij*(q-1)
                    for b = 1:Nj, a = 1:Ni
                        out[ioff+a, joff+b, q] += bufv[oq + b + Nj*(a-1)]
                    end
                end
            end

            # (P|Q) is symmetric, so mirror the block.
            if i != j
                @inbounds for q = 1:3, b = 1:Nj, a = 1:Ni
                    out[joff+b, ioff+a, q] = out[ioff+a, joff+b, q]
                end
            end
        end
    end

    return out
end
