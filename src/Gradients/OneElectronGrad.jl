###########################################################
###########################################################
#                        Overlap                          #
###########################################################
###########################################################

# GaussianBasis-idiomatic shell-pair primitive (hides libcint's shell
# indexing/buffer-layout plumbing the same way `overlap!(out, BS, i, j)`
# does for the plain integrals), derived by symbolically tracing the
# existing, already-validated whole-array gradient loop (matching its
# transpose-mirroring) for an arbitrary shell pair (i,j) in either order,
# rather than re-deriving the underlying integral calculus from scratch --
# lower risk of introducing a NEW error, though this is still checked
# against the whole-array form in the test suite before being trusted.
"""
    ∇overlap_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int)

Libcint call for the overlap gradient block differentiated with respect to
shell `i` (the "μ" AO), for shells `i,j` of `BS`, already sign-flipped to
the nuclear-coordinate convention: libcint's `cint1e_ipovlp_sph!` computes
`∂χ/∂r` (w.r.t. the electron coordinate), but a GTO shell only depends on
`r` and its center `R` through `r-R`, so `∂χ/∂R = -∂χ/∂r` -- this writes
`-cint1e_ipovlp_sph!(...)`, i.e. `∂S/∂R_i` directly, matching every other
`∇`-prefixed function in this module.

> No bounds checking, no zero-block skipping, no output-size validation.
> `i,j` must be valid shell indices (`1:BS.nshells`) and `out` must be
> exactly `(Ni,Nj,3)` -- violating either is undefined behavior at the C
> level, not a catchable Julia error: an out-of-range shell index can
> segfault the whole process, and an undersized `out` can silently write
> past its end into unrelated memory. Prefer `∇overlap!`/`∇overlap` unless
> you're writing your own integral-direct loop and have already validated
> `i,j` and `out` yourself.
"""
function ∇overlap_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    cint1e_ipovlp_sph!(out, @SVector([i, j]), BS.lib)
    out .*= -1.0
    return out
end

"""
    ∇overlap!(out, BS::BasisSet{LCint}, A::Int, i::Int, j::Int; scratch=nothing)
    ∇overlap!(out, BS::BasisSet, A)

Mutating counterpart of [`∇overlap`](@ref): writes into the caller-supplied
`out` instead of allocating. The shell-pair form is the primitive the
full-tensor form builds on -- a direct libcint call, exactly zero (returned
without any libcint call) whenever `i,j` are BOTH on atom `A` or BOTH off
it: `S_ij` doesn't depend on `R_A` at all if neither shell sits there, and
translating `i,j` together with atom `A` (both belong to it) leaves their
relative separation -- and so `S_ij` -- unchanged.

Only implemented for the `LCint` backend -- there is no `ACSint` fallback
for gradients.

# Methods

  - `∇overlap!(out, BS, A, i, j)`: `out` must be `(Ni,Nj,3)`, the block
    for shells `i,j` of `BS` (shell indices, not AO indices).
  - `∇overlap!(out, BS, A)`: `out` must be a dense `nbas × nbas × 3`
    array.

The shell-pair form needs a `3*Ni*Nj` scratch vector whenever the
derivative falls on shell `j` (see below). It allocates one per call unless
you pass `scratch`; in a loop over shell pairs, hand it a buffer sized from
`3*Nmax^2` to make the call allocation-free.
"""
function ∇overlap!(out, BS::BasisSet{LCint}, A::Int, i::Int, j::Int; scratch=nothing)

    i_on_A, j_on_A = GaussianBasis.on_atom_flags(BS, A, i, j)

    if i_on_A == j_on_A
        fill!(out, 0.0)
        return out
    end

    Ni = num_basis(BS[i])
    Nj = num_basis(BS[j])

    if size(out) != (Ni, Nj, 3)
        throw(DimensionMismatch("Size of the output array needs to be ($Ni, $Nj, 3)"))
    end

    if i_on_A
        # Derivative lands on shell i (the one on atom A) directly -- raw
        # libcint output for a [i,j] call is already (Ni,Nj,3) in memory,
        # so write straight into `out`, no intermediate buffer needed.
        ∇overlap_μ!(out, BS, i, j)
    else
        # Derivative lands on shell j. `∇overlap_μ!` always differentiates its
        # FIRST shell argument, so swap them -- that hands back a (Nj,Ni,3)
        # block, transposed here into `out`'s (Ni,Nj,3) AO order. S is
        # symmetric, so ∂S_ij/∂R_A = (∂S_ji/∂R_A)ᵀ over the two AO axes.
        buf = scratch === nothing ? Vector{Cdouble}(undef, 3*Ni*Nj) : scratch
        length(buf) >= 3*Ni*Nj ||
            throw(DimensionMismatch("scratch must hold at least $(3*Ni*Nj) elements"))
        ∇overlap_μ!(buf, BS, j, i)
        @inbounds for k = 1:3, js = 1:Nj, is = 1:Ni
            out[is, js, k] = buf[js + Nj*(is-1) + Nj*Ni*(k-1)]
        end
    end
    return out
end

function ∇overlap(BS::BasisSet, A::Int, i::Int, j::Int)
    Ni = num_basis(BS[i])
    Nj = num_basis(BS[j])
    out = zeros(Ni, Nj, 3)
    return ∇overlap!(out, BS, A, i, j)
end

"""
    ∇overlap(BS::BasisSet, A) -> Array{Float64,3}
    ∇overlap(BS::BasisSet, A::Int, i::Int, j::Int) -> Array{Float64,3}

Gradient of the AO overlap matrix `S` w.r.t. atom `A`'s three Cartesian
coordinates, `∂S/∂R_A`, with `R_A` in bohr (see [Gradients](@ref) for
units).

# Methods

  - `∇overlap(BS, A)`: full, dense `nbas × nbas × 3` array.
  - `∇overlap(BS, A, i, j)`: just the `(Ni,Nj,3)` block for shells `i,j`
    of `BS` (shell indices, not AO indices).

For repeated calls, see `∇overlap!`, which writes into a preallocated
array instead of allocating.
"""
function ∇overlap(BS::BasisSet, A)
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇overlap!(out, BS, A)
end

∇overlap!(out, BS::BasisSet, A) = ∇1e!(∇overlap_μ!, out, BS, A)


###########################################################
###########################################################
#                        Kinetic                           #
###########################################################
###########################################################

"""
    ∇kinetic_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int)

Same as [`∇overlap_μ!`](@ref) (same sign-flip reasoning, same lack of
bounds checking), for the kinetic energy operator instead of overlap.
"""
function ∇kinetic_μ!(out, BS::BasisSet{LCint}, i::Int, j::Int)
    cint1e_ipkin_sph!(out, @SVector([i, j]), BS.lib)
    out .*= -1.0
    return out
end

"""
    ∇kinetic!(out, BS::BasisSet{LCint}, A::Int, i::Int, j::Int; scratch=nothing)
    ∇kinetic!(out, BS::BasisSet, A)

Mutating counterpart of [`∇kinetic`](@ref): writes into the caller-supplied
`out` instead of allocating. The shell-pair form is the primitive the
full-tensor form builds on -- see `∇overlap!` for the identical structure
and zero-block cases, `T` in place of `S`.

Only implemented for the `LCint` backend -- there is no `ACSint` fallback
for gradients.

# Methods

  - `∇kinetic!(out, BS, A, i, j)`: `out` must be `(Ni,Nj,3)`, the block
    for shells `i,j` of `BS` (shell indices, not AO indices).
  - `∇kinetic!(out, BS, A)`: `out` must be a dense `nbas × nbas × 3`
    array.

As with [`∇overlap!`](@ref), the shell-pair form takes an optional
`scratch` vector (`>= 3*Ni*Nj` elements) to stay allocation-free in a loop.
"""
function ∇kinetic!(out, BS::BasisSet{LCint}, A::Int, i::Int, j::Int; scratch=nothing)

    i_on_A, j_on_A = GaussianBasis.on_atom_flags(BS, A, i, j)

    if i_on_A == j_on_A
        fill!(out, 0.0)
        return out
    end

    Ni = num_basis(BS[i])
    Nj = num_basis(BS[j])

    if size(out) != (Ni, Nj, 3)
        throw(DimensionMismatch("Size of the output array needs to be ($Ni, $Nj, 3)"))
    end

    if i_on_A
        # Derivative lands on shell i (the one on atom A) directly -- raw
        # libcint output for a [i,j] call is already (Ni,Nj,3) in memory,
        # so write straight into `out`, no intermediate buffer needed.
        ∇kinetic_μ!(out, BS, i, j)
    else
        # Derivative lands on shell j -- swap the arguments so the
        # differentiated shell is first, then transpose the (Nj,Ni,3) result
        # into `out`. See `∇overlap!` for the full reasoning.
        buf = scratch === nothing ? Vector{Cdouble}(undef, 3*Ni*Nj) : scratch
        length(buf) >= 3*Ni*Nj ||
            throw(DimensionMismatch("scratch must hold at least $(3*Ni*Nj) elements"))
        ∇kinetic_μ!(buf, BS, j, i)
        @inbounds for k = 1:3, js = 1:Nj, is = 1:Ni
            out[is, js, k] = buf[js + Nj*(is-1) + Nj*Ni*(k-1)]
        end
    end
    return out
end

function ∇kinetic(BS::BasisSet, A::Int, i::Int, j::Int)
    Ni = num_basis(BS[i])
    Nj = num_basis(BS[j])
    out = zeros(Ni, Nj, 3)
    return ∇kinetic!(out, BS, A, i, j)
end

"""
    ∇kinetic(BS::BasisSet, A) -> Array{Float64,3}
    ∇kinetic(BS::BasisSet, A::Int, i::Int, j::Int) -> Array{Float64,3}

Gradient of the AO kinetic energy matrix `T` w.r.t. atom `A`'s three
Cartesian coordinates, `∂T/∂R_A`, with `R_A` in bohr (see
[Gradients](@ref) for units).

# Methods

  - `∇kinetic(BS, A)`: full, dense `nbas × nbas × 3` array.
  - `∇kinetic(BS, A, i, j)`: just the `(Ni,Nj,3)` block for shells `i,j`
    of `BS` (shell indices, not AO indices).

For repeated calls, see `∇kinetic!`, which writes into a preallocated
array instead of allocating.
"""
function ∇kinetic(BS::BasisSet, A)
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇kinetic!(out, BS, A)
end

∇kinetic!(out, BS::BasisSet, A) = ∇1e!(∇kinetic_μ!, out, BS, A)


###########################################################
###########################################################
#                        Nuclear                           #
###########################################################
###########################################################

"""
    ∇nuclear_μ!(out, BS::BasisSet{LCint}, charge_atm, i::Int, j::Int)

Same as [`∇overlap_μ!`](@ref) (same sign-flip reasoning, same lack of
bounds checking), for a nuclear-attraction-type potential differentiated
with respect to shell `i` (the "μ" AO) instead -- except here the operator
itself is a potential that can move with a nucleus, not a fixed
multiplicative operator like `1` (overlap) or `-∇²/2` (kinetic).
`charge_atm` fixes WHICH potential: it's a libcint atom-description array
(same layout as `BS.lib.atm`, 6 `Cint` slots per atom) with a chosen subset
of nuclear charges zeroed out, so the raw libcint call sums `Zc/|r-Rc|`
only over the nuclei left un-zeroed. Like `∇overlap_μ!`/`∇kinetic_μ!`,
libcint's `cint1e_ipnuc_sph!` computes `∂/∂r` (w.r.t. the electron
coordinate), so this writes `-cint1e_ipnuc_sph!(...)`, i.e. `∂V^{charge}/∂R_i`
directly for the potential `charge_atm` encodes (see [Nuclear](@ref) for the
derivation).

> No bounds checking, no zero-block skipping, no output-size validation,
> and `charge_atm` isn't validated either -- same segfault/heap-corruption
> risks as [`∇overlap_μ!`](@ref) apply here too.
"""
function ∇nuclear_μ!(out, BS::BasisSet{LCint}, charge_atm, i::Int, j::Int)
    cint1e_ipnuc_sph!(out, @SVector(Cint[i-1, j-1]), charge_atm, BS.lib.natm, BS.lib.bas, BS.lib.nbas, BS.lib.env)
    out .*= -1.0
    return out
end

"""
    nuclear_charge_sets(BS::BasisSet{LCint}, A::Int) -> (only_A, no_A)

The two charge-fudged copies of `BS.lib.atm` that the nuclear gradient needs:
`only_A` keeps atom `A`'s nuclear charge and zeroes every other, `no_A` does
the reverse. Passing either to [`∇nuclear_μ!`](@ref) restricts the potential
`Zc/|r-Rc|` it differentiates to that subset of nuclei.

Depends only on `BS` and `A`, never on the shell pair, so callers looping
over shell pairs for one atom should build it once and pass it through
`∇nuclear!`'s `charges` keyword instead of paying for it per call.
"""
function nuclear_charge_sets(BS::BasisSet{LCint}, A::Int)
    only_A = copy(BS.lib.atm)
    no_A   = copy(BS.lib.atm)
    @inbounds for k in eachindex(BS.atoms)
        if k == A
            no_A[1 + 6*(k-1)] = 0
        else
            only_A[1 + 6*(k-1)] = 0
        end
    end
    return only_A, no_A
end

"""
    ∇nuclear!(out, BS::BasisSet{LCint}, A::Int, i::Int, j::Int; scratch=nothing, charges=nothing)
    ∇nuclear!(out, BS::BasisSet, A)

Mutating counterpart of [`∇nuclear`](@ref): writes into the caller-supplied
`out` instead of allocating. The shell-pair form is the primitive the
full-tensor form builds on. Unlike `∇overlap!`/`∇kinetic!`, no case is ever
a free/skippable zero: `V_ij` sums the potential over every nucleus, so
even a shell pair with neither `i` nor `j` on atom `A` still has a
nonzero derivative through the `Z_A/|r-R_A|` operator term itself moving.

Every case reduces to the same two-term rule. Writing `K_X(p,q)` for
[`∇nuclear_μ!`](@ref) -- the kernel differentiated w.r.t. its FIRST shell
argument, over the nuclei that `X` leaves un-zeroed -- and assigning each
shell a `(sign, charge set)` of `(+, no_A)` when it sits on `A` and
`(-, only_A)` when it does not:

    out = s_i * K_{X_i}(i,j)  +  s_j * K_{X_j}(j,i)ᵀ

(transpose over the two AO axes, per Cartesian direction). Expanding the
four on/off combinations recovers the familiar cases, e.g. `i,j` both on
`A` gives `K_notA(i,j) + K_notA(j,i)ᵀ`, and both off gives
`-K_A(i,j) - K_A(j,i)ᵀ`. Both terms are the *same* primitive, differing only
in argument order and which charge set is passed.

Only implemented for the `LCint` backend -- there is no `ACSint` fallback
for gradients.

# Methods

  - `∇nuclear!(out, BS, A, i, j)`: `out` must be `(Ni,Nj,3)`, the block
    for shells `i,j` of `BS` (shell indices, not AO indices).
  - `∇nuclear!(out, BS, A)`: `out` must be a dense `nbas × nbas × 3`
    array.

In a loop over shell pairs, pass `charges` (from
[`nuclear_charge_sets`](@ref), which depends only on `BS`/`A`) and a
`scratch` vector of at least `3*Ni*Nj` elements; the call is then
allocation-free. Without them it rebuilds both charge arrays and a scratch
buffer every time.
"""
function ∇nuclear!(out, BS::BasisSet{LCint}, A::Int, i::Int, j::Int;
                   scratch=nothing, charges=nothing)
    i_on_A, j_on_A = GaussianBasis.on_atom_flags(BS, A, i, j)
    Ni = num_basis(BS[i])
    Nj = num_basis(BS[j])

    if size(out) != (Ni, Nj, 3)
        throw(DimensionMismatch("Size of the output array needs to be ($Ni, $Nj, 3)"))
    end

    only_A, no_A = charges === nothing ? nuclear_charge_sets(BS, A) : charges

    # Term differentiated w.r.t. shell i -- libcint's [i,j] layout is already
    # out's (Ni,Nj,3), so this writes straight in (and overwrites whatever was
    # there, so a reused `out` needs no zeroing).
    ∇nuclear_μ!(out, BS, i_on_A ? no_A : only_A, i, j)
    i_on_A || (out .*= -1.0)

    # Term differentiated w.r.t. shell j -- same primitive with the arguments
    # swapped, which returns a (Nj,Ni,3) block; transpose it in as we add.
    buf = scratch === nothing ? Vector{Cdouble}(undef, 3*Ni*Nj) : scratch
    length(buf) >= 3*Ni*Nj ||
        throw(DimensionMismatch("scratch must hold at least $(3*Ni*Nj) elements"))
    ∇nuclear_μ!(buf, BS, j_on_A ? no_A : only_A, j, i)
    sj = j_on_A ? 1.0 : -1.0
    @inbounds for k = 1:3, js = 1:Nj, is = 1:Ni
        out[is, js, k] += sj * buf[js + Nj*(is-1) + Nj*Ni*(k-1)]
    end
    return out
end

function ∇nuclear(BS::BasisSet, A::Int, i::Int, j::Int)
    Ni = num_basis(BS[i])
    Nj = num_basis(BS[j])
    out = zeros(Ni, Nj, 3)
    return ∇nuclear!(out, BS, A, i, j)
end

"""
    ∇nuclear(BS::BasisSet, A) -> Array{Float64,3}
    ∇nuclear(BS::BasisSet, A::Int, i::Int, j::Int) -> Array{Float64,3}

Gradient of the AO nuclear attraction matrix `V` w.r.t. atom `A`'s three
Cartesian coordinates, `∂V/∂R_A` (derivative of both the shell centers and
the potential itself, since moving atom `A` also moves its nuclear
charge), with `R_A` in bohr (see [Gradients](@ref) for units).

# Methods

  - `∇nuclear(BS, A)`: full, dense `nbas × nbas × 3` array.
  - `∇nuclear(BS, A, i, j)`: just the `(Ni,Nj,3)` block for shells `i,j`
    of `BS` (shell indices, not AO indices).

For repeated calls, see `∇nuclear!`, which writes into a preallocated
array instead of allocating.
"""
function ∇nuclear(BS::BasisSet, A)
    # Pre allocate output
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇nuclear!(out, BS, A)
end

function ∇nuclear!(out, BS::BasisSet{LCint}, A)

    if size(out) != (BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3)"))
    end

    Nvals = num_basis.(BS.shells)
    ao_offset = cumsum(Nvals) .- Nvals
    Nmax = maximum(Nvals)

    # Depends only on BS/A -- built once here rather than per shell pair.
    only_A, no_A = nuclear_charge_sets(BS, A)

    # Serial, like ∇1e! and for the same reasons (see its comment): a single
    # call is too small to thread profitably, and callers wanting parallelism
    # should thread over atoms, which are independent and write disjoint
    # outputs. Being serial also makes this trivially thread-safe to call
    # from such an outer loop.
    bufi = Vector{Cdouble}(undef, 3*Nmax^2)
    bufj = Vector{Cdouble}(undef, 3*Nmax^2)

    # V is symmetric, so dV_ij/dR_A = (dV_ji/dR_A)^T -- visit each unordered
    # shell pair once and mirror. Every block uses the same two-term rule as
    # the shell-pair form above: one call to the primitive per shell, with the
    # arguments swapped for the second and the charge set chosen by whether
    # that shell sits on A.
    @inbounds for i in 1:BS.nshells
        Ni = Nvals[i]; ioff = ao_offset[i]
        for j in i:BS.nshells
            Nj = Nvals[j]; joff = ao_offset[j]
            i_on_A, j_on_A = GaussianBasis.on_atom_flags(BS, A, i, j)

            ∇nuclear_μ!(bufi, BS, i_on_A ? no_A : only_A, i, j)   # (Ni,Nj,3)
            ∇nuclear_μ!(bufj, BS, j_on_A ? no_A : only_A, j, i)   # (Nj,Ni,3)
            si = i_on_A ? 1.0 : -1.0
            sj = j_on_A ? 1.0 : -1.0

            for k = 1:3
                bi = Ni*Nj*(k-1)
                bj = Nj*Ni*(k-1)
                for js = 1:Nj, is = 1:Ni
                    v = si*bufi[bi + is + Ni*(js-1)] + sj*bufj[bj + js + Nj*(is-1)]
                    out[ioff+is, joff+js, k] = v
                    if i != j
                        out[joff+js, ioff+is, k] = v
                    end
                end
            end
        end
    end

    return out
end


###########################################################
###########################################################
#                     General kernel                       #
###########################################################
###########################################################

# Shared whole-array kernel for overlap/kinetic, parametrized by `callback`
# (the bare, unchecked shell-differentiation primitive, e.g. `∇overlap_μ!`)
# instead of a string, mirroring how `get_1e_matrix!` in
# Integrals/OneElectron.jl takes a callback rather than branching on one.
# Calls the BARE primitive (`∇overlap_μ!`/`∇kinetic_μ!`) rather than the
# safe, bounds-checked `∇overlap!`/`∇kinetic!` shell-pair form: it skips the
# per-pair `on_atom_flags`/size checks this loop has already established, and
# never needs the shell-pair form's transposing branch, because it only ever
# visits pairs with `i` on A -- so the differentiated shell is always
# libcint's first argument and the (j,i) block comes free from the mirror in
# the scatter. Translational invariance also lets it skip same-membership
# pairs entirely and visit Ashells×notAshells once. Safe because `i`/`j` are
# always drawn from this loop's own bounded `Ashells`/`notAshells` and `buf`
# is sized from `Nmax` up front, not because the callback validates anything.
#
# Nuclear isn't wired through here: its per-pair rule needs TWO primitive
# calls with DIFFERENT charge sets (see `∇nuclear!`), which this
# single-callback shape cannot express, and none of its blocks are zero, so
# it has no same-membership pairs to skip either.
function ∇1e(callback, BS::BasisSet, A)
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇1e!(callback, out, BS, A)
end

function ∇1e!(callback, out, BS::BasisSet, A)

    if size(out) != (BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3)"))
    end

    # Blocks with both shells on A (or both off it) are identically zero and
    # are never visited below, so `out` must start clean -- otherwise a reused
    # buffer keeps its old contents there, which is exactly the case this
    # mutating form exists to serve.
    fill!(out, 0.0)

    atomA = BS.atoms[A]

    # Shell indexes for basis in the atom A
    Ashells = Int[]
    notAshells = Int[]
    for i in 1:BS.nshells
        b = BS.shells[i]
        if b.atom == atomA
            push!(Ashells, i)
        else
            push!(notAshells, i)
        end
    end

    Nvals = num_basis.(BS.shells)
    ao_offset = cumsum(Nvals) .- Nvals
    Nmax = maximum(Nvals)

    # Deliberately serial, and thread-safe as a result: a single call is far too
    # small to thread profitably, and its parallelism is capped by |Ashells|
    # (typically 5-15) no matter how many threads exist. Measured medians on 24
    # cores, serial vs a workerpool over Ashells (us):
    #
    #   work     |   216 |   840 |  4872 | 19800 | 25920 | 57240 | 95700
    #   serial   |  12.2 |  30.3 | 156.5 | 297.0 | 489.7 |1657.0 |2102.4
    #   threaded | 132.7 | 104.3 | 251.6 | 333.5 | 490.6 | 862.4 | 958.2
    #
    # (work = output elements written = 6*nbas_on_A*(nbas - nbas_on_A).)
    # Threading only pays past ~350 basis functions, and even then the threaded
    # path is bimodal where the serial one is jitter-free -- minimum-of-N
    # timings flatter it badly. Since a real gradient loops over every atom
    # anyway, and atoms are independent and write disjoint outputs, the useful
    # parallelism lives in that outer loop: callers should thread over atoms
    # and call this per atom.
    buf = Vector{Cdouble}(undef, 3*Nmax^2)
    for i in Ashells
        _scatter_∇1e!(callback, out, BS, i, notAshells, buf, Nvals, ao_offset)
    end
    return out
end

# One shell `i` on atom A against every shell off it: evaluate each pair into
# `buf` and scatter the block plus its (j,i) mirror. Shared by the serial and
# threaded paths above so the two cannot drift apart.
function _scatter_∇1e!(callback, out, BS, i, notAshells, buf, Nvals, ao_offset)
    @inbounds begin
        Ni = Nvals[i]
        ioff = ao_offset[i]
        for j in notAshells
            Nj = Nvals[j]
            joff = ao_offset[j]
            Nij = Ni*Nj
            # Call the bare shell-differentiation primitive. `i` is on A, so
            # it is the differentiated (first) shell and the raw buffer is
            # already (Ni,Nj,3) in memory.
            callback(buf, BS, i, j)

            # Scatter the block and its (j,i) mirror in one scalar pass.
            # The mirror's transpose costs nothing here -- it is just the
            # swapped index expression -- so no permutedims and no
            # range-indexed temporaries are needed.
            for k in 1:3
                base = Nij*(k-1)
                for js = 1:Nj, is = 1:Ni
                    v = buf[base + is + Ni*(js-1)]
                    out[ioff+is, joff+js, k] = v
                    out[joff+js, ioff+is, k] = v
                end
            end
        end
    end
end


###########################################################
###########################################################
#                     General kernel                       #
###########################################################
###########################################################

# Shared whole-array kernel for overlap/kinetic, parametrized by `callback`
# (the bare, unchecked shell-differentiation primitive, e.g. `∇overlap_μ!`)
# instead of a string, mirroring how `get_1e_matrix!` in
# Integrals/OneElectron.jl takes a callback rather than branching on one.
# Calls the BARE primitive (`∇overlap_μ!`/`∇kinetic_μ!`), not the safe,
# bounds-checked `∇overlap!`/`∇kinetic!` shell-pair form -- going through
# the safe form would pay for a fresh `buf` allocation on every call, where
# this loop instead reuses one buffer per worker task across every shell
# pair it handles (translational invariance also lets it skip same-atom
# pairs and only visit Ashells×notAshells once, mirroring the rest via
# transpose). Safe here because `i`/`j` are always drawn from this loop's
# own bounded `Ashells`/`notAshells` and `buf` is sized from `Nmax` up
# front, not because the callback validates anything itself. Nuclear isn't
# wired through here at all: its three-loop structure amortizes the
# `deepcopy`-based fudged nuclear-charge arrays once across every shell
# pair, which this generic per-pair `callback` shape has no way to do.
function ∇1e(callback, BS::BasisSet, A)
    out = zeros(BS.nbas, BS.nbas, 3)
    return ∇1e!(callback, out, BS, A)
end

function ∇1e!(callback, out, BS::BasisSet, A)

    if size(out) != (BS.nbas, BS.nbas, 3)
        throw(DimensionMismatch("Size of the output array needs to be (nbas, nbas, 3)"))
    end

    atomA = BS.atoms[A]

    # Shell indexes for basis in the atom A
    Ashells = Int[]
    notAshells = Int[]
    for i in 1:BS.nshells
        b = BS.shells[i]
        if b.atom == atomA
            push!(Ashells, i)
        else
            push!(notAshells, i)
        end
    end

    Nvals = num_basis.(BS.shells)
    ao_offset = cumsum(Nvals) .- Nvals
    Nmax = maximum(Nvals)

    allocate(body) = body(zeros(Cdouble, 3*Nmax^2))
    workerpool(allocate, Ashells; chunksize = 1) do i, buf
        @inbounds begin
            Ni = Nvals[i]
            ioff = ao_offset[i]
            for j in notAshells
                Nj = Nvals[j]
                joff = ao_offset[j]
                Nij = Ni*Nj
                # Call the bare shell-differentiation primitive
                callback(buf, BS, i, j)
                I = (ioff+1):(ioff+Ni)
                J = (joff+1):(joff+Nj)

                # Get strides for each cartesian
                for k in 1:3
                    r = (1+Nij*(k-1)):(k*Nij)
                    @views ∇k = buf[r]
                    out[I,J,k] .+= reshape(∇k, Ni, Nj)
                end

                # Copy over the transpose
                out[J,I,:] .+= permutedims(out[I,J,:], (2,1,3))
            end
        end #inbounds
    end
    return out
end
