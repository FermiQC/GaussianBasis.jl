export merge_basis, on_atom_flags, atomic_orbital_amplitude

function string_repr(B::SphericalShell)
    # Generate Unicode symbol for sub number
    l_sub = Char(0x2080 + B.l)

    # Unicode for superscript is a bit messier, so gotta use control flow
    l_sup = B.l == 1 ? Char(0x00B9) :
            B.l in [2,3] ? Char(0x00B0 + B.l) :
            Char(0x2070 + B.l)

    nbas = 2*B.l + 1
    mvals = collect(-B.l:B.l)
    nprim = length(B.exp)

    atom_name = Molecules.elements[B.atom.Z].name
    plural_nbas = nbas == 1 ? "" : "s"
    plural_nprim = nprim == 1 ? "" : "s"

    # Reverse Dict(symbol=>num) to get Symbols from B.l
    Lsymbol = Dict(value => key for (key, value) in AMDict)[B.l]
    out = "$(Lsymbol) shell on $atom_name at position $(B.atom.xyz) Å\n"
    out *= "Contains $nbas basis function$plural_nbas built from $nprim primitive gaussian$plural_nprim\n\n"
    for m in mvals
        # Add sub minus sign (0x208B) if necessary
        m_sub = m < 0 ? Char(0x208B)*Char(0x2080 - m) : Char(0x2080 + m)
        out *= format("{:<4s} = ","χ$(l_sub)$m_sub")
        for i in eachindex(B.coef)

            if i > 1
                out *= B.coef[i] > 0 ? "\n     + " : "\n     - "
            end

            #out *= "$(abs(B.coef[i]))⋅Y$(l_sub)$m_sub"
            out *= format("{:>15.10f}⋅Y$(l_sub)$m_sub", abs(B.coef[i]))

            if B.l != 0 
                out *= "⋅r$l_sup"
            end

            out *= "⋅exp(-$(B.exp[i])⋅r²)"
        end
        out *="\n\n"
    end
    return strip(out)
end

function compact_string_repr(B::SphericalShell)
    # Generate Unicode symbol for sub number
    l_sub = Char(0x2080 + B.l)

    # Unicode for superscript is a bit messier, so gotta use control flow
    l_sup = B.l == 1 ? Char(0x00B9) :
            B.l in [2,3] ? Char(0x00B0 + B.l) :
            Char(0x2070 + B.l)

    nbas = 2*B.l + 1
    mvals = collect(-B.l:B.l)
    nprim = length(B.exp)

    atom_name = Molecules.elements[B.atom.Z].name
    plural_nbas = nbas == 1 ? "" : "s"
    plural_nprim = nprim == 1 ? "" : "s"

    # Reverse Dict(symbol=>num) to get Symbols from B.l
    Lsymbol = Dict(value => key for (key, value) in AMDict)[B.l]
    out = "$(Lsymbol) shell on $atom_name ($nbas basis function$plural_nbas, $nprim primitive$plural_nprim)"
    return strip(out)
end

function compact_string_repr(B::CartesianShell)
    nbas = ((B.l + 1) * (B.l + 2)) ÷ 2
    nprim = length(B.exp)

    atom_name = Molecules.elements[B.atom.Z].name
    plural_nbas = nbas == 1 ? "" : "s"
    plural_nprim = nprim == 1 ? "" : "s"

    # Reverse Dict(symbol=>num) to get Symbols from B.l
    Lsymbol = Dict(value => key for (key, value) in AMDict)[B.l]
    out = "$(Lsymbol) shell on $atom_name ($nbas basis function$plural_nbas, $nprim primitive$plural_nprim)"
    return strip(out)
end

function string_repr(B::CartesianShell)
    # Generate Unicode symbol for sub number
    l_sub = Char(0x2080 + B.l)

    # Unicode for superscript is a bit messier, so gotta use control flow
    l_sup(l) = l < 2 ? "" :
            l in [2,3] ? Char(0x00B0 + l) :
            Char(0x2070 + l)

    nbas = ((B.l + 1) * (B.l + 2)) ÷ 2
    mvals = String[]
    for a = B.l:-1:0
        for b = B.l:-1:0
            c = B.l - a - b
            if c < 0 
                continue
            end
            r_str  = a > 0 ? "x"*l_sup(a) : ""
            r_str *= b > 0 ? "y"*l_sup(b) : ""
            r_str *= c > 0 ? "z"*l_sup(c) : ""
            if !isempty(r_str)
                r_str *= "⋅"
            end
            push!(mvals, r_str)
        end
    end
    nprim = length(B.exp)

    atom_name = Molecules.elements[B.atom.Z].name
    plural_nbas = nbas == 1 ? "" : "s"
    plural_nprim = nprim == 1 ? "" : "s"

    # Reverse Dict(symbol=>num) to get Symbols from B.l
    Lsymbol = Dict(value => key for (key, value) in AMDict)[B.l]
    out = "$(Lsymbol) shell on $atom_name at position $(B.atom.xyz) Å\n"
    out *= "Contains $nbas basis function$plural_nbas built from $nprim primitive gaussian$plural_nprim\n\n"
    for m in mvals
        χ = "χ"
        if !isempty(m)
            _m = replace(m, "⋅"=>"")
            χ *= "($_m)"
        end
        out *= format("{:<4s} = ",χ)
        for i in eachindex(B.coef)

            if i > 1
                out *= B.coef[i] > 0 ? "\n     + " : "\n     - "
            end

            #out *= "$(abs(B.coef[i]))⋅Y$(l_sub)$m_sub"
            out *= format("{:>15.10f}⋅$m", abs(B.coef[i]))

            out *= "exp(-$(B.exp[i])⋅r²)"
        end
        out *="\n\n"
    end
    return strip(out)
end

function string_repr(B::BasisSet{T,Y,P}) where {T,Y,P}
    out  =  "$(B.name) Basis Set\n"
    out *= "Type: "*(P <: SphericalShell ? "Spherical" : "Cartesian")
    out *= "   Backend: "
    if T === GaussianBasis.LCint
        out *= "Libcint\n\n"
    else
        out *= replace("$(T)\n\n", "GaussianBasis."=>"")
    end
    out *= "Number of shells: $(B.nshells)\n"
    out *= "Number of basis:  $(B.nbas)\n\n"

    l_to_symbol = Dict(
        0 => "s",
        1 => "p",
        2 => "d",
        3 => "f",
        4 => "g",
        5 => "h",
        6 => "i",
    )

    for i in eachindex(B.atoms)
        A = B.atoms[i]
        # Count how many times s,p,d appears for numbering
        count = zeros(Int16, 7)
        out *= "$(symbol(A)): "
        for b in B.shells
            if b.atom == A
                L = l_to_symbol[b.l]
                count[b.l+1] += 1
                out *= "$(count[b.l+1])$(L) "
            end
        end
        out *="\n"
    end

    return strip(out)
end

# Pretty printing
function show(io::IO, ::MIME"text/plain", X::T) where T<:Union{ShellFunction, BasisSet}
    print(io, string_repr(X))
end

function show(io::IO, X::T) where T<:ShellFunction
    print(io, compact_string_repr(X))
end

"""
    THREADING_THRESHOLD_1E

Work below which the one-electron/multipole full-matrix builders
(`get_1e_matrix!`, `get_multipole_matrix!`) run serially instead of going
through [`workerpool`](@ref). Work is measured as the number of output
elements written (`nbas^2`, times `3^N` Cartesian components for multipoles).

These integrals are cheap enough that spawning a task per thread and routing
shell pairs through a `Channel` can cost more than the integrals themselves.
Measured `overlap!` on a 24-core machine, serial vs threaded (ms):

| `nbas^2` |  49   |  169  |  576  |  3364 | 13456 | 53824 |
|:---------|------:|------:|------:|------:|------:|------:|
| serial   | 0.009 | 0.021 | 0.051 | 0.148 | 0.579 | 2.156 |
| 24 tasks | 0.088 | 0.074 | 0.083 | 0.111 | 0.257 | 0.814 |

The threaded path costs a roughly fixed ~0.07-0.08 ms of pool setup, so it
only pays off once the work itself exceeds that -- crossover sits between
576 and 3364, hence the value here. Below it the pool made a small-molecule
`overlap!` nearly 10x slower than just doing the work.

The two-electron builders are far more expensive per shell tuple and scale
properly, so they always thread.
"""
const THREADING_THRESHOLD_1E = 2_000


# adapted from https://juliafolds.github.io/data-parallelism/tutorials/concurrency-patterns/
function workerpool(work!, allocate, inputs; chunksize,ntasks = Threads.nthreads())
    requests = Channel{Vector{eltype(inputs)}}(Inf)
    for chunk in Iterators.partition(inputs, chunksize)
        put!(requests, chunk)
    end
    close(requests)

    @sync for _ in 1:ntasks
        Threads.@spawn allocate() do resource
            for chunk in requests
                for input in chunk
                    work!(input, resource)
                end
            end
        end
    end
end

"""
    merge_basis(BS1::BasisSet, BS2::BasisSet) -> BasisSet

The concatenation of `BS1` and `BS2` into one `BasisSet` -- regular shells
first, auxiliary second -- as libcint's 3-center kernels require, since they
resolve all three shell indices against a single basis. An auxiliary shell
`k` of `BS2` is therefore addressed as `k + BS1.nshells` in the result.

Depends only on `BS1`/`BS2`, never on atoms or shells, so routines that take
a `Bmerged` keyword should be handed one built once outside the loop.
"""
function merge_basis(BS1::BasisSet, BS2::BasisSet)
    atoms = unique(vcat(BS1.atoms, BS2.atoms))
    basis = vcat(BS1.shells, BS2.shells)
    return BasisSet("$(BS1.name*BS2.name)", atoms, basis)
end

"""
    unique_ij(nshells::Integer) -> Vector{NTuple{2,Int16}}

All shell pairs `(i,j)` with `1 <= i <= j <= nshells` -- the canonical
upper-triangle worklist that `workerpool`-based full-tensor builders
(`get_1e_matrix!`, `get_multipole_matrix!`) iterate over, relying on
`X_ij = X_ji^T` symmetry to compute each unordered pair once and mirror it.

Ordered so that the pair at position `n` is exactly the one with
`index2(i-1, j-1) == n - 1` (i.e. `j` outer, `i` inner). `sparseERI_2e4c`
depends on that: it indexes its Cauchy-Schwarz `σvals` by the same
composite index, so `σvals[n]` pairs with `unique_ij(...)[n]`.

Pre-sized and filled directly rather than `collect`-ed from a lazy
generator, so there's a single allocation instead of the repeated grow/copy
`collect` does when it can't know the final length upfront. Elements are
`Int16` for the same reason as [`unique_ijkl`](@ref) -- 4x smaller than
`Int` at equal length, and shell counts never approach `typemax(Int16)`.
"""
function unique_ij(nshells::Integer)
    pairs = Vector{NTuple{2,Int16}}(undef, (nshells*(nshells+1)) ÷ 2)
    n = 0
    for j = one(Int16):Int16(nshells), i = one(Int16):j
        n += 1
        pairs[n] = (i, j)
    end
    return pairs
end

"""
    unique_ijkl(nshells::Integer) -> Vector{NTuple{4,Int16}}

All shell quartets `(i,j,k,l)` (1-based) that are unique under the 8-fold
permutational symmetry of `(ij|kl)`: `i <= j`, `k <= l`, and the composite
bra index not below the ket index. This is the worklist `ERI_2e4c!`'s
dense full-tensor build iterates over, computing each unique quartet once
and scattering it to its symmetry images.

Like [`unique_ij`](@ref), pre-sized rather than grown with `push!` -- the
exact count is closed-form, `npair*(npair+1)/2` with
`npair = nshells*(nshells+1)/2`, so no counting pass is needed. Elements
are `Int16` (not `Int`): this list has `O(nshells^4)` entries and is the
single largest auxiliary allocation in the dense build, so the narrower
element type matters (4x smaller at equal length).
"""
function unique_ijkl(nshells::Integer)
    npair = (nshells*(nshells+1)) ÷ 2
    quartets = Vector{NTuple{4,Int16}}(undef, (npair*(npair+1)) ÷ 2)
    n = 0
    N = Int16(nshells)
    for i = one(Int16):N, j = i:N, k = one(Int16):N, l = k:N
        index2(i-one(Int16), j-one(Int16)) < index2(k-one(Int16), l-one(Int16)) && continue
        n += 1
        quartets[n] = (i, j, k, l)
    end
    return quartets
end

"""
    on_atom_flags(BS::BasisSet, iA::Int, shells...)

Whether each of `shells` (any number of shell indices) sits on atom `iA`, as
an `NTuple{N,Bool}`.
"""
function on_atom_flags(BS::BasisSet, iA::Int, shells...)
    A = BS.atoms[iA]
    return ntuple(p -> BS.shells[shells[p]].atom === A, length(shells))
end


function legendre_polynomial(m::Integer,l::Integer,x)
    fmm = (-1)^m * prod((2m-1):-2:1)
    Pmm = fmm * (1 - x^2)^(m/2)
    if l == m
        return Pmm
    end
    Pmm1 = x * (2m+1) * Pmm

    if l < m || m < 0
        throw(DomainError("Associated Legendre polynomials are only defined for 0 ≤ m ≤ l, got m = $m, l = $l"))
    end

    Pmm2 = zero(Pmm)
    for i = m+2:l
        Pmm2 = ((2i-1)*x*Pmm1 - (i-1+m) * Pmm)/(i-m)
        Pmm = Pmm1
        Pmm1 = Pmm2
    end
    return Pmm1
end


function atomic_orbital_angular_part(shell::SphericalShell, n::Integer, d, d²)
    l = shell.l
    m = n - l - 1

    # for p shell, ordering is px, py, pz <=> m = 1,-1,0
    if l == 1
        m = (1, -1, 0)[n]
    end

    sqrtd² = sqrt.(d²)
    cosθ = map(selectdim(d, 1, 3), sqrtd²) do z, sqrtd²
        ifelse(iszero(z), zero(z), z / sqrtd²)
    end
    Nlm = (-1)^m * sqrt((2l + 1) / (4π * prod(l - abs(m) + 1:l + abs(m))))
    Plm = legendre_polynomial.(abs(m), l, cosθ)

    v = @. sqrtd²^l * Nlm * Plm
    if m != 0
        φ = [atan(x[2],x[1]) for x in eachslice(d, dims=tuple(2:ndims(d)...))]
        if m > 0
            @. v *= sqrt(2) * cos(m * φ)
        else
            @. v *= sqrt(2) * sin(-m * φ)
        end
    end

    return v
end


function atomic_orbital_amplitude(shell::ShellFunction, n::Integer, r::AbstractArray)
    T = promote_type(eltype(shell.coef), eltype(shell.exp), eltype(r), eltype(shell.atom.xyz))

    if size(r,1) != 3
        throw(DimensionMismatch("First dimension of r must be 3, got $(size(r,1))"))
    end

    v = zeros(T, size(r)[2:end])

    d = r .- shell.atom.xyz ./ Molecules.bohr_to_angstrom

    d² = dropdims(sum(x->x^2, d; dims=1); dims=1)

    for (c, e) in zip(shell.coef, shell.exp)
        @. v += c * exp(-e * d²)
    end

    v .*= atomic_orbital_angular_part(shell, n, d, d²)
    return v
end

# Returns the shell and the index of the basis function within that shell for a given basis function index in the basis set.
# i.e., think for sto-3g water [O, H, H], find_shell_m(bset, 4) returns the O p shell and the index 2.
function find_shell_m(B::BasisSet, idx::Integer)
    if idx <= 0
        throw(BoundsError(B, idx))
    end

    n = idx
    for shell in B.shells
        if n - (2shell.l +1) <= 0
            return shell, n
        end
        n -= 2shell.l + 1
    end

    throw(BoundsError(B, idx))
end

"""
    atomic_orbital_amplitude(B::BasisSet, i::Integer, r::AbstractVector) -> Real
    atomic_orbital_amplitude(B::BasisSet, i::Integer, r::AbstractArray) -> AbstractArray

Returns the amplitude of the `i`th basis function at a set of real space positions. `r` can be either a 3-vector containing a single position or a 3 × … array of positions for evaluating many points at once more efficiently.

# Examples

```julia-repl
julia> bset = BasisSet("sto-3g", "H 0 0 0");
julia> atomic_orbital_amplitude(bset, 1, [0.0,0.0,0.0])
0.6282468778403579
julia> atomic_orbital_amplitude(bset, 1, [0.1;0.2;0.3;;-0.1;0.3;-0.2])
2-element Vector{Float64}:
 0.49840190793869554
 0.49840190793869554
```
"""
function atomic_orbital_amplitude(B::BasisSet, idx::Integer, r::AbstractArray)
    return atomic_orbital_amplitude(find_shell_m(B, idx)..., r)
end
# This dispatch supports a vector instead of matrix, it calls the method above.
function atomic_orbital_amplitude(B::BasisSet, idx::Integer, r::AbstractVector)
    return atomic_orbital_amplitude(B, idx, reshape(r,:,1))[1]
end
