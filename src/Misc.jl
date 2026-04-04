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

    # Reverse Dict(symbol=>num) to get Symbols from B.l
    Lsymbol = Dict(value => key for (key, value) in AMDict)[B.l]
    out = "$(Lsymbol) shell with $nbas basis built from $nprim primitive gaussians\n\n"
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

    # Reverse Dict(symbol=>num) to get Symbols from B.l
    Lsymbol = Dict(value => key for (key, value) in AMDict)[B.l]
    out = "$(Lsymbol) shell with $nbas basis built from $nprim primitive gaussians\n\n"
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
    out *= "Type: "*replace("$(P)", "Shell"=>"")
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
        for b in B.basis
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
function show(io::IO, ::MIME"text/plain", X::T) where T<:Union{BasisFunction, BasisSet}
    print(io, string_repr(X))
end

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


function atomic_orbital_amplitude(shell::BasisFunction, n::Integer, r::AbstractArray)
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

function find_shell_m(B::BasisSet, idx::Integer)
    if idx <= 0
        throw(BoundsError(B, idx))
    end

    n = idx
    for shell in B.basis
        if n - (2shell.l +1) <= 0
            m = n - shell.l - 1

            # for p shell, ordering is px, py, pz <=> m = 1,-1,0
            if shell.l == 1
                m = (1, -1, 0)[n]
            end

            return shell, n
        end
        n -= 2shell.l + 1
    end

    throw(BoundsError(B, idx))
end

"""
    atomic_orbital_value(B::BasisSet, i::Integer, r::AbstractVector) -> Real
    atomic_orbital_value(B::BasisSet, i::Integer, r::AbstractArray) -> AbstractArray

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
function atomic_orbital_amplitude(B::BasisSet, idx::Integer, r::AbstractVector)
    return atomic_orbital_amplitude(B, idx, reshape(r,:,1))[1]
end
