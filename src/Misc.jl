function string_repr(B::SphericalShell)
    # Generate Unicode symbol for sub number
    l_sub = Char(0x2080 + B.l)

    # Unicode for superscript is a bit messier, so gotta use control flow
    l_sup = B.l == 1 ? Char(0x00B9) :
            B.l in [2,3] ? Char(0x00B1 + B.l) :
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