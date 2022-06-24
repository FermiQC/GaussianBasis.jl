const LIBPATH = joinpath(@__DIR__, "../lib")
const AM_pat = r"([SPDFGHI]{1,2})\s+?(\d++)"
const prim_pat = r"([+-]?\d*?\.\d+[D+-]{0,2}\d\d)\s+?([+-]?\d*?\.\d+[D+-]{0,2}+\d\d)"
const prim_pat3 = r"([+-]?\d*?\.\d+[D+-]{0,2}\d\d)\s+?([+-]?\d*?\.\d+[D+-]{0,2}+\d\d)\s+?([+-]?\d*?\.\d+[D+-]{0,2}+\d\d)"
const AMDict = Dict(
        "S" => 0,
        "P" => 1, 
        "D" => 2, 
        "F" => 3, 
        "G" => 4, 
        "H" => 5, 
        "I" => 6, 
    )

"""
    GaussinBasis.read_basisset(bname::String, atom::atom; spherical=true)

Returns an array of BasisFunction objects for the given `Atom` and basis set name (bname).
By default, functions are normalized as Spherical functions (spherical=true). If spherical is set to false,
Cartesian functions are returned instead. 
"""
function read_basisset(bname::String, atom::A; spherical=true) where A <: Atom

    AtomSymbol = Molecules.symbol(atom) 

    # Transform basis name to file name e.g. 6-31g* => 6-31g_st_
    clean_bname = replace(bname, "*"=>"_st_")
    file_path = joinpath(LIBPATH, clean_bname*".gbs")

    if !(isfile(file_path))
        throw(ArgumentError("Basis set file for $bname was not found."))
    end

    # Fetch the block of the basis set file regarding the given atom 
    info = extract_atom_from_bs(file_path, AtomSymbol)

    # Split info into strings for individual basis functions
    BasisStrings = [info[1]*"\n"]
    for line in info[2:end]
        if occursin(AM_pat, line)
            push!(BasisStrings, line*"\n")
        else
            BasisStrings[end] *= line * "\n"
        end
    end

    # Parse each string into a basis function object
    out = BasisFunction[]
    for b in BasisStrings
        r = r"[SPDFGHI]{2}"
        if occursin(r, b)
            bf1, bf2 = two_basis_from_string(b, atom, spherical=spherical)
            push!(out, bf1)
            push!(out, bf2)
        else
            bf = basis_from_string(b, atom, spherical=spherical)
            push!(out, bf)
        end
    end

    return out
end

"""
    GaussinBasis.extract_atom_from_bs(file_path::String, AtomSymbol::String)

Extract the portion of the basis set file (file_path) relevant for the given atom in AtomSymbol.
"""
function extract_atom_from_bs(file_path::String, AtomSymbol::String)

    atom_pat = Regex("$(AtomSymbol)\\s+?0")
    # This flag indicates if the atom was found within the file
    flag_atom = false

    # Output is an array of strings (corresponding to the lines of the block matching the desired atom)
    out = String[]

    # Iterate through the lines of the file until the atom is found
    for line in eachline(file_path)

        if occursin(atom_pat, line)
            flag_atom = true
            continue
        end

        if flag_atom
            occursin(raw"****", line) ? break : nothing
            push!(out, String(line))
        end
    end

    if flag_atom
        return out
    else
        throw(FermiException("Atom $AtomSymbol not found in $file_path."))
    end
end

"""
    GaussinBasis.basis_from_string(bstring::String, atom::A; spherical=true) where A <: Atom

From a String block representing two Basis Function (e.g. SP blocks) in gbs format, returns a
`BasisFunction` object. For the special case of two basis being described within the 
same block (e.g. SP blocks) see `two_basis_from_string`
"""
function basis_from_string(bstring::String, atom::A; spherical=true) where A <: Atom
    lines = split(strip(bstring), "\n")
    head = lines[1]
    m = match(AM_pat, head)

    if m === nothing
        throw(ArgumentError("cannot parse basis set file line: $head"))
    end

    AMsymbol = String(m.captures[1])
    l = AMDict[AMsymbol]
    nprim = parse(Int, m.captures[2])
    coef = zeros(nprim)
    exp = zeros(nprim)

    for i in eachindex(exp)
        m = match(prim_pat, lines[i+1])

        if m === nothing
            throw(ArgumentError("cannot parse basis set file line: $(lines[i+1])"))
        end

        e = replace(m.captures[1], "D"=>"E")
        c = replace(m.captures[2], "D"=>"E")
        coef[i] = parse(Float64, c)
        exp[i] = parse(Float64, e)
    end

    if spherical
        normalize_spherical!(coef, exp, l)
        coef = SVector{nprim}(coef)
        exp = SVector{nprim}(exp)
        return SphericalShell(l, coef, exp, atom)
    else
        normalize_cartesian!(coef, exp, l)
        coef = SVector{nprim}(coef)
        exp = SVector{nprim}(exp)
        return CartesianShell(l, coef, exp, atom)
    end
end

"""
    GaussinBasis.two_basis_from_string(bstring::String, atom::A; spherical=true) where A <: Atom

From a String block representing two Basis Function (e.g. SP blocks) in gbs format, returns a
`BasisFunction` object. For the case of a single basis being described within the block see `basis_from_string`
"""
function two_basis_from_string(bstring::String, atom::A; spherical=true) where A <: Atom
    lines = split(strip(bstring), "\n")
    head = lines[1]
    m = match(AM_pat, head)
    if m === nothing
        throw(ArgumentError("cannot parse basis set file line: $head"))
    end
    AMsymbol = String(m.captures[1])

    length(AMsymbol) == 2 || throw(FermiException("cannot extract two basis from $AMsymbol function"))

    l1 = AMDict[AMsymbol[1]*""]
    l2 = AMDict[AMsymbol[2]*""]

    nprim = parse(Int, m.captures[2])
    coef1 = zeros(nprim)
    coef2 = zeros(nprim)
    exp = zeros(nprim)

    for i in eachindex(exp)
        m = match(prim_pat3, lines[i+1])

        if m === nothing
            throw(ArgumentError("cannot parse basis set file line: $(lines[i+1])"))
        end

        e = replace(m.captures[1], "D"=>"E")
        c1 = replace(m.captures[2], "D"=>"E")
        c2 = replace(m.captures[3], "D"=>"E")
        coef1[i] = parse(Float64, c1)
        coef2[i] = parse(Float64, c2)
        exp[i] = parse(Float64, e)
    end

    if spherical
        normalize_spherical!(coef1, exp, l1)
        normalize_spherical!(coef2, exp, l2)
        coef1 = SVector{nprim}(coef1)
        coef2 = SVector{nprim}(coef2)
        exp   = SVector{nprim}(exp)
        return SphericalShell(l1, coef1, exp, atom), SphericalShell(l2, coef2, exp, atom)
    else
        normalize_cartesian!(coef1, exp, l1)
        normalize_cartesian!(coef2, exp, l2)
        coef1 = SVector{nprim}(coef1)
        coef2 = SVector{nprim}(coef2)
        exp   = SVector{nprim}(exp)
        return CartesianShell(l1, coef1, exp, atom), CartesianShell(l2, coef2, exp, atom)
    end
end