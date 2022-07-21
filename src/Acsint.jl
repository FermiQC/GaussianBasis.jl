module Acsint

using GaussianBasis
using GaussianBasis: index2, num_basis
using Molecules
using SpecialFunctions
using Combinatorics: doublefactorial
using LinearAlgebra: norm, eigen
using StaticArrays

export generate_ERI_quartet!, generate_S_pair!, generate_T_pair!, generate_V_pair!

const ang2bohr = 1.8897261246257702

function cumulative_cart_dim(L)
    # The number of Cartesian components in all shell with angular momentum L and lower """
    ((L + 1) * (L + 2) * (L + 3)) ÷ 6
end

function cart_dim(L)
    # The number of Cartesian components in a shell with angular momentum L
    ((L + 1) * (L + 2)) ÷ 2
end

function addr(t::Int, u::Int, v::Int, l::Int)
    # The address in a supermatrix with angular momentum containing all components from zero
    # to t+u+v=l, where the summation runs with t then u then v being fastest running
    ( v * (11 + 3 * l * (4 + l) - 6 * v - 3 * l * v + v * v) + 3 * u * (3 + 2 * l - 2 * v - u)) ÷ 6 + t + 1
end

function generate_S_pair!(out, BS::BasisSet, s1, s2)
    generate_S_pair!(out, BS.basis[s1], BS.basis[s2])
end

function generate_S_pair!(out, B1::CartesianShell, B2::CartesianShell)
    # Generates the block of the overlap matrix corresponding to a given shell pair

    out .= 0.0

    am1 = B1.l
    am2 = B2.l

    A = B1.atom.xyz .* ang2bohr
    B = B2.atom.xyz .* ang2bohr

    N1 = num_basis(B1)

    Ex = similar(A, am1 + 1, am2 + 1, am1 + am2 + 2)
    Ey = similar(A, am1 + 1, am2 + 1, am1 + am2 + 2)
    Ez = similar(A, am1 + 1, am2 + 1, am1 + am2 + 2)
    Ex .= 0.0
    Ey .= 0.0
    Ez .= 0.0

    for (ca, a) in zip(B1.coef, B1.exp)
        for (cb, b) in zip(B2.coef, B2.exp)
            p = a + b
            P = (a * A + b * B) / p
            prefac = ca * cb * (π / p) ^ 1.5

            generate_E_matrix!(Ex, Ey, Ez, am1, am2, P, A, B, a, b)

            index1 = 1
            for l1 in am1 : -1 : 0
                for n1 in 0 : am1 - l1
                    m1 = am1 - l1 - n1
                    index2 = 1
                    for l2 in am2 : -1 : 0
                        for n2 in 0 : am2 - l2
                            m2 = am2 - l2 - n2
                            out[index1 + N1*(index2-1)] += Ex[l1+1, l2+1, 1] * Ey[m1+1, m2+1, 1] * Ez[n1+1, n2+1, 1] * prefac
                            index2 += 1
                        end
                    end
                    index1 += 1
                end
            end
        end
    end
    return out
end

function generate_T_pair!(out, BS::BasisSet, s1, s2)
    generate_T_pair!(out, BS.basis[s1], BS.basis[s2])
end

function generate_T_pair!(out, B1::CartesianShell, B2::CartesianShell)

    # Generates the block of the overlap matrix corresponding to a given shell pair
    out .= 0.0

    am1 = B1.l
    am2 = B2.l

    A = B1.atom.xyz .* ang2bohr
    B = B2.atom.xyz .* ang2bohr

    N1 = num_basis(B1)

    Ex = similar(A, am1 + 3, am2 + 1, am1 + am2 + 4)
    Ey = similar(A, am1 + 3, am2 + 1, am1 + am2 + 4)
    Ez = similar(A, am1 + 3, am2 + 1, am1 + am2 + 4)

    Ex .= 0.0
    Ey .= 0.0
    Ez .= 0.0
    for (ca, a) in zip(B1.coef, B1.exp)
        for (cb, b) in zip(B2.coef, B2.exp)
            p = a + b
            P = (a * A + b * B) / p
            prefac = -0.5 * ca * cb * (π / p) ^ 1.5

            generate_E_matrix!(Ex, Ey, Ez, am1+2, am2, P, A, B, a, b)

            index1 = 1
            for l1 in am1 : -1 : 0
                for n1 in 0 : am1 - l1
                    m1 = am1 - l1 - n1
                    index2 = 1
                    for l2 in am2 : -1 : 0
                        for n2 in 0 : am2 - l2
                            m2 = am2 - l2 - n2
                            Sx0 = Ex[l1+1, l2+1, 1]
                            Sy0 = Ey[m1+1, m2+1, 1]
                            Sz0 = Ez[n1+1, n2+1, 1]
                            # 9.3.31
                            Dx2 = 2 * a * (2 * a * Ex[l1+3, l2+1, 1] - (2 * l1 + 1) * Ex[l1+1, l2+1, 1])
                            if l1 > 1
                                Dx2 += l1 * (l1 - 1) * Ex[l1-1, l2+1, 1]
                            end
                            Dy2 = 2 * a * (2 * a * Ey[m1+3, m2+1, 1] - (2 * m1 + 1) * Ey[m1+1, m2+1, 1])
                            if m1 > 1
                                Dy2 += m1 * (m1 - 1) * Ey[m1-1, m2+1, 1]
                            end
                            Dz2 = 2 * a * (2 * a * Ez[n1+3, n2+1, 1] - (2 * n1 + 1) * Ez[n1+1, n2+1, 1])
                            if n1 > 1
                                Dz2 += n1 * (n1 - 1) * Ez[n1-1, n2+1, 1]
                            end
                            # eqn. 9.3.37
                            out[index1 + N1*(index2-1)] += prefac * (Dx2 * Sy0 * Sz0
                                                             + Sx0 * Dy2 * Sz0
                                                             + Sx0 * Sy0 * Dz2)
                            index2 +=1
                        end
                    end
                    index1 += 1
                end
            end
        end
    end
    return out
end

function generate_V_pair!(out, BS::BasisSet, s1, s2)
    generate_V_pair!(out, BS.basis[s1], BS.basis[s2], BS.atoms)
end

function generate_V_pair!(out, B1::CartesianShell, B2::CartesianShell, atoms::Vector{<:Atom})

    out .= 0.0

    am1 = B1.l
    am2 = B2.l
    am = am1 + am2

    A = B1.atom.xyz .* ang2bohr
    B = B2.atom.xyz .* ang2bohr

    N1 = num_basis(B1)

    Rmat = similar(A, am + 1, cumulative_cart_dim(am))
    ints = similar(A, cart_dim(am1), cart_dim(am2))
    Rmat .= 0.0
    ints .= 0.0

    Ex = similar(A, am1 + 1, am2 + 1, am1 + am2 + 2)
    Ey = similar(A, am1 + 1, am2 + 1, am1 + am2 + 2)
    Ez = similar(A, am1 + 1, am2 + 1, am1 + am2 + 2)

    Ex .= 0.0
    Ey .= 0.0
    Ez .= 0.0
    for At in 1:length(atoms)
        Z = atoms[At].Z
        C = atoms[At].xyz .* ang2bohr
        for (ca, a) in zip(B1.coef, B1.exp)
            for (cb, b) in zip(B2.coef, B2.exp)
                p = a + b
                P = (a * A + b * B) / p
                generate_E_matrix!(Ex, Ey, Ez, am1, am2, P, A, B, a, b)
                prefac = -Z * ca * cb * 2 * π / p

                generate_R_matrix!(Rmat, am, p, P, C)
                index1 = 1
                for l1 in am1 : -1 : 0
                    for n1 in 0 : am1 - l1
                        m1 = am1 - l1 - n1
                        index2 = 1
                        for l2 in am2 : -1 : 0
                            for n2 in 0 : am2 - l2
                                m2 = am2 - l2 - n2
                                val = 0.0
                                for t = 0 : l1 + l2
                                    for u = 0 : m1 + m2
                                        for v = 0 : n1 + n2
                                            Eab = Ex[l1+1, l2+1, t+1] * Ey[m1+1, m2+1, u+1] * Ez[n1+1, n2+1, v+1]
                                            # eqn. 9.9.32
                                            val += Eab * Rmat[1, addr(t, u, v, am)]
                                        end
                                    end
                                end
                                out[index1 + N1*(index2-1)] += prefac * val
                                index2 += 1
                            end
                        end
                        index1 += 1
                    end
                end
            end
        end
    end
    return out
end

function generate_quanta_list(am::Int)
    """ A cached version of the generator, because generators seem to be painfully slow """
    vals = Vector{Tuple{Int, Int, Int, Int}}(undef, cart_dim(am))
    index = 0
    for l in am : -1 : 0
        for n in 0 : am - l
            m = am - l - n
            index += 1
            vals[index] = (l, m, n, index)
        end
    end
    return vals 
end

function generate_ERI_quartet!(out, BS::BasisSet, s1, s2, s3, s4, α::Float64=1.0, β::Float64=0.0, ω::Float64=0.0)
    generate_ERI_quartet!(out, BS.basis[s1], BS.basis[s2], BS.basis[s3], BS.basis[s4], α, β, ω)
end

function generate_ERI_quartet!(out, B1::CartesianShell, B2::CartesianShell, B3::CartesianShell, B4::CartesianShell,
                              α::Float64=1.0, β::Float64=0.0, ⍵::Float64=0.0)

    """ Generates generalized Coulomb integrals of the form
        α + β(erf ⍵ R)
        --------------
              R       
        for a given shell quartet """

    out .= 0.0

    atom1 = B1.atom
    atom2 = B2.atom
    atom3 = B3.atom
    atom4 = B4.atom

    T = eltype(atom1.xyz)
    am1 = B1.l
    am2 = B2.l
    am3 = B3.l
    am4 = B4.l
    am = am1 + am2 + am3 + am4

    N1   = num_basis(B1)
    N12  = N1*num_basis(B2)
    N123 = N12*num_basis(B3)
    #ints = zeros(T, cart_dim(am1), cart_dim(am2), cart_dim(am3), cart_dim(am4))


    A = atom1.xyz .* ang2bohr
    B = atom2.xyz .* ang2bohr
    C = atom3.xyz .* ang2bohr
    D = atom4.xyz .* ang2bohr

    Rmat = zeros(T, am + 1, cumulative_cart_dim(am))

    bra_Ex = zeros(T, am1 + 1, am2 + 1, am1 + am2 + 2)
    bra_Ey = zeros(T, am1 + 1, am2 + 1, am1 + am2 + 2)
    bra_Ez = zeros(T, am1 + 1, am2 + 1, am1 + am2 + 2)
    ket_Ex = zeros(T, am3 + 1, am4 + 1, am3 + am4 + 2)
    ket_Ey = zeros(T, am3 + 1, am4 + 1, am3 + am4 + 2)
    ket_Ez = zeros(T, am3 + 1, am4 + 1, am3 + am4 + 2)

    iter1 = generate_quanta_list(am1)
    iter2 = generate_quanta_list(am2)
    iter3 = generate_quanta_list(am3)
    iter4 = generate_quanta_list(am4)

    for (ca, a) in zip(B1.coef, B1.exp)
        for (cb, b) in zip(B2.coef, B2.exp)
            p = a + b
            P = (a * A + b * B) / p
            generate_E_matrix!(bra_Ex, bra_Ey, bra_Ez, am1, am2, P, A, B, a, b)
            for (cc, c) in zip(B3.coef, B3.exp)
                for (cd, d) in zip(B4.coef, B4.exp)
                    q = c + d
                    Q = (c * C + d * D) / q
                    generate_E_matrix!(ket_Ex, ket_Ey, ket_Ez, am3, am4, Q, C, D, c, d)

                    θ = p * q / (p + q)
                    prefac = 2 * π^(5/2) * ca * cb * cc * cd / (p * q * sqrt(p + q))

                    generate_R_matrix!(Rmat, am, θ, P, Q, α, β, ⍵)
                    for (l1, m1, n1, index1) in iter1
                        for (l2, m2, n2, index2) in iter2
                            _ij = index1 + N1*(index2-1) 
                            for t = 0 : l1 + l2
                                for u = 0 : m1 + m2
                                    for v = 0 : n1 + n2
                                        @inbounds Eab = bra_Ex[l1+1, l2+1, t+1] * bra_Ey[m1+1, m2+1, u+1] * bra_Ez[n1+1, n2+1, v+1]
                                        for (l3, m3, n3, index3) in iter3
                                            _ijk = _ij + N12*(index3-1)
                                            for (l4, m4, n4, index4) in iter4
                                                idx = _ijk + N123*(index4-1)
                                                val = 0.0
                                                for τ = 0 : l3 + l4
                                                    for ν = 0 : m3 + m4
                                                        for φ = 0 : n3 + n4
                                                            @inbounds @fastmath Ecd = ket_Ex[l3+1, l4+1, τ+1] * ket_Ey[m3+1, m4+1, ν+1] * ket_Ez[n3+1, n4+1, φ+1]
                                                            # eqn. 9.9.33
                                                            if isodd(τ+ν+φ)
                                                                @inbounds @fastmath val -= Ecd * Rmat[1, addr(t + τ, u + ν, v + φ, am)]
                                                            else
                                                                @inbounds @fastmath val += Ecd * Rmat[1, addr(t + τ, u + ν, v + φ, am)]
                                                            end
                                                        end
                                                    end
                                                end
                                                @inbounds out[idx] += prefac * Eab * val
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return out
end

@inline function generate_E_matrix!(Ex, Ey, Ez, maxam1, maxam2,
                           P, A, B,
                           a::Float64, b::Float64)
    # Makes the Hermite->Cartesian conversion factors by recursion
    Ex .= 0
    Ey .= 0
    Ez .= 0

    # 9.2.11
    p = a + b
    # 9.2.14
    AB = A - B
    PA = P - A
    PB = P - B
    # 9.2.12
    µ = a * b / p
    # 9.2.15
    @inbounds begin
        Ex[1, 1, 1] = exp(-µ * AB[1]^2)
        Ey[1, 1, 1] = exp(-µ * AB[2]^2)
        Ez[1, 1, 1] = exp(-µ * AB[3]^2)
    end
    oo2p = 1 / (2 * p)
    @inbounds for j = 1 : maxam2 + 1
        for i = 1 : maxam1 + 1
            # handle t = 1 case
            if i > 1
                Ex[i, j, 1] += PA[1] * Ex[i - 1, j, 1] + Ex[i - 1, j, 2]
                Ey[i, j, 1] += PA[2] * Ey[i - 1, j, 1] + Ey[i - 1, j, 2]
                Ez[i, j, 1] += PA[3] * Ez[i - 1, j, 1] + Ez[i - 1, j, 2]
            elseif j > 1
                Ex[i, j, 1] += PB[1] * Ex[i, j - 1, 1] + Ex[i, j - 1, 2]
                Ey[i, j, 1] += PB[2] * Ey[i, j - 1, 1] + Ey[i, j - 1, 2]
                Ez[i, j, 1] += PB[3] * Ez[i, j - 1, 1] + Ez[i, j - 1, 2]
            end
            for t = 2 : i + j - 1
                if i > 1
                    # 9.5.6
                    Ex[i, j, t] += PA[1] * Ex[i - 1, j, t] + t * Ex[i - 1, j, t + 1] + oo2p * Ex[i - 1, j, t - 1]
                    Ey[i, j, t] += PA[2] * Ey[i - 1, j, t] + t * Ey[i - 1, j, t + 1] + oo2p * Ey[i - 1, j, t - 1]
                    Ez[i, j, t] += PA[3] * Ez[i - 1, j, t] + t * Ez[i - 1, j, t + 1] + oo2p * Ez[i - 1, j, t - 1]
                elseif j > 1
                    # 9.5.7
                    Ex[i, j, t] += PB[1] * Ex[i, j - 1, t] + t * Ex[i, j - 1, t + 1] + oo2p * Ex[i, j - 1, t - 1]
                    Ey[i, j, t] += PB[2] * Ey[i, j - 1, t] + t * Ey[i, j - 1, t + 1] + oo2p * Ey[i, j - 1, t - 1]
                    Ez[i, j, t] += PB[3] * Ez[i, j - 1, t] + t * Ez[i, j - 1, t + 1] + oo2p * Ez[i, j - 1, t - 1]
                end
            end
        end
    end
end

@inline function Fn(nmax::Int, T, ⍵::Float64, α::Float64)
    """ Computes the Boys function from the Kummer confluent
        hypergeometric function, eq. 9.8.39, modified per eq. 63 of
        dx.doi.org/10.1063/1.2956507 to make it an erf integral"""
    Fnvals = zeros(typeof(T), nmax + 1)
    δ = ⍵^2 / (⍵^2 + α)
    TT = T * δ

    Fnmax = if T < 1e-20    # T ≈ 0
        1.0 / (2 * nmax + 1)
    else 
        gamma(nmax + 0.5) * gamma_inc(nmax + 0.5, TT)[1] / (2 * TT ^ (nmax + 0.5))
    end

    Fnvals[nmax + 1] = Fnmax
    for n in nmax-1 : -1 : 0
        Fnvals[n + 1] = (2 * TT * Fnvals[n + 2] + exp(-TT)) / (2 * n + 1)
    end

    fac = sqrt(δ)
    for n in 0 : nmax
        Fnvals[n + 1] *= fac
        fac *= δ
    end

    return  Fnvals
end

@inline function generate_R_matrix!(R, maxam::Int, p::Float64, P, C,
                           α::Float64=1.0, β::Float64=0.0, ⍵::Float64=0.0)
    """ Makes the Hermite Coulomb integrals from derivatives of the Boys function, by recursion,
        defined as eqn. 9.9.13.  Generalized to deliver integrals of the form
        α + β(erf ⍵ R)
        --------------
              R       """
    PC = P .- C
    RPC = norm(PC)
    T = p * RPC ^ 2
    R .= 0

    # 1/R like integrals, which are just the erf in the limit ⍵->∞
    if α != 0.0
        Fnvals = Fn(maxam, T, 1e50, p)
        for n = 0 : maxam
            R[n+1, 1] += α * (-2 * p) ^ n * Fnvals[n + 1]
        end
    end

    # erf(⍵ R)/R integrals
    if β != 0.0
        Fnvals = Fn(maxam, T, ⍵, p)
        for n = 0 : maxam
            R[n+1, 1] += β * (-2 * p) ^ n * Fnvals[n + 1]
        end
    end

    for v = 0 : maxam
        for u = 0 : maxam - v
            for t = 0 : maxam - v - u
                for n = 1 : maxam
                    val = 0.0
                    if t == u == v == 0
                        # The Boys function has already been populated up above
                        continue
                    elseif t == u == 0
                        if v > 1
                            val += (v - 1) * R[n + 1, addr(t, u, v - 2, maxam)]
                        end
                        val += PC[3] * R[n + 1, addr(t, u, v - 1, maxam)]
                    elseif t == 0
                        if u > 1
                            val += (u - 1) * R[n + 1, addr(t, u - 2, v, maxam)]
                        end
                        val += PC[2] * R[n + 1, addr(t, u - 1, v, maxam)]
                    else
                        if t > 1
                            val += (t - 1) * R[n + 1, addr(t - 2, u, v, maxam)]
                        end
                        val += PC[1] * R[n + 1, addr(t - 1, u, v, maxam)]
                    end
                    R[n, addr(t, u, v, maxam)] = val
                end
            end
        end
    end
    R
end

end #module