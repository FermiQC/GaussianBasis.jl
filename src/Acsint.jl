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

    E = similar(A, 3, am1 + am2 + 2, am1 + 1, am2 + 1)
    E .= 0.0

    for (ca, a) in zip(B1.coef, B1.exp)
        for (cb, b) in zip(B2.coef, B2.exp)
            p = a + b
            P = (a * A + b * B) / p
            prefac = ca * cb * (π / p) ^ 1.5

            generate_E_matrix!(E, am1, am2, P, A, B, a, b)

            index1 = 1
            for l1 in am1 : -1 : 0
                for n1 in 0 : am1 - l1
                    m1 = am1 - l1 - n1
                    index2 = 1
                    for l2 in am2 : -1 : 0
                        for n2 in 0 : am2 - l2
                            m2 = am2 - l2 - n2
                            out[index1 + N1*(index2-1)] += E[1, 1, l1+1, l2+1] * E[2, 1, m1+1, m2+1] * E[3, 1, n1+1, n2+1] * prefac
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

    E = similar(A, 3, am1 + am2 + 4, am1 + 3, am2 + 1)
    for (ca, a) in zip(B1.coef, B1.exp)
        for (cb, b) in zip(B2.coef, B2.exp)
            p = a + b
            P = (a * A + b * B) / p
            prefac = -0.5 * ca * cb * (π / p) ^ 1.5

            generate_E_matrix!(E, am1+2, am2, P, A, B, a, b)

            index1 = 1
            for l1 in am1 : -1 : 0
                for n1 in 0 : am1 - l1
                    m1 = am1 - l1 - n1
                    index2 = 1
                    for l2 in am2 : -1 : 0
                        for n2 in 0 : am2 - l2
                            m2 = am2 - l2 - n2
                            Sx0 = E[1, 1, l1+1, l2+1]
                            Sy0 = E[2, 1, m1+1, m2+1]
                            Sz0 = E[3, 1, n1+1, n2+1]
                            # 9.3.31
                            Dx2 = 2 * a * (2 * a * E[1, 1, l1+3, l2+1] - (2 * l1 + 1) * E[1, 1, l1+1, l2+1])
                            if l1 > 1
                                Dx2 += l1 * (l1 - 1) * E[1, 1, l1-1, l2+1]
                            end
                            Dy2 = 2 * a * (2 * a * E[2, 1, m1+3, m2+1] - (2 * m1 + 1) * E[2, 1, m1+1, m2+1])
                            if m1 > 1
                                Dy2 += m1 * (m1 - 1) * E[2, 1, m1-1, m2+1]
                            end
                            Dz2 = 2 * a * (2 * a * E[3, 1, n1+3, n2+1] - (2 * n1 + 1) * E[3,1,n1+1, n2+1])
                            if n1 > 1
                                Dz2 += n1 * (n1 - 1) * E[3, 1, n1-1, n2+1]
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

struct FnBuffer{T}
    Γnmax::T
    nmax::Int
    vals::Vector{T}
end
FnBuffer{T}(nmax) where {T} = FnBuffer{T}(gamma(nmax+0.5), nmax, zeros(T, nmax+1))


@inline function Fn!(buffer::FnBuffer, nmax::Int, T, ⍵::Float64, α::Float64)
    """ Computes the Boys function from the Kummer confluent
        hypergeometric function, eq. 9.8.39, modified per eq. 63 of
        dx.doi.org/10.1063/1.2956507 to make it an erf integral"""
    δ = ⍵^2 / (⍵^2 + α)
    TT = T * δ

    Fnmax = if T < 1e-20    # T ≈ 0
        1.0 / (2 * nmax + 1)
    else 
        buffer.Γnmax * gamma_inc(nmax + 0.5, TT)[1] / (2 * TT ^ (nmax + 0.5))
    end

    buffer.vals[nmax + 1] = Fnmax
    for n in nmax-1 : -1 : 0
        @inbounds buffer.vals[n + 1] = (2 * TT * buffer.vals[n + 2] + exp(-TT)) / (2 * n + 1)
    end

    fac = sqrt(δ)
    for n in 0 : nmax
        @inbounds buffer.vals[n + 1] *= fac
        fac *= δ
    end

    return buffer.vals
end

@inline function generate_R_matrix!(R, fnbuf::FnBuffer, p::Float64, P, C,
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

    maxam = fnbuf.nmax

    # 1/R like integrals, which are just the erf in the limit ⍵->∞
    if α != 0.0
        Fn!(fnbuf, maxam, T, 1e50, p)
        for n = 0 : maxam
            @inbounds R[1,1,1, n + 1] += α * (-2 * p) ^ n * fnbuf.vals[n + 1]
        end
    end

    # erf(⍵ R)/R integrals
    if β != 0.0
        Fn!(fnbuf, maxam, T, ⍵, p)
        for n = 0 : maxam
            @inbounds R[1,1,1, n+1] += β * (-2 * p) ^ n * fnbuf.vals[n + 1]
        end
    end

    @inbounds for v = 0 : maxam
        for u = 0 : maxam - v
            for t = 0 : maxam - v - u
                for n = 1 : maxam
                    val = 0.0
                    if t == u == v == 0
                        # The Boys function has already been populated up above
                        continue
                    elseif t == u == 0
                        if v > 1
                            val += (v - 1) * R[t+1, u+1, v - 1, n + 1]
                        end
                        val += PC[3] * R[t+1, u+1, v, n + 1]
                    elseif t == 0
                        if u > 1
                            val += (u - 1) * R[t+1, u - 1, v+1, n + 1]
                        end
                        val += PC[2] * R[t+1, u, v+1, n + 1]
                    else
                        if t > 1
                            val += (t - 1) * R[t - 1, u+1, v+1, n + 1]
                        end
                        val += PC[1] * R[t, u+1, v+1, n + 1]
                    end
                    R[t+1, u+1, v+1, n] = val
                end
            end
        end
    end
    R
end

function generate_V_pair!(out, BS::BasisSet, s1, s2)
    generate_V_pair!(out, BS.basis[s1], BS.basis[s2], BS.atoms)
end

function generate_V_pair!(out, B1::CartesianShell, B2::CartesianShell, atoms::Vector{<:Atom})

    out .= 0.0

    am1 = B1.l
    am2 = B2.l
    am = am1 + am2
    fnbuf = FnBuffer{eltype(out)}(am)

    A = B1.atom.xyz .* ang2bohr
    B = B2.atom.xyz .* ang2bohr

    N1 = num_basis(B1)

    Rmat = zeros(eltype(A), am + 1, am + 1, am + 1, am + 1)
    ints = zeros(eltype(A), cart_dim(am1), cart_dim(am2))

    E = similar(A, 3, am1 + am2 + 2, am1 + 1, am2 + 1)
    E .= 0.0
    for At in 1:length(atoms)
        Z = atoms[At].Z
        C = atoms[At].xyz .* ang2bohr
        for (ca, a) in zip(B1.coef, B1.exp)
            for (cb, b) in zip(B2.coef, B2.exp)
                p = a + b
                P = (a * A + b * B) / p
                generate_E_matrix!(E, am1, am2, P, A, B, a, b)
                prefac = -Z * ca * cb * 2 * π / p

                generate_R_matrix!(Rmat, fnbuf, p, P, C)
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
                                            Eab = E[1, t+1, l1+1, l2+1] * E[2, u+1, m1+1, m2+1] * E[3, v+1, n1+1, n2+1]
                                            # eqn. 9.9.32
                                            val += Eab * Rmat[t+1, u+1, v+1, 1]
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
    
function generate_quanta_list2(am2, am1, am, N1)
    dim = (am1 + 1)*(am1 + 2)*(am2 + 1)*(am2 + 2)*(am1+am2+3)*(am1^2+am2^2+4*am1*am2+9*(am1+am2)+20)÷240
    vals = Vector{NTuple{6,Int}}(undef, dim)
    L1 = am1+am2+2
    L12 = L1*(am1+1)
    idx2 = 0
    i = 1
    for l2 in am2:-1:0
        for n2 in 0:am2 - l2
            m2 = am2 - l2 - n2
            idx1 = 0
            for l1 in am1:-1:0
                for n1 in 0:am1 - l1
                    m1 = am1 - l1 - n1
                    for v = 0:n1 + n2, u = 0:m1 + m2, t = 0:l1 + l2
                        E1idx = 1 + 3*(t + L1*l1 + L12*l2)
                        E2idx = 2 + 3*(u + L1*m1 + L12*m2)
                        E3idx = 3 + 3*(v + L1*n1 + L12*n2)
                        Ridx = t + (am+1)*(u + (am+1)*v)
                        s = 1-2*isodd(t+u+v)
                        vals[i] = (idx2*N1 + idx1, E1idx, E2idx, E3idx, Ridx, s)
                        i += 1
                    end
                    idx1 += 1
                end
            end
            idx2 += 1
        end
    end
    return vals
end

function generate_ERI_quartet!(out, BS::BasisSet, s1, s2, s3, s4, α::Float64=1.0, β::Float64=0.0, ω::Float64=0.0)
    generate_ERI_quartet!(out, BS.basis[s1], BS.basis[s2], BS.basis[s3], BS.basis[s4], α, β, ω)
end

function generate_ERI_quartet!(out, B1::CartesianShell, B2::CartesianShell, B3::CartesianShell, B4::CartesianShell,
                              α::Float64=1.0, β::Float64=0.0, ⍵::Float64=0.0; cutoff = 1e-16)

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
    fnbuf = FnBuffer{T}(am)

    N1   = num_basis(B1)
    N12  = N1*num_basis(B2)


    A = atom1.xyz .* ang2bohr
    B = atom2.xyz .* ang2bohr
    C = atom3.xyz .* ang2bohr
    D = atom4.xyz .* ang2bohr

    Rmat = zeros(T, am + 1,am + 1,am + 1, am + 1)

    bra_E = zeros(T, 3, am1 + am2 + 2, am1 + 1, am2 + 1)
    ket_E = zeros(T, 3, am3 + am4 + 2, am3 + 1, am4 + 1)

    quanta1 = generate_quanta_list2(am2,am1,am, num_basis(B1))
    quanta2 = generate_quanta_list2(am4,am3,am, num_basis(B3))

    for (ca, a) in zip(B1.coef, B1.exp)
        for (cb, b) in zip(B2.coef, B2.exp)
            p = a + b
            P = (a * A + b * B) / p
            generate_E_matrix!(bra_E, am1, am2, P, A, B, a, b)
            @inbounds for (cc, c) in zip(B3.coef, B3.exp)
                for (cd, d) in zip(B4.coef, B4.exp)
                    q = c + d
                    Q = (c * C + d * D) / q

                    θ = p * q / (p + q)
                    prefac = 2 * π^(5/2) * ca * cb * cc * cd / (p * q * sqrt(p + q))
                    if abs(prefac) < cutoff
                        continue
                    end

                    generate_E_matrix!(ket_E, am3, am4, Q, C, D, c, d)
                    generate_R_matrix!(Rmat, fnbuf, θ, P, Q, α, β, ⍵)
                    for (idx43, ikE1, ikE2, ikE3, Ridx43, s) in quanta2
                        Ecd = prefac * s * ket_E[ikE1] * ket_E[ikE2] * ket_E[ikE3]
                        if abs(Ecd) < cutoff
                            continue
                        end
                        midx43 = idx43 * N12
                        for (idx21, ibE1, ibE2, ibE3, Ridx21, _) in quanta1
                            Eab = bra_E[ibE1] * bra_E[ibE2] * bra_E[ibE3]
                            idx = midx43 + idx21 + 1
                            # eqn. 9.9.33
                            out[idx] += Ecd * Eab * Rmat[Ridx21 + Ridx43 + 1]
                        end
                    end
                end
            end
        end
    end
    return out
end

@inline function generate_E_matrix!(E, maxam1, maxam2,
                           P, A, B,
                           a::Float64, b::Float64)
    # Makes the Hermite->Cartesian conversion factors by recursion
    E .= 0

    # 9.2.11
    p = a + b
    # 9.2.14
    AB = A - B
    PA = P - A
    PB = P - B
    # 9.2.12
    µ = a * b / p
    # 9.2.15
    @inbounds @fastmath @. E[:, 1, 1, 1] = exp(-µ * AB^2)
    oo2p = 1 / (2 * p)
    @inbounds for j = 1 : maxam2 + 1
        for i = 1 : maxam1 + 1
            # handle t = 1 case
            if i > 1
                @views E[:, 1, i, j] .+= PA .* E[:, 1, i - 1, j] .+ E[:, 2, i - 1, j]
            elseif j > 1
                @views E[:, 1, i, j] .+= PB .* E[:, 1, i, j - 1] .+ E[:, 2, i, j - 1]
            end
            for t = 2 : i + j - 1
                if i > 1
                    # 9.5.6
                    @. @views E[:, t, i, j] += PA * E[:, t, i - 1, j] + t * E[:, t + 1, i - 1, j] + oo2p * E[:, t - 1, i - 1, j]
                elseif j > 1
                    # 9.5.7
                    @. @views E[:, t, i, j] += PB * E[:, t, i, j - 1] + t * E[:, t + 1, i, j - 1] + oo2p * E[:, t - 1, i, j - 1]
                end
            end
        end
    end
end

end #module
