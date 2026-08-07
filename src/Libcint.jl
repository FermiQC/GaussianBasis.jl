"""
    Libcint

Minimal wrap around the integral library libcint. This module exposes
libcint functions to the Julia interface. 
"""
module Libcint

using GaussianBasis: BasisSet, LCint

export cint1e_kin_sph!, cint1e_nuc_sph!, cint1e_ovlp_sph!, cint2c2e_sph!, cint2e_sph!, cint3c2e_sph!
export cint1e_ipkin_sph!, cint1e_ipnuc_sph!, cint1e_ipovlp_sph!, cint1e_ipipovlp_sph!, cint2e_ip1_sph!, cint1e_r_sph!
export cint1e_r_sph!, cint1e_rr_sph!, cint1e_rrr_sph!, cint1e_rrrr_sph!
export cint3c2e_ip2_sph!, cint3c2e_ip1_sph!, cint2c2e_ip1_sph!
export cint1e_ipipkin_sph!, cint1e_ipipnuc_sph!, cint1e_ipovlpip_sph!, cint1e_ipkinip_sph!, cint1e_ipnucip_sph!
export cint2e_ipip1_sph!, cint2e_ip1ip2_sph!
export cint1e_rinv_sph!, cint1e_iprinv_sph!, cint1e_ipiprinv_sph!, cint1e_iprinvip_sph!
export cint2e_ipvip1_sph!

using libcint_jll
const LIBCINT = libcint

function CINTcgtos_spheric(bas_id, bas)
    @ccall LIBCINT.CINTcgtos_spheric(bas_id::Cint, bas::Ptr{Cint})::Cint
end

function cint1e_ovlp_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ovlp_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                )::Cvoid
end
function cint1e_ovlp_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ovlp_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_kin_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_kin_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                )::Cvoid
end
function cint1e_kin_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_kin_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_nuc_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_nuc_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                )::Cvoid
end
function cint1e_nuc_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_nuc_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint2e_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint2e_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                )::Cvoid
end
function cint2e_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint2e_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint2c2e_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint2c2e_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                )::Cvoid
end
function cint2c2e_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint2c2e_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint3c2e_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint3c2e_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                )::Cvoid
end
function cint3c2e_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint3c2e_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_ipovlp_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipovlp_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipovlp_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipovlp_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_ipipovlp_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipipovlp_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipipovlp_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipipovlp_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end


function cint1e_ipkin_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipkin_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipkin_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipkin_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_ipnuc_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipnuc_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipnuc_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipnuc_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

# --- Second-derivative one-electron kernels (Hessians) ---
# "ipip" = same-center double derivative, "Xip Y ip" (e.g. ipovlpip) = one
# derivative on each of the two shell centers (cross-center).

function cint1e_ipipkin_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipipkin_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipipkin_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipipkin_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_ipipnuc_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipipnuc_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipipnuc_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipipnuc_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_ipovlpip_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipovlpip_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipovlpip_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipovlpip_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_ipkinip_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipkinip_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipkinip_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipkinip_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_ipnucip_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipnucip_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipnucip_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipnucip_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

# --- Repositionable point-charge (rinv) kernels ---
# 1/|r-Rc| operator with Rc read from env[PTR_RINV_ORIG] (env[5:7], 1-indexed)
# instead of a fixed nuclear position. Caller is responsible for writing Rc
# into a *copy* of env before calling these (see NuclearHess.jl).

function cint1e_rinv_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_rinv_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_rinv_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_rinv_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_iprinv_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_iprinv_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_iprinv_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_iprinv_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_ipiprinv_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_ipiprinv_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_ipiprinv_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_ipiprinv_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_iprinvip_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_iprinvip_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_iprinvip_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_iprinvip_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint2e_ip1_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint2e_ip1_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                   )::Cvoid
end
function cint2e_ip1_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint2e_ip1_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

# --- Second-derivative two-electron kernels (Hessians) ---
# ipip1 = double derivative on the first shell center (same-center); ip1ip2 =
# one derivative each on the first and second shell centers (cross-center).
# Together with shell-argument permutation and translational invariance,
# these two are enough to build every atom-pair ERI Hessian block.

function cint2e_ipip1_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint2e_ipip1_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                   )::Cvoid
end
function cint2e_ipip1_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint2e_ipip1_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint2e_ip1ip2_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint2e_ip1ip2_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                   )::Cvoid
end
function cint2e_ip1ip2_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint2e_ip1ip2_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint2e_ipvip1_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint2e_ipvip1_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                   )::Cvoid
end
function cint2e_ipvip1_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint2e_ipvip1_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint2c2e_ip1_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint2c2e_ip1_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                )::Cvoid
end
function cint2c2e_ip1_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint2c2e_ip1_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint3c2e_ip1_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint3c2e_ip1_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                )::Cvoid
end
function cint3c2e_ip1_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint3c2e_ip1_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint3c2e_ip2_sph!(buf, shls, atm, natm, bas, nbas, env)
    opt = Ptr{UInt8}(C_NULL)
    @ccall LIBCINT.cint3c2e_ip2_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt :: Ptr{UInt8},
                                )::Cvoid
end
function cint3c2e_ip2_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint3c2e_ip2_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_r_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_r_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_r_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_r_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_rr_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_rr_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_rr_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_rr_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_rrr_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_rrr_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_rrr_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_rrr_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

function cint1e_rrrr_sph!(buf, shls, atm, natm, bas, nbas, env)
    @ccall LIBCINT.cint1e_rrrr_sph(
                                    buf  :: Ptr{Cdouble},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble}
                                   )::Cvoid
end
function cint1e_rrrr_sph!(buf::Array{Cdouble}, shls::AbstractArray{<:Integer}, lib::LCint)
    cint1e_rrrr_sph!(buf, Cint.(shls.-1), lib.atm, lib.natm, lib.bas, lib.nbas, lib.env)
end

end #module
