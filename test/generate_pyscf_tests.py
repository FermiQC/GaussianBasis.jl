import pyscf
import h5py

atom = """
C   -2.131551124300    2.286168823700    0.000000000000
H   -1.061551124300    2.286168823700    0.000000000000
H   -2.488213906200    1.408104616400    0.496683911300
H   -2.488218762100    2.295059432700   -1.008766153900
H   -2.488220057000    3.155340844300    0.512081313000
"""
basis = "sto-3g"

m = pyscf.M(atom=atom, basis=basis)

print(m.bas_exps())

with h5py.File("pyscf_test.h5", "w") as f:
    f["/atom"] = atom
    f["/basis"] = basis
    f["/kinetic"] = m.intor("cint1e_kin_sph")
    f["/dipole"] = m.intor("cint1e_r_sph")
    dim = f["/dipole"].shape[-1]
    f["/quadrupole"] = m.intor("cint1e_rr_sph").reshape(3, 3, dim, dim)
    f["/octupole"] = m.intor("cint1e_rrr_sph").reshape(3, 3, 3, dim, dim)
    f["/hexadecapole"] = m.intor("cint1e_rrrr_sph").reshape(3, 3, 3, 3, dim, dim)
