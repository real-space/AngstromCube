#! /usr/bin/env python3.11
######################################################
### Python interface for the LiveAtom library
######################################################
import ctypes as ct
import pathlib
import numpy as np

c_int_p = ct.POINTER(ct.c_int32)
Float64 = ct.c_double
c_double_p = ct.POINTER(Float64)
Int32 = ct.c_int32


if __name__ == "__main__":
    libname = pathlib.Path().absolute() / "libliveatom.so"
    print("# Try to open dynamic library", libname)
    LA = ct.CDLL(libname)
    status = np.zeros(1, np.int32)
    status = status.astype(np.int32)
    stat = status.ctypes.data_as(c_int_p)

    LA.live_atom_init_env_(ct.c_char_p(b"control.sh"), stat)
    print("### init_env returns", status[0])

    LA.live_atom_set_env_(ct.c_char_p(b"called.from"), ct.c_char_p(b"python"), stat)
    print("### set_env returns", status[0])

    # the total number of atoms
    natoms = 2
    ar2 = 16.
    nr2 = 2**12
    na = Int32(natoms)

    # pointers = ???

    Z_core = np.ones(natoms, dtype=Float64)
    atom_id = np.zeros(natoms, dtype=Int32)
    numax = np.ones(natoms, dtype=Int32)
    sigma = np.ones(natoms, dtype=Float64)
    rcut = np.ones(natoms, dtype=Float64)
    nn = np.zeros((natoms, 8), dtype=ct.c_byte) ### is this the right array order?
    ionization = np.zeros(natoms, dtype=Float64)
    magnetization = np.zeros(natoms, dtype=Float64)
    xc_key = b"LDA" # .encode('utf-8')
    stride = Int32(0)
    lmax_qlm = np.ones(natoms, dtype=Int32)
    lmax_vlm = np.ones(natoms, dtype=Int32)
    n_valence_electrons = np.zeros(natoms, dtype=Float64)

    # numax .*= -9 # if we want to signal to load pawxml files, HOW IS THAT DONE IN PYTHON?

    for ia in range(natoms):
        Z_core[ia] = ia + 1
        atom_id[ia] = ia
        # pointers[ia] = np.zeros(nr2, dtype=Float64)

    LA.live_atom_initialize_(ct.byref(na)
                             , Z_core.ctypes.data_as(c_double_p)
                             , atom_id.ctypes.data_as(c_int_p)
                             , numax.ctypes.data_as(c_int_p)
                             , sigma.ctypes.data_as(c_double_p)
                             , rcut.ctypes.data_as(c_double_p)
                             , nn.ctypes.data_as(ct.POINTER(ct.c_byte))
                             , ionization.ctypes.data_as(c_double_p)
                             , magnetization.ctypes.data_as(c_double_p)
                             , ct.c_char_p(xc_key)
                             , ct.byref(stride)
                             , lmax_qlm.ctypes.data_as(c_int_p)
                             , lmax_vlm.ctypes.data_as(c_int_p)
                             , n_valence_electrons.ctypes.data_as(c_double_p)
                             , stat)
    print("### live_atom_initialize_ returns", status[0])

    # show number of partial wave per ell-channel of each atom
    for ia in range(natoms):
        print("### atom#"+str(ia),"\thas nn=", nn[ia])
    else:
        print("###")
        print("")

    LA.live_atom_finalize_(ct.byref(na), stat)
    print("### live_atom_finalize_ returns", status[0])

