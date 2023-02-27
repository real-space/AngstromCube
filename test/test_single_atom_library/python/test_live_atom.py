#! /usr/bin/env python3.11
######################################################
### interface the LiveAtom library
######################################################
import ctypes as ct
import pathlib
import numpy as np

c_int_p = ct.POINTER(ct.c_int)

if __name__ == "__main__":
    libname = pathlib.Path().absolute() / "libliveatom.so"
    print("# Try to open dynamic library", libname)
    LA = ct.CDLL(libname)
    status = np.zeros(1, np.int32)
    status = status.astype(np.int32)
    stat = status.ctypes.data_as(c_int_p)
    LA.live_atom_init_env_(ct.c_char_p(b"control.sh"), stat)
    print("# init_env returns", status[0])
    LA.live_atom_set_env_(ct.c_char_p(b"called.from"), ct.c_char_p(b"python"), stat)
    print("# set_env returns", status[0])

    # the total number of atoms
    natoms = 1
    ar2 = 16.
    nr2 = 2**12
