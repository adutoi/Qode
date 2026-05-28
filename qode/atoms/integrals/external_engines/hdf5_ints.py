#    (C) Copyright 2026 Marco Bauer
# 
#    This file is part of Qode.
# 
#    Qode is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    Qode is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with Qode.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np
import h5py
import timeit
import os

from .psi4_ints import brabraketket

class supplementary(object):
    def __init__(self, sourcefile):
        self.sourcefile = sourcefile


def _parse_geometry(geometry):
    atoms = []
    for line in geometry.strip().split("\n"):
        if not line.strip():
            continue
        el, x, y, z = line.split()
        atoms.append((el, float(z)))
    return atoms


def _geom_to_filename(geometry):
    """
    Convert geometry string to aos_*.h5 filename.
    """
    atoms = _parse_geometry(geometry)

    tags = []
    for el, z in atoms:
        ghost = el.startswith("@")
        element = el.replace("@", "").lower()

        tag = ""
        if ghost:
            tag += "gh"

        tag += element

        if abs(z) > 1e-10:
            tag += str(int(round(z * 10)))  # dirty hack for simpler filenames

        tags.append(tag)

    name = "_".join(tags)
    return f"/home/marco/hummr_tests/aos_{name}.h5"


def AO_ints(geometry, basis, max_mem=1e9, NucPotentialOnly=False, printout=print):

    start_time = timeit.default_timer()

    filename = _geom_to_filename(geometry)

    printout(f"Loading AO integrals from {filename}", flush=True)

    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    with h5py.File(filename, "r") as f:

        g = f["AOs"]

        if not NucPotentialOnly:
            S = np.array(g["S"])
            T = np.array(g["T"])
            V = brabraketket(np.array(g["ERI"]))
            #print("currently no brabraketket transformation is applied in the hdf5 import!!!!!!!!!!!!!!!!!!!!!!!")
            #V = np.array(g["ERI"])
            print("apply 1/2 definition on ERI in hdf5 load")
            V *= 0.5
        else:
            S = None
            T = None
            V = None

        U = np.array(g["U"])
        print("U shape ", U.shape)

    elapsed = timeit.default_timer() - start_time
    printout("elapsed time (HDF5 integrals load) =", elapsed)

    return S, T, U, V, supplementary(filename)
