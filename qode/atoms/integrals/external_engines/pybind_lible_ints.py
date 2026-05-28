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

import timeit
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pylible_ints as pli

from .psi4_ints import brabraketket

# Lible is the integral backend used by hummr

class supplementary(object):
    def __init__(self, sourcefile):
        self.sourcefile = sourcefile


def _parse_geometry(geometry):
    """
    Parse the geometry string used by the rest of the XR code.

    Expected format (one atom per line):
        symbol  x  y  z
    Ghost atoms are recognised by a leading "@" on the symbol
    (e.g. "@O", "@H"), which is the same convention used in the
    original hdf5_ints / json_ints code.

    Returns a list of (symbol, x, y, z, is_ghost) tuples.
    Coordinates are in Angstrom.
    """
    atoms = []
    for line in geometry.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split()
        symbol = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        is_ghost = symbol.startswith("@")
        atoms.append((symbol, x, y, z, is_ghost))
    return atoms


def _geom_to_key(geometry):
    """
    Build a short human-readable key for the geometry (used in diagnostic
    messages only, mirrors the old filename logic).
    """
    atoms = _parse_geometry(geometry)
    tags = []
    for symbol, x, y, z, is_ghost in atoms:
        element = symbol.replace("@", "").lower()
        tag = ("gh" if is_ghost else "") + element
        if abs(z) > 1e-10:
            tag += str(int(round(z)))
        tags.append(tag)
    return "_".join(tags)


def AO_ints(geometry, basis, max_mem=1e9, NucPotentialOnly=False, printout=print):
    """
    Compute AO integrals by calling lible directly via the pylible_ints
    pybind11 module.  Ghost atoms (symbol prefixed with "@") contribute
    basis functions but zero nuclear charge.

    Parameters
    ----------
    geometry : str
        Multi-line string with one atom per line: "symbol  x  y  z"
        Coordinates in Angstrom.  Ghost atoms: prefix symbol with "@".
    basis : str
        Basis set name understood by lible (e.g. "6-31g", "def2-svp").
    max_mem : float
        Ignored (kept for API compatibility with the HDF5/JSON path).
    NucPotentialOnly : bool
        If True, only the nuclear attraction matrix U is computed;
        S, T, and ERI are returned as None.
    printout : callable
        Logging function (default: print).

    Returns
    -------
    S, T, U, V, supplementary
        S, T : numpy arrays (nao x nao) or None
        U    : numpy array  (nao x nao)
        V    : numpy array  (nao x nao x nao x nao) or None
               Note: the 0.5 prefactor from the original ERI definition
               is applied here (matching the HDF5/JSON path convention).
        supplementary : object carrying a .sourcefile attribute (the
               geometry key string, for traceability).
    """
    start_time = timeit.default_timer()

    atoms = _parse_geometry(geometry)
    key   = _geom_to_key(geometry)

    printout(f"Computing AO integrals via pylible_ints for geometry: {key}", flush=True)

    if NucPotentialOnly:
        U = pli.nuclear_attraction(atoms, basis)
        S = None
        T = None
        V = None
    else:
        S, T, U, ERI = pli.ao_integrals(atoms, basis)

        V = brabraketket(ERI)
        print("apply 1/2 definition on ERI in pylible_ints load")
        V = V * 0.5

    print("U shape ", U.shape)

    elapsed = timeit.default_timer() - start_time
    printout("elapsed time (pylible_ints AO integrals) =", elapsed)

    return S, T, U, V, supplementary(key)
