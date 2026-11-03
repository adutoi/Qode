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
import numpy
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


def AO_ints(geometry, basis, max_mem=1e9, NucPotentialOnly=False, printout=print, with_eri=True):
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
    elif not with_eri:
        # one-electron integrals only (the two-electron integrals come from RI_factors)
        S, T = pli.overlap_and_kinetic(atoms, basis)
        U = pli.nuclear_attraction(atoms, basis)
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


def lible_basis_name(name, aux=False):
    """
    The lible name of a basis set as hummr resolves it (hummr::ints::translateBasis2Lible): the lower-case name if
    lible knows it, otherwise hummr's correspondence table (e.g. def2-JK -> def2-universal-jkfit); None if neither.
    """
    known = set(pli.available_basis_sets()) | set(pli.available_aux_basis_sets())
    if name.lower() in known:
        return name.lower()
    table = {"DEF2-SV(P)": "def2-sv(p)", "DEF2-SVPD": "def2-svpd", "DEF2-TZVPP": "def2-tzvpp",
             "DEF2-TZVPD": "def2-tzvpd", "DEF2-TZVPPD": "def2-tzvppd", "DEF2-QZVPP": "def2-qzvpp",
             "DEF2-QZVPD": "def2-qzvppd", "DEF2-QZVPPD": "def2-qzvppd", "DEF2-JK": "def2-universal-jkfit",
             "DEF2-SVP-C": "def2-svp-rifit", "DEF2-SVPD-C": "def2-svpd-rifit", "DEF2-TZVP-C": "def2-tzvp-rifit",
             "DEF2-TZVPD-C": "def2-tzvpd-rifit", "DEF2-QZVP-C": "def2-qzvp-rifit", "DEF2-QZVPP-C": "def2-qzvpp-rifit",
             "DEF2-QZVPPD-C": "def2-qzvppd-rifit", "SV": "ahlrichs-vdz", "SVP": "ahlrichs-pvdz"}
    return table.get(name.upper())


def RI_factors(geometry, basis, aux_basis, C, max_batch_bytes=5e8, metric_thresh=1e-10, printout=print):
    """
    The RI (density-fitting) factors of the two-electron integrals of the orbitals C (nao x n, AO order of the
    geometry string, as in AO_ints),
        B[P, p, q] = sum_Q (pq|Q) [J^-1/2]_QP,        (pq|rs) ~ sum_P B[P, p, q] B[P, r, s],
    with the Coulomb metric J_PQ = (P|Q) of the auxiliary basis aux_basis (eigenvalues below metric_thresh are
    dropped, i.e. P runs over the kept eigenvectors).  The three-index AO integrals (mu nu|Q) are computed by lible in
    batches of auxiliary shells of at most about max_batch_bytes and transformed to the orbitals right away, so that the
    AO three-index tensor (nao^2 naux) is never held; the memory is that of B (naux n^2) twice.
    Ghost atoms ("@") carry basis and auxiliary functions, but no charge (irrelevant here).
    """
    start_time = timeit.default_timer()
    atoms = _parse_geometry(geometry)
    engine = pli.RIEngine(atoms, basis, aux_basis)
    nao, naux, n = engine.num_aos(), engine.num_aux(), C.shape[1]
    if C.shape[0] != nao:
        raise ValueError(f"RI_factors: {C.shape[0]} AO rows of C, but {nao} AOs of {basis} for the geometry")
    C = numpy.ascontiguousarray(C, dtype=float)
    L = numpy.empty((naux, n, n))
    n_shells, first = engine.num_aux_shells(), 0
    while first < n_shells:
        last, funcs = first, []
        while last < n_shells and (not funcs or 8.0 * (len(funcs) + len(engine.aux_shell_functions(last))) * nao * nao <= max_batch_bytes):
            funcs += engine.aux_shell_functions(last)
            last += 1
        batch = engine.eri3_batch(first, last)                       # (n_functions, nao, nao)
        L[funcs] = numpy.matmul(numpy.matmul(C.T, batch), C)         # (pq|Q) of the batch
        batch = None
        first = last
    w, v = numpy.linalg.eigh(engine.eri2())
    keep = w > metric_thresh
    Jm = v[:, keep] / numpy.sqrt(w[keep])                            # J^-1/2 restricted to the kept eigenvectors
    B = (Jm.T @ L.reshape(naux, n * n)).reshape(-1, n, n)
    printout(f"RI factors ({basis} / {aux_basis}): {nao} AOs, {naux} auxiliary functions ({keep.sum()} kept), "
             f"{n} orbitals, elapsed time {timeit.default_timer() - start_time:.2f} s", flush=True)
    return B
