import json
import numpy as np
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
    Convert geometry string to aos_*.json filename.
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
            tag += str(int(round(z * 10)))

        tags.append(tag)

    name = "_".join(tags)
    return f"/home/marco/hummr_tests/aos_{name}.json"


def _load_tensor(obj):
    """
    Reconstruct a numpy array from the { "dims": [...], "data": [...] } dict
    that export_tensors.cpp writes for every tensor.
    """
    dims = obj["dims"]
    data = np.array(obj["data"], dtype=np.float64)
    return data.reshape(dims, order="C")


def _load_arma_mat(obj):
    """
    Armadillo stores matrices column-major, so the flat buffer in "data" is
    col-major.  We reshape as (n_cols, n_rows) in C order and transpose to get
    the correct (n_rows, n_cols) matrix – identical to what h5py + numpy gave.
    """
    n_rows, n_cols = obj["dims"]
    data = np.array(obj["data"], dtype=np.float64)
    # col-major: reshape as (n_cols, n_rows) then transpose
    return data.reshape((n_cols, n_rows), order="C").T


def AO_ints(geometry, basis, max_mem=1e9, NucPotentialOnly=False, printout=print):

    start_time = timeit.default_timer()

    filename = _geom_to_filename(geometry)

    printout(f"Loading AO integrals from {filename}", flush=True)

    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    with open(filename, "r") as f:
        root = json.load(f)

    g = root["AOs"]

    if not NucPotentialOnly:
        S = _load_arma_mat(g["S"])
        T = _load_arma_mat(g["T"])
        V = brabraketket(_load_tensor(g["ERI"]))
        print("apply 1/2 definition on ERI in json load")
        V *= 0.5
    else:
        S = None
        T = None
        V = None

    U = _load_arma_mat(g["U"])
    print("U shape ", U.shape)

    elapsed = timeit.default_timer() - start_time
    printout("elapsed time (JSON integrals load) =", elapsed)

    return S, T, U, V, supplementary(filename)
