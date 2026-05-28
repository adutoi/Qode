/*   (C) Copyright 2026 Marco Bauer
 *
 *   This file is part of Qode.
 *
 *   Qode is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Qode is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Qode.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * pylible_ints.cpp
 *
 * Pybind11 module that exposes lible AO integrals directly to Python,
 * with full support for ghost atoms.
 *
 * Python API (see docstrings below):
 *
 *   import pylible_ints as pli
 *
 *   # Atom list: each entry is (symbol, x_ang, y_ang, z_ang, is_ghost)
 *   atoms = [
 *       ("O",  0.0,    0.0,    0.119748, False),
 *       ("H",  0.0,    0.756950, -0.478993, False),
 *       ("H",  0.0,   -0.756950, -0.478993, False),
 *       ("@O", 3.0,    0.0,    0.119748, True),   # ghost oxygen
 *   ]
 *
 *   S, T, U, ERI = pli.ao_integrals(atoms, "6-31g")
 *
 * All matrices/tensors are returned as numpy arrays.
 * S, T, U  : 2-D (nao x nao) arrays.
 * ERI      : 4-D (nao x nao x nao x nao) array, row-major (C order).
 *
 * Ghost atom convention
 * ---------------------
 * A ghost atom contributes basis functions at its position, but contributes
 * ZERO to the nuclear attraction integrals.  This is the standard behaviour
 * of "Bq" / "@X" atoms as used in counterpoise corrections.
 *
 * Implementation note: we construct the lible::ints::Structure via the
 * basis_atoms_t constructor.  For every atom (real or ghost) we look up
 * the basis functions by the real atomic number, but we store atomic_nr = 0
 * for ghost atoms.  Because Structure::getZ(iatom) returns 0 for ghosts, the
 * nuclearAttraction() helper will multiply their contribution by zero, which
 * is exactly the ghost-atom convention.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <lible/ints/defs.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/structure.hpp>
#include <lible/ints/shell.hpp>
#include <lible/types.hpp>

#include <array>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace py    = pybind11;
namespace lints = lible::ints;

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

/// Look up the atomic number for a symbol, stripping a leading "@" for ghost atoms.
static int resolveAtomicNr(const std::string &symbol)
{
    std::string s = symbol;
    if (!s.empty() && s[0] == '@')
        s = s.substr(1);

    auto it = lints::atomic_numbers.find(s);
    if (it == lints::atomic_numbers.end())
        throw std::runtime_error("pylible_ints: unknown element symbol '" + symbol + "'");

    return it->second;
}

/// Convert a lible::vec2d to a 2-D numpy array (row-major copy).
static py::array_t<double> vec2d_to_numpy(const lible::vec2d &m)
{
    size_t rows = m.dim<0>();
    size_t cols = m.dim<1>();

    py::array_t<double> arr({rows, cols});
    auto buf = arr.mutable_unchecked<2>();

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            buf(i, j) = m(i, j);

    return arr;
}

/// Convert a lible::vec4d to a 4-D numpy array (row-major copy).
static py::array_t<double> vec4d_to_numpy(const lible::vec4d &t)
{
    size_t d0 = t.dim<0>();
    size_t d1 = t.dim<1>();
    size_t d2 = t.dim<2>();
    size_t d3 = t.dim<3>();

    py::array_t<double> arr({d0, d1, d2, d3});
    auto buf = arr.mutable_unchecked<4>();

    for (size_t a = 0; a < d0; a++)
    for (size_t b = 0; b < d1; b++)
    for (size_t c = 0; c < d2; c++)
    for (size_t d = 0; d < d3; d++)
        buf(a, b, c, d) = t(a, b, c, d);

    return arr;
}

// ---------------------------------------------------------------------------
// Build a lible::ints::Structure that handles ghost atoms correctly.
//
// atoms: list of (symbol, x_ang, y_ang, z_ang, is_ghost)
//   symbol   – element symbol, optionally prefixed with "@" (ignored if
//              is_ghost is also True; either convention is accepted)
//   x/y/z    – coordinates in Angstrom
//   is_ghost – if True the atom contributes basis functions but Z = 0
// ---------------------------------------------------------------------------
static lints::Structure buildStructure(
    const std::vector<std::tuple<std::string, double, double, double, bool>> &atoms,
    const std::string &basis_set)
{
    std::vector<int>                       atomic_nrs;
    std::vector<std::array<double, 3>>     coords_ang;
    lints::basis_atoms_t                   basis_atoms;

    for (const auto &[symbol, x, y, z, is_ghost] : atoms)
    {
        int real_Z = resolveAtomicNr(symbol);
        // For the Structure's atomic number list: 0 for ghosts (nuclr. attraction = 0),
        // real Z for true atoms.
        int structure_Z = is_ghost ? 0 : real_Z;

        atomic_nrs.push_back(structure_Z);
        coords_ang.push_back({x, y, z});

        // Basis functions are looked up by the REAL atomic number regardless of ghost status.
        lints::BasisAtom ba = lints::basisForAtom(real_Z, basis_set);
        // Override the stored atomic number to match what the Structure will use.
        ba.atomic_nr_ = structure_Z;
        basis_atoms.push_back(ba);
    }

    // Use the basis_atoms_t constructor so we control exactly which basis shells
    // are placed on each atom (enabling different element bases on ghost positions).
    return lints::Structure(basis_atoms, atomic_nrs, coords_ang);
}

// ---------------------------------------------------------------------------
// Primary exposed function
// ---------------------------------------------------------------------------

/**
 * Compute S, T, U, ERI for the given molecular geometry and basis set.
 *
 * Parameters
 * ----------
 * atoms : list of (symbol, x, y, z, is_ghost)
 *     symbol  : element symbol (str), may be prefixed with "@" for ghosts.
 *     x, y, z : coordinates in Angstrom (float).
 *     is_ghost: bool – if True the nuclear charge is set to zero.
 * basis_set : str
 *     Basis set name as recognised by lible (e.g. "6-31g", "def2-svp").
 *
 * Returns
 * -------
 * S   : np.ndarray, shape (nao, nao)  – overlap matrix
 * T   : np.ndarray, shape (nao, nao)  – kinetic energy matrix
 * U   : np.ndarray, shape (nao, nao)  – nuclear attraction matrix (ghost contributions = 0)
 * ERI : np.ndarray, shape (nao, nao, nao, nao) – two-electron repulsion integrals, C order
 */
static py::tuple ao_integrals(
    const std::vector<std::tuple<std::string, double, double, double, bool>> &atoms,
    const std::string &basis_set)
{
    lints::Structure structure = buildStructure(atoms, basis_set);

    py::array_t<double> S   = vec2d_to_numpy(lints::overlap(structure));
    py::array_t<double> T   = vec2d_to_numpy(lints::kineticEnergy(structure));
    py::array_t<double> U   = vec2d_to_numpy(lints::nuclearAttraction(structure));
    py::array_t<double> ERI = vec4d_to_numpy(lints::eri4(structure));

    return py::make_tuple(S, T, U, ERI);
}

/**
 * Like ao_integrals(), but computes S and T only.
 * Useful when U and ERI are not needed (e.g. for overlap-only calculations).
 */
static py::tuple overlap_and_kinetic(
    const std::vector<std::tuple<std::string, double, double, double, bool>> &atoms,
    const std::string &basis_set)
{
    lints::Structure structure = buildStructure(atoms, basis_set);

    py::array_t<double> S = vec2d_to_numpy(lints::overlap(structure));
    py::array_t<double> T = vec2d_to_numpy(lints::kineticEnergy(structure));

    return py::make_tuple(S, T);
}

/**
 * Return only the nuclear attraction matrix U.
 * Ghost atoms (is_ghost=True or symbol prefixed with "@") contribute zero.
 */
static py::array_t<double> nuclear_attraction(
    const std::vector<std::tuple<std::string, double, double, double, bool>> &atoms,
    const std::string &basis_set)
{
    lints::Structure structure = buildStructure(atoms, basis_set);
    return vec2d_to_numpy(lints::nuclearAttraction(structure));
}

/**
 * Return only the ERI4 tensor.
 */
static py::array_t<double> eri4_integrals(
    const std::vector<std::tuple<std::string, double, double, double, bool>> &atoms,
    const std::string &basis_set)
{
    lints::Structure structure = buildStructure(atoms, basis_set);
    return vec4d_to_numpy(lints::eri4(structure));
}

/**
 * Return the number of (spherical) AOs for the given atom list and basis set.
 * Useful for pre-allocating arrays before calling the integral functions.
 */
static size_t num_aos(
    const std::vector<std::tuple<std::string, double, double, double, bool>> &atoms,
    const std::string &basis_set)
{
    lints::Structure structure = buildStructure(atoms, basis_set);
    return structure.getDimAO();
}

/**
 * Set the path to the main basis set directory at runtime.
 * Needed when lible was not compiled with LIBLE_MAIN_BASIS_DIR pointing at
 * the correct location (e.g. after installing into a non-standard prefix).
 */
static void set_basis_path(const std::string &path)
{
    lints::BasisPaths::setMainBasisSetsPath(path);
}

/**
 * Set the path to the auxiliary basis set directory at runtime.
 */
static void set_aux_basis_path(const std::string &path)
{
    lints::BasisPaths::setAuxBasisSetsPath(path);
}

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------

PYBIND11_MODULE(pylible_ints, m)
{
    m.doc() = R"doc(
        pylible_ints
        ============
        Direct Python bindings to the lible integral library.

        Ghost atoms are supported: set ``is_ghost=True`` (or prefix the symbol
        with "@") for atoms that contribute basis functions but not nuclear
        charge.

        All coordinate inputs are expected in Angstrom.  All returned matrices
        and tensors are C-contiguous (row-major) NumPy arrays of float64.
    )doc";

    m.def("ao_integrals", &ao_integrals,
        py::arg("atoms"),
        py::arg("basis_set"),
        R"doc(
            Compute S, T, U, ERI for the given geometry and basis set.

            Parameters
            ----------
            atoms : list of (symbol: str, x: float, y: float, z: float, is_ghost: bool)
            basis_set : str

            Returns
            -------
            (S, T, U, ERI) : tuple of np.ndarray
                S   – overlap, shape (nao, nao)
                T   – kinetic energy, shape (nao, nao)
                U   – nuclear attraction, shape (nao, nao)
                ERI – two-electron integrals, shape (nao, nao, nao, nao)
        )doc");

    m.def("overlap_and_kinetic", &overlap_and_kinetic,
        py::arg("atoms"),
        py::arg("basis_set"),
        R"doc(
            Compute S and T only (cheaper when ERI/U are not needed).

            Returns
            -------
            (S, T) : tuple of np.ndarray
        )doc");

    m.def("nuclear_attraction", &nuclear_attraction,
        py::arg("atoms"),
        py::arg("basis_set"),
        R"doc(
            Compute the nuclear attraction matrix U.
            Ghost atoms contribute zero nuclear charge.

            Returns
            -------
            U : np.ndarray, shape (nao, nao)
        )doc");

    m.def("eri4_integrals", &eri4_integrals,
        py::arg("atoms"),
        py::arg("basis_set"),
        R"doc(
            Compute the full ERI4 tensor.

            Returns
            -------
            ERI : np.ndarray, shape (nao, nao, nao, nao)
        )doc");

    m.def("num_aos", &num_aos,
        py::arg("atoms"),
        py::arg("basis_set"),
        R"doc(
            Return the number of spherical AOs for the given geometry and basis.
        )doc");

    m.def("set_basis_path", &set_basis_path,
        py::arg("path"),
        "Override the main basis set directory path at runtime.");

    m.def("set_aux_basis_path", &set_aux_basis_path,
        py::arg("path"),
        "Override the auxiliary basis set directory path at runtime.");
}
