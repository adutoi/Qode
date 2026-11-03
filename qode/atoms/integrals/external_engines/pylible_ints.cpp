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
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <algorithm>
#include <memory>

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
// RI (density fitting) in batches of auxiliary shells
// ---------------------------------------------------------------------------

/// Structure with the main and the auxiliary basis set (ghost atoms as in buildStructure: basis functions of the
/// real element, no nuclear charge; the auxiliary functions are placed on ghost atoms as well).
static lints::Structure buildStructureRI(
    const std::vector<std::tuple<std::string, double, double, double, bool>> &atoms,
    const std::string &basis_set, const std::string &aux_basis_set)
{
    std::vector<int>                       atomic_nrs;
    std::vector<std::array<double, 3>>     coords_ang;
    lints::basis_atoms_t                   basis_atoms, aux_atoms;

    for (const auto &[symbol, x, y, z, is_ghost] : atoms)
    {
        int real_Z = resolveAtomicNr(symbol);
        bool ghost = is_ghost || (!symbol.empty() && symbol[0] == '@');
        int structure_Z = ghost ? 0 : real_Z;
        atomic_nrs.push_back(structure_Z);
        coords_ang.push_back({x, y, z});

        lints::BasisAtom ba = lints::basisForAtom(real_Z, basis_set);
        ba.atomic_nr_ = structure_Z;
        basis_atoms.push_back(ba);
        lints::BasisAtom aux = lints::basisForAtomAux(real_Z, aux_basis_set);
        aux.atomic_nr_ = structure_Z;
        aux_atoms.push_back(aux);
    }
    return lints::Structure(basis_atoms, aux_atoms, atomic_nrs, coords_ang);
}

/**
 * Three-index integrals (mu nu|P) in batches of auxiliary shells, so that the AO three-index tensor is never held
 * as a whole: the caller transforms every batch to its orbitals right away (e.g. (pq|P) = C^T (..|P) C) and fits it
 * with the Coulomb metric (P|Q) of eri2().  The auxiliary shells are numbered 0 ... n_aux_shells()-1 over all
 * angular momenta; aux_shell_functions(i) gives the global auxiliary-function indices of shell i.
 */
class RIEngine
{
public:
    RIEngine(const std::vector<std::tuple<std::string, double, double, double, bool>> &atoms,
             const std::string &basis_set, const std::string &aux_basis_set)
        : structure_(buildStructureRI(atoms, basis_set, aux_basis_set))
    {
        sh_datas_ = lints::shellDataAux(structure_);
        sp_datas_ = lints::shellPairData(true, structure_);
        for (size_t g = 0; g < sh_datas_.size(); g++)
            for (size_t s = 0; s < sh_datas_[g].n_shells_; s++)
                shells_.push_back({g, s});
        kernels_.resize(sp_datas_.size() * sh_datas_.size());
    }

    size_t numAOs() const { return structure_.getDimAO(); }
    size_t numAux() const { return structure_.getDimAOAux(); }
    size_t numAuxShells() const { return shells_.size(); }

    std::vector<size_t> auxShellFunctions(size_t i) const
    {
        const auto &[g, s] = shells_.at(i);
        const lints::ShellData &sh = sh_datas_[g];
        size_t ofs = sh.offsets_sph_[s];
        std::vector<size_t> out;
        for (int c = 0; c < lints::numSphericals(sh.l_); c++)
            out.push_back(ofs + c);
        return out;
    }

    py::array_t<double> eri2() const
    {
        lible::vec2d m;
        {
            py::gil_scoped_release release;
            m = lints::eri2(structure_);
        }
        return vec2d_to_numpy(m);
    }

    /// (mu nu|P) for the auxiliary shells first ... last-1, as an array (n_functions, nao, nao) whose first index
    /// runs over the functions of the shells in order (see auxShellFunctions).
    py::array_t<double> eri3Batch(size_t first, size_t last)
    {
        if (first > last || last > shells_.size())
            throw std::runtime_error("RIEngine.eri3_batch: shell range out of bounds");
        size_t nao = numAOs();
        std::vector<size_t> local_ofs;    // first local function index of every shell of the batch
        size_t n_func = 0;
        for (size_t i = first; i < last; i++)
        {
            local_ofs.push_back(n_func);
            n_func += lints::numSphericals(sh_datas_[shells_[i].first].l_);
        }
        py::array_t<double> arr({n_func, nao, nao});
        double *buf = arr.mutable_data();
        {
            py::gil_scoped_release release;
            std::fill(buf, buf + n_func * nao * nao, 0.0);
            for (size_t isp = 0; isp < sp_datas_.size(); isp++)
            {
                const lints::ShellPairData &sp = sp_datas_[isp];
                for (size_t i = first; i < last; i++)
                {
                    const size_t g = shells_[i].first, s = shells_[i].second;
                    const lints::ShellData &sh = sh_datas_[g];
                    const lints::ERI3Kernel &kernel = getKernel(isp, g);
                    const size_t lofs = local_ofs[i - first];
#pragma omp parallel for schedule(dynamic)
                    for (size_t ipair = 0; ipair < sp.n_pairs_; ipair++)
                    {
                        lible::vec3d batch = kernel(ipair, s, sp, sh);
                        size_t ofs_a = sp.offsets_sph_[2 * ipair];
                        size_t ofs_b = sp.offsets_sph_[2 * ipair + 1];
                        for (size_t ia = 0; ia < batch.dim<0>(); ia++)
                            for (size_t ib = 0; ib < batch.dim<1>(); ib++)
                                for (size_t ic = 0; ic < batch.dim<2>(); ic++)
                                {
                                    double v = batch(ia, ib, ic);
                                    size_t P = lofs + ic, mu = ofs_a + ia, nu = ofs_b + ib;
                                    buf[(P * nao + mu) * nao + nu] = v;
                                    buf[(P * nao + nu) * nao + mu] = v;
                                }
                    }
                }
            }
        }
        return arr;
    }

private:
    const lints::ERI3Kernel &getKernel(size_t isp, size_t g)
    {
        auto &k = kernels_[isp * sh_datas_.size() + g];
        if (!k)
            k = std::make_unique<lints::ERI3Kernel>(sp_datas_[isp], sh_datas_[g]);
        return *k;
    }

    lints::Structure structure_;
    std::vector<lints::ShellData> sh_datas_;
    std::vector<lints::ShellPairData> sp_datas_;
    std::vector<std::pair<size_t, size_t>> shells_;    // (shell data group, shell in group)
    std::vector<std::unique_ptr<lints::ERI3Kernel>> kernels_;
};

static std::vector<std::string> available_basis_sets()
{
    auto s = lints::availableBasisSets();
    return std::vector<std::string>(s.begin(), s.end());
}

static std::vector<std::string> available_aux_basis_sets()
{
    auto s = lints::availableBasisSetsAux();
    return std::vector<std::string>(s.begin(), s.end());
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

    py::class_<RIEngine>(m, "RIEngine",
        "Three-index integrals (mu nu|P) in batches of auxiliary shells and the metric (P|Q) (RI / density fitting)")
        .def(py::init<const std::vector<std::tuple<std::string, double, double, double, bool>> &,
                      const std::string &, const std::string &>(),
             py::arg("atoms"), py::arg("basis_set"), py::arg("aux_basis_set"))
        .def("num_aos", &RIEngine::numAOs)
        .def("num_aux", &RIEngine::numAux)
        .def("num_aux_shells", &RIEngine::numAuxShells)
        .def("aux_shell_functions", &RIEngine::auxShellFunctions, py::arg("shell"),
             "global auxiliary-function indices of an auxiliary shell")
        .def("eri2", &RIEngine::eri2, "the Coulomb metric (P|Q), shape (naux, naux)")
        .def("eri3_batch", &RIEngine::eri3Batch, py::arg("first"), py::arg("last"),
             "(mu nu|P) of the auxiliary shells first ... last-1, shape (n_functions, nao, nao)");

    m.def("available_basis_sets", &available_basis_sets, "names of the main basis sets known to lible");
    m.def("available_aux_basis_sets", &available_aux_basis_sets, "names of the auxiliary basis sets known to lible");

    m.def("set_aux_basis_path", &set_aux_basis_path,
        py::arg("path"),
        "Override the auxiliary basis set directory path at runtime.");
}
