import numpy as np
import veloxchem as vlx
import timeit

def brabraketket(braketbraket, printout=print):
    start_time = timeit.default_timer()
    n_orb = braketbraket.shape[0]
    V = np.zeros((n_orb, n_orb, n_orb, n_orb))
    for p in range(n_orb):
        for q in range(n_orb):
            for r in range(n_orb):
                for s in range(n_orb):
                    V[p,q,r,s] = braketbraket[p,r,q,s]
    elapsed = timeit.default_timer() - start_time
    printout('elapsed time (brabraketket function) =', elapsed)
    return V

def get_rev_map_of_bfs(molecule, basis):  # adapted from veloxchem get_basis_function_indices_of_atoms
    natoms = molecule.number_of_atoms()
    aoinds_atoms = [[] for atomidx in range(natoms)]

    max_angl = basis.max_angular_momentum()

    aoidx = 0

    for angl in range(max_angl + 1):
        indices = [[] for atomidx in range(natoms)]

        for s in range(-angl, angl + 1):

            for atomidx in range(natoms):
                indices[atomidx].append([])
                nao = basis.number_of_basis_functions([atomidx], angl)

                for i in range(nao):
                    indices[atomidx][-1].append(aoidx)
                    aoidx += 1

        for atomidx in range(natoms):

            reordered_indices = []
            for s in range(len(indices[atomidx])):
                reordered_indices.append(indices[atomidx][s])

            aoinds_atoms[atomidx] += list(
                np.array(reordered_indices).reshape(-1))

    flat_inds = []
    for atomidx in range(natoms):
        flat_inds += aoinds_atoms[atomidx]

    return np.array(flat_inds)

class supplementary(object):
    def __init__(self, mol, bas, ints):
        self.mol, self.bas, self.ints = mol, bas, ints

n_calls = 0

def AO_ints(geometry, basis_str, max_mem=1e9, NucPotentialOnly=False, printout=print):
    print(geometry)
    start_time = timeit.default_timer()
    global n_calls
    n_calls += 1
    printout("Integrals call #:", n_calls, flush=True)
    if "@" in geometry:
        # veloxchem is a bit different with ghost atoms, since one creates the basis and molecule classes separately anyway,
        # and it provides the option for a placeholder atom (BQ) in the molecule object, which, combined with
        # the basis object for the "full" molecule, yields the desired ghost atom.
        mol_with_ghost = vlx.Molecule.read_str(geometry.replace("@", ""))
        mol_with_placeholder = []
        for at in geometry.split("\n"):
            if "@" in at:
                tmp_list = []
                for elem in at.split(" "):
                    if "@" in elem:
                        tmp_list.append("BQ")
                    else:
                        tmp_list.append(elem)
                mol_with_placeholder.append(" ".join(tmp_list))
            else:
                mol_with_placeholder.append(at)
        mol = vlx.Molecule.read_str("\n".join(mol_with_placeholder))
        basis = vlx.MolecularBasis.read(mol_with_ghost, basis_str)
    else:
        mol = vlx.Molecule.read_str(geometry)
        basis = vlx.MolecularBasis.read(mol, basis_str)
    naos = basis.get_dimensions_of_basis()
    mem = naos**4 * 8	# assuming 8 bytes per integral
    if mem>max_mem:  raise ValueError("Raise max_mem to above {} (bytes) to store this many integrals.".format(mem))
    U = vlx.compute_nuclear_potential_integrals(mol, basis)
    if NucPotentialOnly:
        S = None
        T = None
        V = None
        # since this is the only time, where multiple fragments are involved, the following
        # reordering only needs to be applied here.
        # veloxchem orders its AOs in "ascending ang mom order"
        # type S S S S S S S S ... P-1 P-1 ... P-1 P-1 P_0 P_0 ... P_0 P_0 P+1 P+1 ... P+1 P+1 D-2 D-2 ...
        # atom A A A A B B B B ... A   A   ... B   B   A   A   ... B   B   A   A   ... B   B   A   A   ...
        # but we need atom A A A A ... A A B B ... B B B B
    else:
        S = vlx.compute_overlap_integrals(mol, basis)
        T = vlx.compute_kinetic_energy_integrals(mol, basis)
        ERIs = vlx.FockDriver().compute_eri(mol, basis)
        V = brabraketket(ERIs, printout)
    if "@" in geometry:
        rev_map = get_rev_map_of_bfs(mol, basis)
        if NucPotentialOnly:
            mats = [U]
            U = U[:, rev_map]
            U = U[rev_map]
        else:
            mats = [S,T,U]
            for i in range(len(mats)):
                mats[i] = mats[i][:, rev_map]
                mats[i] = mats[i][rev_map]
            S,T,U = mats
            V = V[:, :, :, rev_map]
            V = V[:, :, rev_map]
            V = V[:, rev_map]
            V = V[rev_map]
    elapsed = timeit.default_timer() - start_time
    printout('elapsed time (call to integrals {}) ='.format(n_calls), elapsed)
    return S,T,U,V,supplementary(mol,basis,None)
