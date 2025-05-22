import numpy
import psi4
import timeit

def brabraketket(braketbraket, printout=print):
    start_time = timeit.default_timer()
    n_orb = braketbraket.shape[0]
    V = numpy.zeros((n_orb, n_orb, n_orb, n_orb))
    for p in range(n_orb):
        for q in range(n_orb):
            for r in range(n_orb):
                for s in range(n_orb):
                    V[p,q,r,s] = braketbraket[p,r,q,s]
    elapsed = timeit.default_timer() - start_time
    printout('elapsed time (brabraketket function) =', elapsed)
    return V

class supplementary(object):
    def __init__(self, mol, bas, ints):
        self.mol, self.bas, self.ints = mol, bas, ints

n_calls = 0

def AO_ints(geometry, basis, max_mem=1e9, NucPotentialOnly=False, printout=print):
    start_time = timeit.default_timer()
    global n_calls
    n_calls += 1
    printout("Integrals call #:", n_calls, flush=True)
    mol  = psi4.geometry(geometry)
    mol.update_geometry()
    if basis.startswith("CUSTOM\n"):
        basis = basis[7:]
        psi4.basis_helper("assign custombasis\n[custombasis]\n"+basis)	# same syntax as input file:  http://www.psicode.org/psi4manual/1.2/basissets.html#user-defined-basis-sets
        bas = psi4.core.BasisSet.build(mol, quiet=True)
    else:
        psi4.basis_helper("assign "+basis+"\n")
        bas = psi4.core.BasisSet.build(mol, quiet=True)
        #bas = psi4.core.BasisSet.build(mol, target=basis)
    ints = psi4.core.MintsHelper(bas)
    mem = ints.nbf()**4 * 8	# assuming 8 bytes per integral
    if mem>max_mem:  raise ValueError("Raise max_mem to above {} (bytes) to store this many integrals.".format(mem))
    if NucPotentialOnly:
        S = None
        U = numpy.array(ints.ao_potential())
        T = None
        V = None
    else:
        S = numpy.array(ints.ao_overlap())
        U = numpy.array(ints.ao_potential())
        T = numpy.array(ints.ao_kinetic())
        V = brabraketket(numpy.array(ints.ao_eri()), printout)
    elapsed = timeit.default_timer() - start_time
    printout('elapsed time (call to integrals {}) ='.format(n_calls), elapsed)
    return S,T,U,V,supplementary(mol,bas,ints)
