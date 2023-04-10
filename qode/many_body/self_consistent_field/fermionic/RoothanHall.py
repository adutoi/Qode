#    (C) Copyright 2018 Anthony D. Dutoi and Boris D. Gutierrez-Cortes
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
import numpy
from ....util import sort_eigen



# The function compute() below performs the most basic Roothan-Hall logic, delegating the evaluation of the energy and the 
# Fock matrix to a class given in the argument 'energetics'.  All three of these tasks can be performed without knowledge
# of the representation of the integral tensors (orthogonal, nonorthogonal, biorthogonal).  The energetics evaluations do
# need to know the functional and whether the closed or open-shell convention is being used, however (see the choices in
# rSCF_uSCF.py in this directory).  The open/closed convention also affects how the density matrix is constructed from MO
# coefficients, but in a trivial way.  So the argument 'n_MOs' here should be half the number of electrons in the closed-shell
# case (the number of spatial orbitals under consideration) and it should be the number of electrons in the open-shell case
# (the number of spin orbitals under consideration).  Conversely, the algorithms used to obtain a coefficient matrix from a
# Fock matrix and a density matrix from the coefficient matrix depend on the integrals representation (orthogonal,
# nonorthogonal, biorthogonal), but not on the open/closed convention. This dependence is encapsulated in the 'diagonalizer'
# argument given to compute().  The diagonalizer is an instance of one of the classes below that follows compute(), for which
# two variants require initialization using an overlap matrix (to keep these routines internally independent of the open/closed
# convention, they must be set up with an overlap matrix in the same spatial- or spin-orbital basis).  Optionally, a (valid!)
# initial guess 'C' can be given for the MO coefficient matrix (each MO is a column, with appropriate spatial- or spin-orbital
# dimensions, and the first 'n_MOs' columns will be interpreted as occupied -- doubly or singly, depending on the convention).
# If no guess for C is given, the core Hamiltonian will be diagonalized.  The user can also specify the energy convergence
# threshold as 'thresh'.  Note that the density matrix here has maximum eigenvalues of 1 also for the restricted case here,
# as a consequence of insulating most of the code from knowing which convention is used.  Since it is never returned, this is ok.

def compute(n_MOs, diagonalizer, energetics, C, thresh):
    energy, Fock = energetics.energy, energetics.Fock
    iteration = 0
    if C is None:  e,C = diagonalizer(Fock.core())
    D = diagonalizer.rho_from_MOs(C, n_MOs)
    F = Fock(D, iteration)
    E, E_diff = energy(F, D, zeroth_iteration=True)
    while (abs(E_diff) > thresh):       # check also commutator of F and D? (make sure metric independent of rep!)
        iteration += 1
        e, C = diagonalizer(F)
        D = diagonalizer.rho_from_MOs(C, n_MOs)
        F = Fock(D, iteration)
        E, E_diff = energy(F, D)
        #print("Iter=%3d  E = % 16.12f  E_diff = % 8.4e" % (iteration, E, E_diff))
    #print("SCF has finished!\n")
    return E, e, C



# There are three basic representations in which one can imagine performing an SCF, all differing in where and how overlaps
# are handled.
# Orthonormal:     In this case the integrals are given in an orthonormal basis, and the representation of all internal
#                  quantities is in that basis.  There is not need for an overlap matrix.
# Non-orthonormal: In this case the basis need be neither orthogonal nor normalized, but it is assumed to be linearly
#                  indpendent.  All internal quantities are represented as symmetric matrices in this basis, and the
#                  diagonalizer is aware of the overlap matrix.
# Biorthogonal:    In this case the basis need be neither orthogonal nor normalized, and, depending on the sub-case, it
#                  might be linearly dependent (really?).  All internal quantities are represented as non-symmetric matrices
#                  in the biorthogonal bases.  Depending on the subcase, the diagonalizer may or may not be aware of
#                  the overlap matrix (for the right-hand basis).  Here are the subcases:
#                  1) In the more general case, the diagonalizer need not be aware of the overlap matrix, but the price
#                     to be paid for this is that the general diagonalizer used cannot be forced to return eigenvectors
#                     that correspond to an orthonormal set, nor which are even real.  For non-degenerate eigenvalues
#                     the eigenvectors will build orbitals that are orthogonal (but not necessarily normalized).  Eigenvectors
#                     within a degenerate eigenspace, however, may not be orthogonal (but they should be linearly
#                     independent) and they may be fundamentally complex (not just a real vector times a complex scalar),
#                     due to such mixings.  That said, the density matrix will still be real (though not symmetric),
#                     as will any physical quantities derived from it, like the energy, albeit perhaps with an imaginary part
#                     on the order of machine precision.  It will be up to the calling routine to deal with these
#                     difficulties, as any automated solution at this level has the potential to be fragile.
#                     It is this case which can be made to handle linearly-dependent problems, since no overlap matrix
#                     needs to be inverted.  However, in this case the integrals are not well-defined, and all that
#                     complexity is left to the calling routine (which will need to apply a pseudo-inverse to the transforms).
#                  2) If the biorthogonal representation is desired, but the user insists on obtaining a coefficient
#                     matrix that corresponds to real, orthonormal orbitals, the overlap matrix must be supplied.  
# It is worth noting that the second S-dressed biorthogonal option can be thought of as a "wrapper" for the
# non-orthogonal case (in the algebraic sense, not coded as such), and the coefficient matrix will be virtually
# indistinguishable, from what is obtained in that way.  For that matter, the non-orthogonal case is an algebraic wrapper
# over the orthogonal case.  The difference between "algebraic wrapping" and an external code wrapper that a user
# might easily write is that the user's integrals are never transformed with the algebraic wrapping.  The user should
# choose an option based on how the integrals are stored.  In the biorthogonal case, S-dressed is probably the best
# option, if you can supply an overlap matrix.

def inv_rt(M):
    """ takes the -1/2 power of a symmetric matrix M """
    vals, vecs = numpy.linalg.eigh(M)
    vals = numpy.diag([1/numpy.sqrt(v) for v in vals])
    return vecs @ vals @ vecs.T

class orthonormal_diagonalizer:
    def __init__(self): pass
    def __call__(self, F):
        return sort_eigen(numpy.linalg.eigh(F))
    def rho_from_MOs(self, C, n_MOs):
        return C[:,:n_MOs] @ C.T[:n_MOs,:]

class nonorthogonal_diagonalizer:
    def __init__(self, S):
        self.A = inv_rt(S)
    def __call__(self, F):
        A = self.A
        e, U = sort_eigen(numpy.linalg.eigh( A.T @ F @ A ))
        return e, A@U
    def rho_from_MOs(self, C, n_MOs):
        return C[:,:n_MOs] @ C.T[:n_MOs,:]

class biorthogonal_diagonalizer_general:
    def __init__(self):  pass
    def __call__(self, F):
        return sort_eigen(numpy.linalg.eig(F))
    def rho_from_MOs(self, C, n_MOs):
        C_inv = numpy.linalg.inv(C)
        return C[:,:n_MOs] @ C_inv[:n_MOs,:]

class biorthogonal_diagonalizer_Sdressed:
    def __init__(self, S):
        A = inv_rt(S)
        self.args = A, A@S
    def __call__(self, F):
        A, A_inv = self.args
        e, U = sort_eigen(numpy.linalg.eigh( A_inv @ F @ A ))
        return e, A@U
    def rho_from_MOs(self, C, n_MOs):
        C_inv = numpy.linalg.inv(C)
        return C[:,:n_MOs] @ C_inv[:n_MOs,:]
