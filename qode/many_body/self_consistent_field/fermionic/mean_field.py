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



# the stuff in here defined what is meant by the "mean field" (restricted or not, DFT/HF flavor)

# These are the energy evaluation and Fock build classes that correspond to restricted Hartree--Fock (RHF).  Each is initialized
# with the molecular integrals in the spatial-orbital basis, but in any representation (orthogonal, nonorthogonal,
# biorthogonal), and each takes the density matrix (and Fock matrix for the energy) in that same representation as
# a call argument.  The 2e- integrals are supposed to be raw (not antisymmetrized) with the first two indices being
# those for ket and the latter two indices being those for the bra (with the electrons in the same order).  The Fock
# builder may also be initialized with a pair of parameters in a tuple, where the first is the iteration number at
# which damping starts (cannot be <1) and the second is the fraction of the most recent prior Fock build to add to
# each iteration.  The HartreeFock class just wraps these two into one object for convenience.
#
# I envision a number of other of these for DFT flavors, where the functional is a constructor argument (on par with the integrals)

class HF_Energy(object):
    def __init__(self, h, E_factor):
        self.args = h, E_factor
    def __call__(self, F, D, zeroth_iteration=False):
        h, E_factor = self.args
        E = E_factor * numpy.tensordot( h+F, D, axes=(([1,0],[0,1])) ) 	# A little convoluted see Szabo and Ostlund Eq. (3.184)
        if zeroth_iteration:  E_diff = float('inf')
        else:                 E_diff = E - self.E_old
        self.E_old = E
        return E, E_diff

class HF_FockMatrix(object):
    def __init__(self, h, v, damp, J_factor):
        self.args = h, v, damp, J_factor
    def __call__(self, D, iteration):
        h, v, damp, J_factor = self.args
        J = numpy.tensordot(v, D, axes=([3,1],[0,1]))
        K = numpy.tensordot(v, D, axes=([2,1],[0,1]))
        F_new = h + J_factor*J - K
        F = F_new
        if damp is not None:
            damp_start, damp_value = damp
            if (iteration>0) and (iteration>=damp_start):	# even if the user says to start damping at 0, there is nothing there (or it is garbage from a previous usage)
                F *= 1-damp_value
                F += damp_value*self.F_old
        self.F_old = F_new
        return F
    def core(self):
        h, v, damp, J_factor = self.args
        return h

class HartreeFock(object):
    def __init__(self, h, v, damp, E_factor, J_factor):
        self.energy = HF_Energy(h, E_factor)
        self.Fock =   HF_FockMatrix(h, v, damp, J_factor)
    @staticmethod
    def restricted(  h, v, damp=None):  return HartreeFock(h, v, damp, E_factor=1,    J_factor=2)
    @staticmethod
    def unrestricted(h, v, damp=None):  return HartreeFock(h, v, damp, E_factor=1/2., J_factor=1)


