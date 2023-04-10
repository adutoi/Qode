#    (C) Copyright 2019 Anthony D. Dutoi and Simeng Zhang
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
from .operator import operator
from .amp_eqs import amp_eqs1, amp_eqs2, energy



class BakerCampbellHausdorff(object):
    """ This class is a wrapper for the operator-specific dirty-work.  It has an interface needed by the coupled_cluster module, but hides the internal nature of the operator.  """
    def __init__(self, n_occ, n_vrt, resources=None):	# temporary  hack to get things working ... ignore resources
        self.dims = n_occ, n_vrt
        self.resources = resources
    def computeOmega(self, H, T, textlog):  # call must have this signature!   Computes excitation part (0th, 1st and 2nd order) of HBar
        n_occ, n_vrt = self.dims
        F = H.scalars[1]
        Foo = F[:n_occ,:n_occ]
        Fov = F[:n_occ,n_occ:]
        Fvo = F[n_occ:,:n_occ]
        Fvv = F[n_occ:,n_occ:]
        V = H.scalars[2]
        Voooo = V[:n_occ,:n_occ,:n_occ,:n_occ]
        Vooov = V[:n_occ,:n_occ,:n_occ,n_occ:]
        Voovo = V[:n_occ,:n_occ,n_occ:,:n_occ]
        Voovv = V[:n_occ,:n_occ,n_occ:,n_occ:]
        Vovoo = V[:n_occ,n_occ:,:n_occ,:n_occ]
        Vovov = V[:n_occ,n_occ:,:n_occ,n_occ:]
        Vovvo = V[:n_occ,n_occ:,n_occ:,:n_occ]
        Vovvv = V[:n_occ,n_occ:,n_occ:,n_occ:]
        Vvooo = V[n_occ:,:n_occ,:n_occ,:n_occ]
        Vvoov = V[n_occ:,:n_occ,:n_occ,n_occ:]
        Vvovo = V[n_occ:,:n_occ,n_occ:,:n_occ]
        Vvovv = V[n_occ:,:n_occ,n_occ:,n_occ:]
        Vvvoo = V[n_occ:,n_occ:,:n_occ,:n_occ]
        Vvvov = V[n_occ:,n_occ:,:n_occ,n_occ:]
        Vvvvo = V[n_occ:,n_occ:,n_occ:,:n_occ]
        Vvvvv = V[n_occ:,n_occ:,n_occ:,n_occ:]
        Omega = operator(T.order, n_vrt, n_occ, empty=True)
        Omega.scalars[0] = energy(Fov, Voovv, T.scalars[1], T.scalars[2])
        Omega.scalars[1] = amp_eqs1(Foo, Fov, Fvo, Fvv,
             Voooo,Vooov,Voovo,Vovoo,Vvooo,Voovv,Vovov,Vvoov,Vovvo,Vvovo,Vvvoo,Vovvv,Vvovv,Vvvov,Vvvvo, Vvvvv,
             T.scalars[1], T.scalars[2])
        Omega.scalars[2] = amp_eqs2(Foo, Fov, Fvo, Fvv,
             Voooo,Vooov,Voovo,Vovoo,Vvooo,Voovv,Vovov,Vvoov,Vovvo,Vvovo,Vvvoo,Vovvv,Vvovv,Vvvov,Vvvvo, Vvvvv,
             T.scalars[1], T.scalars[2])
        return Omega
    @staticmethod
    def Energy(Omega, textlog):  # call must have this signature!   Retrieves energy from excitation part (0th order component) of HBar
        return Omega.scalars[0]
