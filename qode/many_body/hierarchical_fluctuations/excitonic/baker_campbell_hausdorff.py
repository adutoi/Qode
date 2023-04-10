#    (C) Copyright 2018, 2019 Yuhong Liu and Anthony D. Dutoi
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
from ....util.PyC import import_C, BigInt
from .operator    import RSD_operator
BCH = import_C("bch", include=["partition_storage.h"], flags="-O3")



class BakerCampbellHausdorff(object):
    """ This class is a wrapper for the operator-specific dirty-work.  It has an interface needed by the coupled_cluster module, but hides the internal nature of the operator.  """
    def __init__(self, states_per_frag, resources):
        self.states_per_frag = states_per_frag
        self.resources = resources
    def computeOmega(self, H, T, textlog):  # call must have this signature!   Computes excitation part (0th, 1st and 2nd order) of HBar
        states_per_frag = self.states_per_frag
        N_frag = len(states_per_frag)
        Nvrt = numpy.array([n_st-1 for n_st in states_per_frag], dtype=BigInt.numpy)
        Omega = RSD_operator(states_per_frag)
        Hb = H.blocks
        Tb = T.blocks
        Ob = Omega.blocks
        BCH.execute(self.resources.n_cores, N_frag, Nvrt, \
            Hb["K"], Hb["Ex"], Hb["Fo"], Hb["Fv"], Hb["Dx"], Hb["ExEx"], Hb["ExFo"], Hb["ExFv"], Hb["ExDx"], Hb["FoFo"], Hb["FoFv"], Hb["FoDx"], Hb["FvFv"], Hb["FvDx"], Hb["DxDx"], \
            Tb["Ex"], Tb["ExEx"], \
            Ob["K"], Ob["Ex"], Ob["ExEx"])
        return Omega
    @staticmethod
    def Energy(Omega, textlog):  # call must have this signature!   Retrieves energy from excitation part (0th order component) of HBar
        return Omega.blocks["K"][0]
