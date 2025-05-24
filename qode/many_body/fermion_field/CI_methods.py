#    (C) Copyright 2025 Anthony D. Dutoi
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
from ...util import struct, indented
from ...math import linear_inner_product_space
from ...math.lanczos import lowest_eigen
from .CI_space_traits import CI_space_traits
from .field_op_ham import Hamiltonian
from . import configurations



def monomer_configs(frag, Sz=0):
    n_spatial = frag.basis.n_spatial_orb
    n_elec_dn, n_elec_up = configurations.dn_up_elec(frag.n_elec_ref, Sz)
    dn_configs, up_configs = configurations.Sz_configs(n_spatial, n_elec_dn, n_elec_up, frag.basis.core)
    return configurations.tensor_product_configs([dn_configs,up_configs], [n_spatial,n_spatial])

def dimer_configs(frag0, frag1, Sz=0):
    n_spatial_0 = frag0.basis.n_spatial_orb
    n_spatial_1 = frag1.basis.n_spatial_orb
    n_spatial   = n_spatial_0 + n_spatial_1
    core        = configurations.combine_orb_lists(frag0.basis.core, frag1.basis.core, n_spatial_0)

    n_elec_dn,  n_elec_up  = configurations.dn_up_elec(frag0.n_elec_ref + frag1.n_elec_ref, Sz)
    dn_configs, up_configs = configurations.Sz_configs(n_spatial, n_elec_dn, n_elec_up, core)
    dn_configs_decomp = configurations.decompose_configs(dn_configs, [n_spatial_0, n_spatial_1])
    up_configs_decomp = configurations.decompose_configs(up_configs, [n_spatial_0, n_spatial_1])

    nested = []
    combine_configs = configurations.config_combination([n_spatial_1, n_spatial_1])
    for dn_config_1,dn_configs_0 in dn_configs_decomp:
        for up_config_1,up_configs_0 in up_configs_decomp:
            config_1  = combine_configs([up_config_1, dn_config_1])
            configs_0 = configurations.tensor_product_configs([up_configs_0, dn_configs_0], [n_spatial_0, n_spatial_0])
            nested += [(config_1, configs_0)]
    configs = configurations.recompose_configs(nested, [2*n_spatial_0, 2*n_spatial_1])

    return configs, nested



def lanczos_ground(integrals, configs, occupied, thresh=None, printout=print, n_threads=1):
    options = struct(printout=indented(printout))
    if thresh is not None:         # if not defined/passed forward ...
        options.thresh = thresh    # ... default from lanczos takes over

    N, h, V = integrals
    CI_space = linear_inner_product_space(CI_space_traits(configs))

    H     = CI_space.lin_op(Hamiltonian(h,V, n_threads=n_threads))
    guess = CI_space.member(CI_space.aux.basis_vec(occupied))

    printout("Energy of guess =", (guess|H|guess) + N)
    (Eval,Evec), = lowest_eigen(H, [guess], **options)
    printout("Ground-state energy = ", Eval+N)

    return Eval, Evec
