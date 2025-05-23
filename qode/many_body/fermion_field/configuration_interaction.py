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
from ..util import struct, indented
from ..math import linear_inner_product_space
from ..math.lanczos import lowest_eigen
from ..fermion_field import CI_space_traits, field_op_ham

def lanczos_lowest(integrals, configs, occupied, thresh=None, printout=print, n_threads=1):
    options = struct(printout=indented(printout))
    if thresh is not None:         # if not defined/passed forward ...
        options.thresh = thresh    # ... default from lanczos takes over

    N, h, V = integrals
    CI_space = linear_inner_product_space(CI_space_traits(configs))

    H     = CI_space.lin_op(field_op_ham.Hamiltonian(h,V, n_threads=n_threads))
    guess = CI_space.member(CI_space.aux.basis_vec(occupied))

    printout("Energy of guess =", (guess|H|guess) + N)
    (Eval,Evec), = lowest_eigen(H, [guess], **options)
    printout("Ground-state energy = ", Eval+N)

    return Eval, Evec




#    (C) Copyright 2023, 2024, 2025 Anthony D. Dutoi
# 
#    This file is part of QodeApplications.
# 
#    QodeApplications is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    QodeApplications is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with QodeApplications.  If not, see <http://www.gnu.org/licenses/>.
#
from qode.fermion_field import configurations



def _dn_up_elec(n_elec, Sz):
    n_extra_up = round(2*Sz)
    if abs(n_extra_up - 2*Sz)>1e-12:     # hard-coded threshhold here probably ok
        raise("Sz value given was not a half-integer")
    n_elec_even = n_elec - n_extra_up    # yes, this all works if Sz<0, too
    if n_elec_even%2!=0:
        raise("Sz value not compatible with number of electrons")
    n_elec_dn =              n_elec_even//2
    n_elec_up = n_extra_up + n_elec_even//2
    return n_elec_dn, n_elec_up

def _Sz_configs(n_spatial, n_elec_dn, n_elec_up, core):
    dn_configs = configurations.all_configs(n_spatial, n_elec_dn-len(core), frozen_occ_orbs=core)
    up_configs = configurations.all_configs(n_spatial, n_elec_up-len(core), frozen_occ_orbs=core)
    return dn_configs, up_configs



def combine_orb_lists(orbs_A, orbs_B, n_orbs_A):
    return list(orbs_A) + [n_orbs_A+p for p in orbs_B]

def monomer_configs(frag, Sz=0):
    n_spatial = frag.basis.n_spatial_orb
    n_elec_dn, n_elec_up = dn_up_elec(frag.n_elec_ref, Sz)
    dn_configs, up_configs = _Sz_configs(n_spatial, n_elec_dn, n_elec_up, frag.basis.core)
    return configurations.tensor_product_configs([dn_configs,up_configs], [n_spatial,n_spatial])

def dimer_configs(frag0, frag1, Sz=0):
    n_spatial_0 = frag0.basis.n_spatial_orb
    n_spatial_1 = frag1.basis.n_spatial_orb
    n_spatial   = n_spatial_0 + n_spatial_1
    core        = combine_orb_lists(frag0.basis.core, frag1.basis.core, n_spatial_0)

    n_elec_dn,  n_elec_up  = dn_up_elec(frag0.n_elec_ref + frag1.n_elec_ref, Sz)
    dn_configs, up_configs = _Sz_configs(n_spatial, n_elec_dn, n_elec_up, core)
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
