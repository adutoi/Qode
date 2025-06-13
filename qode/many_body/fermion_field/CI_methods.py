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
from ...math.lanczos import lowest_eigen, lowest_eigen_one_by_one
from .CI_space_traits import CI_space_traits
from .field_op_ham import Hamiltonian
from . import configurations



def monomer_configs(frag, Sz=0):
    configs, _ = configurations.Sz_configs(frag.basis.n_spatial_orb, frag.n_elec_ref, Sz, frag.basis.core)
    return configs    # just need complete spin-orbital configs

def dimer_configs(frag0, frag1, Sz=0):
    n_spatial_0 = frag0.basis.n_spatial_orb
    n_spatial_1 = frag1.basis.n_spatial_orb
    n_spatial   = n_spatial_0 + n_spatial_1
    core        = configurations.combine_orb_lists(frag0.basis.core, frag1.basis.core, n_spatial_0)    # core in dimer spatial-orb representation

    _, (dn_configs, up_configs) = configurations.Sz_configs(n_spatial, frag0.n_elec_ref+frag1.n_elec_ref, Sz, core)    # get individual spin-dn/up configs
    dn_configs_decomp = configurations.decompose_configs(dn_configs, [n_spatial_0, n_spatial_1])    # break spin-dn/up configurations into fragments ...
    up_configs_decomp = configurations.decompose_configs(up_configs, [n_spatial_0, n_spatial_1])    # ... Gives nested representation, with frag1 as slow index

    # To understand the order/structure of this looping, it may be helpful to see the notes in __init__.py
    nested = []
    combine_configs = configurations.config_combination([n_spatial_1, n_spatial_1])
    for up_configs_0,up_config_1 in dn_configs_decomp:        # loop over all (frag-decomposed) combinations of spin-up ...
        for dn_configs_0,dn_config_1 in up_configs_decomp:    # ... and spin-down configurations (spin-up first because these are high-order bits)
            configs_0 = configurations.tensor_product_configs([dn_configs_0, up_configs_0], [n_spatial_0, n_spatial_0])    # a list of configs
            config_1  = combine_configs([dn_config_1, up_config_1])                                                        # a single config
            nested += [(configs_0, config_1)]
    configs = configurations.recompose_configs(nested, [2*n_spatial_0, 2*n_spatial_1])

    return configs, nested



def _act_Sops(Sops, state_in):
    state = 1 * state_in
    for Sop in Sops:
        Sop_state = Sop|state
        state += Sop_state
    return state

def lanczos_ground(integrals, configs, occupied, n_states=1, thresh=None, printout=print, n_threads=1, full_ints=None):
    options = struct(printout=indented(printout))
    if thresh is not None:         # if not defined/passed forward ...
        options.thresh = thresh    # ... default from lanczos takes over

    N, h, V  = integrals("N h V")
    CI_space = linear_inner_product_space(CI_space_traits(configs))
    H        = CI_space.lin_op(Hamiltonian(h,V, n_threads=n_threads))
    guess    = CI_space.member(CI_space.aux.basis_vec(occupied))

    energy = (guess|H|guess)
    printout(f"Projective energy of guess = {energy} -> {energy+N}")
    if full_ints is not None:
        Hfull = CI_space.lin_op(Hamiltonian(full_ints.h, full_ints.V, n_threads=n_threads))
        Sops = [CI_space.lin_op(Hamiltonian(sig1, sig2, n_threads=n_threads)) for sig1,sig2 in zip(full_ints.sigma_1, full_ints.sigma_2)]
        Sguess  = _act_Sops(Sops, guess)
        SHguess = _act_Sops(Sops, Hfull|guess)
        norm   = (guess|Sguess)
        energy = (guess|SHguess) / norm
        printout(f"Proper norm of guess = {norm}")
        printout(f"Proper expectation energy of guess = {energy} -> {energy+N}")

    if n_states==1:
        results = lowest_eigen(H, [guess], **options)    # in theory, can just go to "else" directly, but out of caution, leave old functionality alone for now
    else:
        results = lowest_eigen_one_by_one(H, [guess]*n_states, **options)    # a little weird to demand multiple states from one guess, but can refine later

    for i,(Eval,Evec) in enumerate(results):
        printout(f"Projective Energy eigenvalue for state {i} = {Eval} -> {Eval+N}")
        if full_ints is not None:
            SEvec  = _act_Sops(Sops, Evec)
            SHEvec = _act_Sops(Sops, Hfull|Evec)
            norm   = (Evec|SEvec)
            energy = (Evec|SHEvec) / norm
            printout(f"Proper norm of state {i} = {norm}")
            printout(f"Proper expectation energy of state {i} = {energy} -> {energy+N}")

    printout("Raw (projective) overlap matrix:")
    for _1,bra in results:
        for _2,ket in results:
            printout(f"{(bra|ket): .10e}", end="  ")
        printout()

    return results
