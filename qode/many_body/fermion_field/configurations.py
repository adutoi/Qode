#    (C) Copyright 2023, 2025 Anthony Dutoi
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



def dn_up_elec(n_elec, Sz):
    """ for a given number of electrons and Sz value returns the number that are spin down and spin up """
    n_extra_up = round(2*Sz)
    if abs(n_extra_up - 2*Sz)>1e-12:     # hard-coded threshhold here probably ok in the real world
        raise("Sz value given was not a half-integer")
    n_elec_even = n_elec - n_extra_up    # yes, this all works if Sz<0, too
    if n_elec_even%2!=0:
        raise("Sz value not compatible with number of electrons")
    n_elec_dn =              n_elec_even//2
    n_elec_up = n_extra_up + n_elec_even//2
    return n_elec_dn, n_elec_up

def combine_orb_lists(orbs_A, orbs_B, n_orbs_A):
    """ given lists of indices referencing two different sets of orbitals, return the combined indices referencing a combined set """
    return list(orbs_A) + [n_orbs_A+p for p in orbs_B]



# For understanding the insides of the functions below, it may be helpful to see the notes in __init__.py

def all_configs(num_tot_orb, num_active_elec, frozen_occ_orbs=None, frozen_vrt_orbs=None):
    """ all possible confifigurations of the given number of active electrons in the given orbitals with some orbitals having frozen occupancy """
    def recur_all_configs(active_orbs, num_active_elec, static_config):
        configs = []
        for p in range(num_active_elec-1, len(active_orbs)):
            config = static_config + 2**active_orbs[p]
            if num_active_elec==1:  configs += [config]
            else:                   configs += recur_all_configs(active_orbs[:p], num_active_elec-1, config)
        return configs
    if frozen_occ_orbs is None:  frozen_occ_orbs = []
    if frozen_vrt_orbs is None:  frozen_vrt_orbs = []
    frozen_occ_orbs = set(frozen_occ_orbs)
    frozen_vrt_orbs = set(frozen_vrt_orbs)
    static_config = 0
    active_orbs = []
    for p in range(num_tot_orb):
        if p in frozen_occ_orbs:
            static_config += 2**p
        elif p not in frozen_vrt_orbs:
            active_orbs += [p]
    return recur_all_configs(active_orbs, num_active_elec, static_config)

def Sz_configs(n_spatial, n_elec, Sz, core):
    """ all possible configurations of specified Sz for the given number of electrons in the given number of *spatial* orbitals, also spin down and up parts  """
    # generalize as needed
    n_elec_dn, n_elec_up = dn_up_elec(n_elec, Sz)
    dn_configs = all_configs(n_spatial, n_elec_dn-len(core), frozen_occ_orbs=core)
    up_configs = all_configs(n_spatial, n_elec_up-len(core), frozen_occ_orbs=core)
    configs = tensor_product_configs([dn_configs,up_configs], [n_spatial,n_spatial])
    return configs, (dn_configs, up_configs)



def decompose_configs(configs, orb_counts):
    """ decompose configs into multiply nested lists for subsystems with numbers of orbitals as given, lower-numbered systems (earlier in list) are more deeply nested (ie, fast index)"""
    # has only been tested for two systems
    # keep in mind that the highest-indexed system (Z) is stored in the first bits of a configuration (which is why the higher systems are the "slow" index)
    nested = []
    orb_counts = orb_counts[:-1]
    shift = 2**sum(orb_counts)
    configZ_prev = None
    for config in configs:
        configZ  = config // shift
        configAY = config %  shift
        if configZ!=configZ_prev:
            if configZ_prev is not None:
                if len(orb_counts)>1:  nested += [(decompose_configs(configsAY, orb_counts), configZ_prev)]
                else:                  nested += [(configsAY, configZ_prev)]
            configsAY = []
            configZ_prev = configZ
        configsAY += [configAY]
    if len(orb_counts)>1:  nested += [(decompose_configs(configsAY, orb_counts), configZ_prev)]
    else:                  nested += [(configsAY, configZ_prev)]
    return nested

def config_combination(orb_counts):
    """ returns a function that combines configs from subsystems into a supersystem config, given the numbers of orbitals for each subsystem """
    shifts = [1]
    orb_count_tot = 0
    for orb_count in orb_counts[:-1]:
        orb_count_tot += orb_count
        shifts += [2**orb_count_tot]
    def combine_configs(configsX):
        config = 0
        for configX,shift in zip(configsX,shifts):
            config += configX*shift
        return config
    return combine_configs

def recompose_configs(nested, orb_counts):
    """ given a multiply nested representation of a list of supersystem configurations, given the flattened list of explicit supersystem configurations """
    # has only been tested for two systems
    # keep in mind that the highest-indexed system (Z) is stored in the first bits of a configuration (which is why the higher systems are the "slow" index)
    orb_count  = orb_counts[-1]
    orb_counts = orb_counts[:-1]
    combine_configs = config_combination([sum(orb_counts), orb_count])
    configs = []
    for nestedAY,configZ in nested:
        if len(orb_counts)>1:  configsAY = recompose_configs(nestedAY, orb_counts)
        else:                  configsAY = nestedAY
        for configAY in configsAY:
            configs += [combine_configs([configAY, configZ])]
    return configs

def tensor_product_configs(configsX, orb_counts):
    """ given lists of configurations for each subsystem, given a flattened list of explicit supersystem configurations corresponding to the tensor-product basis """
    def recur_tensor_product_configs(configsX):
        configsZ = configsX[-1]
        configsX = configsX[:-1]
        if len(configsX)==0:
            nested = configsZ
        else:
            nested = []
            for configZ in configsZ:
                nested += [(recur_tensor_product_configs(configsX), configZ)]
        return nested
    return recompose_configs(recur_tensor_product_configs(configsX), orb_counts)



def print_configs(nested, orb_counts, printout=print, _indent=""):
    """ given a list of configurations (maybe in multiply nested representation), print the bit-string representation for each configuration """
    num_orb = orb_counts[-1]
    orb_counts = orb_counts[:-1]
    formatter = "{{:0{}b}}".format(num_orb)
    if len(orb_counts)==0:
        for config in nested:
            printout(_indent + formatter.format(config))
    else:
        for nestedAY,configZ in nested:
            printout(_indent + formatter.format(configZ))    # remember high-index system stored in first (high-order) bits
            print_configs(nestedAY, orb_counts, printout, _indent+" "*num_orb)



if __name__ == "__main__":
    #configs = _all_configs([0,2,4,6,8], 3, 0b0010101010)
    #print_configs(configs, [10])
    #print()
    configs = all_configs(10, 4, frozen_occ_orbs=[0,5], frozen_vrt_orbs=[1,6])
    print_configs(configs, [10])
    print()
    configs = all_configs(6, 3)
    print_configs(configs, [6])
    print()
    nested = decompose_configs(configs, [3,3])
    print_configs(nested, [3,3])
    print()
    configs = recompose_configs(nested, [3,3])
    print_configs(configs, [6])
    print()
    configs = all_configs(4, 2)
    print_configs(configs, [4])
    print()
    #nested = _tensor_pdt_nested([configs, configs])
    #print_configs(nested, [4,4])
    #print()
    configs = tensor_product_configs([configs, configs], [4,4])
    print_configs(configs, [8])
