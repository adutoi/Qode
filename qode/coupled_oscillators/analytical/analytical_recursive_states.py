#    (C) Copyright 2018 Anthony D. Dutoi
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
from qode.coupled_oscillators.analytical.analytical_routines import unitless_HO_system_energy
from qode.coupled_oscillators.pair_dipole_matrix import sort_by_energy
# Declaration of sort_by_energy is,  def sort_by_energy( list_of_index, list_of_energy )
import copy




def one_body_low_n_states(
                          n_states,
                          k_value
                          ):
    return [ [i] for i in range(n_states)], [ unitless_HO_system_energy([i],[k_value]) for i in range(n_states) ]




def one_body_below_E_states(
                            energy_limit,
                            k_value
                            ):
    states = []
    energy = []
    n = 0
    tmp_E = unitless_HO_system_energy( [n], [k_value] )
    while tmp_E < energy_limit:
        energy += [ tmp_E ]
        states += [[n]]
        n += 1
        tmp_E = unitless_HO_system_energy( [n], [k_value] )
    
    return states, energy




def nbody_plus_one_get_low_n_states(
                                    n_states,
                                    next_states,
                                    current_states,
                                    next_energy,
                                    current_energy
                                    ):
    new_states = []
    new_energy = []
    for i in range( len(current_energy) ):
        for j in range( len(next_energy) ):
            new_states += [ copy.deepcopy( current_states[i] ) + copy.deepcopy( next_states[j] ) ]
            new_energy += [ current_energy[i]+next_energy[j] ]
    new_states , new_energy = sort_by_energy( new_states, new_energy )
    return copy.deepcopy( new_states[:n_states] ), copy.deepcopy( new_energy[:n_states] )




def get_lowest_n_states(
                        n_states,
                        eff_k_list
                        ):
    num_frag = len( eff_k_list )
    
    # states, energies, energy limit at the very beginning. ( Initialization: ONLY ONE FRAGMENT )
    current_states, current_energy = one_body_low_n_states( n_states, eff_k_list[0] )
    energy_limit = unitless_HO_system_energy( current_states[-1], [eff_k_list[0]] )
    
    for i in range( num_frag - 1 ):
        #print(i,i+1)
        next_states, next_energy = one_body_below_E_states( energy_limit, eff_k_list[i+1] )
        current_states, current_energy = nbody_plus_one_get_low_n_states(
                                                                         n_states,
                                                                         next_states,
                                                                         current_states,
                                                                         next_energy,
                                                                         current_energy
                                                                         )
        
        energy_limit = unitless_HO_system_energy( current_states[-1], eff_k_list[:i+2] )
                                                                         
    
    
    return current_states, current_energy









if __name__ == "__main__":
    k = [1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1]
    num_states = 16
    states , energy = get_lowest_n_states( num_states, k )
    print("FIRST", num_states, "STATES:")
    for i in range( num_states ):
        print( states[i], "=", energy[i] )



