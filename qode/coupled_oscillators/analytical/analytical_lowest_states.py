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
from qode.coupled_oscillators import pair_dipole_matrix
import copy

class analytical_indices( object ):
    def __init__(self,
                 num_oscillators,
                 states_required,
                 eff_k_list
                 ):
        self.num_oscillators = num_oscillators
        self.states_required = states_required
        self.eff_k_list      = eff_k_list
        self.current_states  = [ 0 for i in range(self.num_oscillators) ]
        self.thresh_config   = [ 0 for i in range(self.num_oscillators) ]
        self.thresh_config[0]= self.states_required - 1
        self.continue_flag   = True
    
    def get_continue_flag( self ):
        return self.continue_flag
    
    def set_continue_flag( self, bool_type ):
        self.continue_flag = bool_type
    
    def round_up_digit( self, digit_i ):
        self.current_states[digit_i]      =  0
        self.current_states[digit_i + 1] +=  1
    
    def round_up( self ):
        n = 0
        check_zero_flag = True
        while check_zero_flag and self.current_states[n] == 0:
            n += 1
            if n > self.num_oscillators - 1:
                check_zero_flag = False
        
        if n >= self.num_oscillators - 1:
            self.set_continue_flag( False )
        else:
            self.round_up_digit(n)
    
    
    def increase_state( self ):
        self.current_states[0] += 1
            
            
    def generate_states( self, energy_func ):
        energy_thresh = energy_func( self.thresh_config, self.eff_k_list )
        print("THRESHOLD =",energy_thresh)
        states        = []
        eigenvalues   = []

        while self.get_continue_flag():
            energy_tmp = energy_func( self.current_states, self.eff_k_list )
            if energy_tmp > energy_thresh:
                self.round_up()
            else:
                states      += [ copy.deepcopy( self.current_states ) ]
                eigenvalues += [ energy_tmp ]
                self.increase_state()
            
        return states, eigenvalues

            
    def __call__(self, energy_func):
        states, eigenvalues = self.generate_states( energy_func )
        return states, eigenvalues



def get_lowest_n_states(
                        eff_k_list,
                        lowest_num_states
                        ):
    # eff_k_list comes from diagonalization of coupling matrix,
    # it's in sorted order but put a guard to ensure this.
    # lowest_num_states is the state number requested for analytical solution
    k_list          = sorted( eff_k_list )
    num_oscillators = len( k_list )
    #
    # Make a function that increases the indices.
    #
    analytical_calc     = analytical_indices( num_oscillators, lowest_num_states, k_list )
    states, eigenvalues = analytical_calc( unitless_HO_system_energy )
    
    # Call the sort_by_energy fucntion in qode.coupled_oscillators.pair_dipole_matrix
    states, eigenvalues = pair_dipole_matrix.sort_by_energy( states, eigenvalues )

    states_required = copy.deepcopy( states[0:lowest_num_states] )
    eigenvalues_required = copy.deepcopy( eigenvalues[0:lowest_num_states] )
    
    return states_required, eigenvalues_required


    



if __name__ == "__main__":
    k = [1.,2.,3.] #[1.1, 1.2, 1.3, 1.4, 1,5, 1.6, 1.7, 1.8, 1.9, 2.0]
    num_states = 10
    states, eigenvalues =  get_lowest_n_states( k, num_states )
    for line in states:
        print(line, unitless_HO_system_energy( line , k ))
    #print(eigenvalues)

