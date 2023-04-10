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
import sys
import numpy as np
import copy
from math import sqrt
from qode.coupled_oscillators.analytical import analytical_solution
from qode.coupled_oscillators import pair_dipole_matrix

class hamiltonian(object):
    "Base Class for Any Hamiltonian"
    pass


#
#
#  top_level_rec_instruct: top level recouple calculation instructions
#
#     This is mandatory for each recouple calculation. User has to supply
#     a 1D list of either int or float type positive numbers.
#       1) int   : indicates   NUMBER OF STATES    for this fragment
#       2) float : indicates   ENERGY UPPER BOUND  for this fragment.



class RecSys_HO_hamiltonian( hamiltonian ):
    def __init__( self, 
                  total_coupling_mat     ,
                  list_of_groups         ,
                  top_level_coupling     ,
                  top_level_rec_instruct ,   
                  unit_length=1.0        ,
                  unit_mass  =1.0        ,
                  unit_hbar  =1.0
               ):
        # Set up parameters: 
        self.list_of_groups     = list_of_groups
        self.total_coupling_mat = total_coupling_mat
        self.top_level_coupling = top_level_coupling
        self.num_fragment       = len( self.top_level_coupling )
        self.all_grp_masses     = [ grp.list_of_masses() for grp in self.list_of_groups ]
        self.rec_instruction    = top_level_rec_instruct
        self.unit_length        = unit_length
        self.unit_mass          = unit_mass
        self.unit_hbar          = unit_hbar
        # Initialize variables to hold analytical results
        self.all_grp_effect_k      = []
        self.all_grp_states_num    = []
        self.all_grp_eigen_energy  = []
        self.all_grp_eigen_vector  = []
        # Loop Over Groups ( Analytical Solutions )
        loop_index = 0        
        for grp in self.list_of_groups:
            if type( self.rec_instruction[ loop_index ] ) == int:
                #print("\nFragment",loop_index,"requested for", self.rec_instruction[ loop_index ],"States")
                # Number of States is Required in rec_instruction
                # Make a guess for the energy_upper_bound
                # Do analytical calculation and check if it got enough states
                # If not then slightly increase the energy_upper_bound.
                rec_states = self.rec_instruction[ loop_index ]
                eff_k_mat = analytical_solution.build_k_mat(
                                grp.get_k_mat()          ,
                                grp.list_of_masses()     ,
                                self.unit_length         ,
                                self.unit_mass           ,
                                self.unit_hbar              )
                eff_k_list, vecs = np.linalg.eigh( eff_k_mat ) 

                effect_k_list, states, eigen_energy, vectors =\
                    analytical_solution.analytical_main(
                        grp.list_of_masses()     ,
                        grp.get_k_mat()          ,
                        rec_states               ,
                        self.unit_length         ,
                        self.unit_mass           ,
                        self.unit_hbar            )
                   
            
            #elif type( self.rec_instruction[ loop_index ] ) == float:
                        #    print("\nFragment",loop_index,"requested for", self.rec_instruction[ loop_index ],"energy limit")
                        #energy_upper_bound = self.rec_instruction[ loop_index ]
                        #effect_k_list, states, eigen_energy, vectors =\
                      #analytical_solution.analytical_main( 
                        #grp.list_of_masses()     ,
                        # grp.get_k_mat()          ,
                        #     energy_upper_bound       ,
                        #     self.unit_length         ,
                        #     self.unit_mass           ,
                        #     self.unit_hbar            )
            else:
                print("MUST SUPPLY NUMBER OF STATES AS INT TYPE FOR RECOUPLE CALCULATION!")
                sys.exit(1)

            self.all_grp_effect_k      += [ effect_k_list ] 
            self.all_grp_states_num    += [ states        ] 
            self.all_grp_eigen_energy  += [ eigen_energy  ]
            self.all_grp_eigen_vector  += [ vectors       ] 
            loop_index += 1

        # Loop Over Pairs of Groups
        # store all Hamiltonian Pairs 

        # This is for debugging...
        from Applications.component_tests.transition_dipole import transition_dipole

        self.all_pair_trans_dipole_mat = []
        for i in range( self.num_fragment ):
            for j in  range( i+1 , self.num_fragment ):
                self.all_pair_trans_dipole_mat += \
                       [ transition_dipole.paired_transition_dipole_mat(\
                                self.all_grp_states_num[i]        ,
                                self.all_grp_eigen_energy[i]      ,
                                self.all_grp_eigen_vector[i]      ,
                                self.all_grp_masses[i]            ,
                                self.all_grp_effect_k[i]          ,
                                self.all_grp_states_num[j]        ,
                                self.all_grp_eigen_energy[j]      ,
                                self.all_grp_eigen_vector[j]      , 
                                self.all_grp_masses[j]            ,
                                self.all_grp_effect_k[j]          ,
                                self.top_level_coupling[ i ][ j ] ,
                                self.unit_length                  ,
                                self.unit_mass                    ,
                                self.unit_hbar                     ) ]
                           
    def get_pair_dipole_mat(self):
        return copy.deepcopy( self.all_pair_trans_dipole_mat )

    def get_unsorted_eigen_energy(self):
        return copy.deepcopy( self.all_grp_eigen_energy )

    def get_sorted_eigen_energy(self):
        # This function is the SAME as get_unsorted_eigen_energy() 
        # b/c they are already sorted but just to be safe.
        return [ sorted( item ) for item in self.all_grp_eigen_energy ] 
    '''\
    def energy_bound_guess(self, rec_states, k_list, hbar):
        # rec_states is the state number requested for this fragment in a recouple calc 
        # guess of energy_upper_bound = hbar * sqrt(k) * rec_states
        k_list = sorted( k_list )
        num_HO     = len( k_list )
        if ( num_HO % 2 == 1 ):
            k = k_list[ (num_HO - 1) // 2 ]  # Median if odd
        else:
            k = k_list[ num_HO // 2 ]        # right median if even
        return hbar * sqrt( k ) * rec_states '''

    def get_hamiltonian_dimension(self):
        dim = 1
        for i in range( self.num_fragment ):
            dim *= len( self.all_grp_eigen_energy[i] )
        return dim










