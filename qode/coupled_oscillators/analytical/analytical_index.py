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
import copy

#
# INPUTS: 1) To initialize the class:
#                    number of HOs 
#         2) To Run the generator: 
#                    energy upper bound  
#                    energy function for one configuraton
#                    other dependencies...
#
# OUTPUTS: list of sorted EigenEnergies and corresponding state numbers.
#

#
# USAGE:
#         one_calculation = analytical_index.index_generator( num_oscillators )
#         states, analytical_energy = index_generator( 
#                                                      energy_upper_bound,
#                                                      energy_function   ,    
#                                                      dependencies ...
#                                                     )
#
#



class index_generator(object):
    def __init__( self, num_ho ):
        self.index_list = [ 0 for i in range( num_ho ) ]
        self.num_ho = num_ho
        self.run = True

    def increase_index( self ):
        self.index_list[ 0 ] += 1

    def round_up_index( self ):
        # end the loop by setting self.run = False
        # from left ( or vice versa ) to right, find zeros, 
        # and round up the first non-zero index.
        # if ONLY the LAST index is non-zero, round_up == end the loop
        i = 0
        check_zero = True
        while check_zero and self.index_list[ i ] == 0:
            i += 1
            if i > self.num_ho - 1:
                check_zero = False
        if i >= self.num_ho - 1: # ONLY the LAST IS NON-ZERO
            self.run = False     # index: i == number_HO - 1
        else: 
            self.index_list[ i ] = 0
            self.index_list[ i+1 ] += 1

    def get_index( self ):
        return self.index_list

    def continue_loop_flag( self ):
        return self.run

    def __call__( self, energy_upper_bound, energy_func, *func_dependencies ):
        all_index = []
        eigen_energy = []
        while self.continue_loop_flag() == True:
            energy = energy_func( self.get_index(), *func_dependencies )
            if energy > energy_upper_bound:
                self.round_up_index()
            else:
                all_index += [ copy.deepcopy( self.get_index() ) ]
                eigen_energy += [ energy ]
                self.increase_index()
        return all_index, eigen_energy        




		
if __name__ =='__main__':
    from math import sqrt
    def get_eigen( state, hbar ):
        return (state + 0.5)*hbar
    def get_energy( index_combo, effect_k_list, hbar ):
        energy = 0.0
        for i in range(len(index_combo)):
            energy += sqrt( effect_k_list[i] ) * get_eigen( index_combo[i], hbar )
        return energy
    
    num_ho = 4 
    energy_upper_bound = 5.0
    effect_k_list = [ 1.0 for i in range( num_ho ) ]
    hbar = 1.0

    index1 = index_generator( num_ho )
    all_index, eigen_energy = index1( energy_upper_bound, get_energy, effect_k_list, hbar )
    print("Energy Upper-Bound =",energy_upper_bound)
    for i in range(len(all_index)):
        print("Index: ", all_index[i]," Eigenvalue =", eigen_energy[i] )





