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
import copy
from math import sqrt

#
#  THIS MODULE IS FOR BUILDING PAIR-WISE HARMONIC OSCILLATOR CI MATRIX 
#
#  THE PAIR-WISE HAMILTONIAN IS BUILT IN EACH GROUP'S ENERGY LOW TO HIGH ORDER.
# 
#


#
#  Numpy Array/Matrix to Python List function
#
def np2py( a_numpy_matrix ):
    return copy.deepcopy( a_numpy_matrix.tolist() ) 



#
#  SORTING FUNCTIONS FOR ENERGY
#
def sort_by_energy( list_of_index, list_of_energy ):
    # This function sorts the eigen energies and those corresponding
    # oscillator state numbers. 
    # INPUTS :  all states in a list of list , all eigen energies
    # OUTPUTS:  deepcopy of both sorted lists
    sorted_index   =  copy.deepcopy( list_of_index  )
    sorted_energy  =  copy.deepcopy( list_of_energy )
    if len( sorted_index ) != len( sorted_energy ):
        print("SORT_BY_ENERGY FUNCTION GOT INVALID INPUTS!")
        sys.exit(1)
    else:
        total_pairs = len( sorted_index ) - 1 
        for i in range( total_pairs ):
            for j in range( total_pairs - i ):
                #print("compare",j,"=",sorted_energy[j],"and",j+1,"=",sorted_energy[j+1]) 
                if sorted_energy[j] > sorted_energy[j+1]:
                    tmp_E               =  sorted_energy[j] 
                    sorted_energy[j]    =  sorted_energy[j+1]
                    sorted_energy[j+1]  =  tmp_E
                    tmp_i               =  sorted_index[j]
                    sorted_index[j]     =  sorted_index[j+1]
                    sorted_index[j+1]   =  tmp_i   

    return sorted_index, sorted_energy

#  Assertation of deltas for pairs of state numbers.
#
#
#
#
#

def delta_assert( list_of_states1 , list_of_states2 , num_diff_allowed ):
    # One-Dimensional list of states for 
    #    bra:  list_of_states1
    #    ket:  list_of_states2
    # Number of different states Allowed.
    if len(list_of_states1) == len(list_of_states2):
        num_diff = 0
        i = 0                 # index
        while (i < len(list_of_states1) and num_diff <= num_diff_allowed ):
            if list_of_states1[i] != list_of_states2[i]:
                num_diff += 1
            i += 1
        if num_diff == num_diff_allowed:
            return True
        else:
            return False

    else:
        print("Number Of States Not Alignned. Terrible Error!")
        sys.exit(1)

#
#
# Get the SINGLE DIFFERENT STATE NUMBER ( must differ by one state or exit )
#
#
#
def one_diff_state_index( list_of_states1 , list_of_states2 ):
    num_diff_allowed = 1
    if delta_assert( list_of_states1 , list_of_states2 , num_diff_allowed ):
        return_index = 0
        seek = True
        while return_index < len( list_of_states1 ) and seek:
            if list_of_states1[ return_index ] == list_of_states2[ return_index ]:
                return_index += 1
            else:
                seek = False
        if return_index < len( list_of_states1 ):
            return return_index
        else:
            print("Terrible Error, main function didn't phase out same states ZERO TERMS")
    else:
        print("Terrible Error, main function didn't phase out different state ZERO TERMS")
        sys.exit(1)


def dipole_delta_assert( first_state_num, second_state_num ):
    if first_state_num == second_state_num - 1 \
        or first_state_num == second_state_num + 1:
        return True 
    else:
        return False



#
#  INTEGRAL ENGINES:
#
#  Non-rotated unitless coordinate transition dipole integral.
#  Rotated transition dipole integral for each pair of oscillators. 
#  Rotated Pair-Wise transition dipole integral for each state configuraton.
#
#
#
def core_dipole_integral( first_state_num, second_state_num ):
    # Transition dipole integral in unitless original coordinates (NOT ROTATTED)
    # integral = <i''|x|i'> =  sqrt( i' ) delta_i'',i'-1 + sqrt( i'+1 ) delta_i'',i'+1 
    if first_state_num == second_state_num - 1:
        return pow( second_state_num , 0.5 ) 
    elif first_state_num == second_state_num + 1:
        return pow( second_state_num + 1  , 0.5 )
    else:
        return 0.0





def dipole_integral_engine(     bra_first_states     ,
                                ket_first_states     ,
                                first_vectors        ,
                                bra_second_states    ,
                                ket_second_states    ,
                                second_vectors       ,
                                first_mass_list      ,
                                second_mass_list     ,
                                first_rotated_k_list ,
                                second_rotated_k_list,
                                unit_length          , 
                                unit_mass            , 
                                unit_hbar           ):
    "LOOPING OVER OSCILLATOR PAIRS"
    # states have to be different for ONLY ONE PAIR for each.
    num_diff_allowed = 1 
    if delta_assert( bra_first_states, ket_first_states, num_diff_allowed ) \
           and delta_assert( bra_second_states , ket_second_states, num_diff_allowed ):
        first_diff_index  = one_diff_state_index( bra_first_states  , ket_first_states  ) 
        second_diff_index = one_diff_state_index( bra_second_states , ket_second_states )
        if dipole_delta_assert( bra_first_states[ first_diff_index ] ,             \
                                ket_first_states[ first_diff_index ] )             \
               and dipole_delta_assert( bra_second_states[ second_diff_index ] ,   \
                                        ket_second_states[ second_diff_index ] ): 
            lda_first_ho  = len( bra_first_states )
            lda_second_ho = len( bra_second_states )
            integral      = 0.0
            for i in range( lda_first_ho ):
                for j in range( lda_second_ho ):
                    integral += 0.5 * \
                                pow( unit_mass * unit_length**2, 2 ) / (unit_hbar**2 * sqrt( first_mass_list[i] * second_mass_list[j] ) )  * \
                                first_vectors[ i ][ first_diff_index ]                            * \
                                sqrt( 1.0 / sqrt( first_rotated_k_list[ first_diff_index ] ) ) *\
                                core_dipole_integral( bra_first_states[ first_diff_index ] ,    \
                                                          ket_first_states[ first_diff_index ] )  * \
                                second_vectors[ j ][ second_diff_index ]                          * \
                                sqrt( 1.0 / sqrt( second_rotated_k_list[ second_diff_index ] ) ) *\
                                core_dipole_integral( bra_second_states[ second_diff_index ] ,  \
                                                          ket_second_states[ second_diff_index ] )
            return integral

        else:
            return 0.0

    else:
        return 0.0 





#
#
#
#
#  MODULE MAIN FUNCTION
#
#
#
#
def pair_dipole_matrix_main( 
                           first_index           , 
                           first_energy          , 
                           first_vectors         ,
                           second_index          ,
                           second_energy         ,
                           second_vectors        ,
                           pair_coupling_const   ,
                           first_mass_list       ,
                           second_mass_list      ,
                           first_rotated_k_list  ,
                           second_rotated_k_list ,
                           unit_length           , 
                           unit_mass             , 
                           unit_hbar
                          ):
    # vectors come sorted due to the use of numpy.linalg.eigh function.
    # Index and energy are in the order of looping indecies but not in 
    # the order of energy, therefore at first sorting is needed.

    first_vectors  = np2py( first_vectors  )
    second_vectors = np2py( second_vectors )


    sorted_index1, sorted_energy1 = sort_by_energy( first_index, first_energy )
    sorted_index2, sorted_energy2 = sort_by_energy( second_index, second_energy )
    

    first_lda  = len( first_index )
    second_lda = len( second_index )
    pair_H_lda = first_lda * second_lda
    #print("\n\n")
    #print("##########################################################")
    #print("FROM PAIR TRANSITION DIPOLE MATRIX MODULE")
    #print( "LD1       =", first_lda )
    #print( "LD2       =", second_lda )
    #print( "LD_TOTAL  =", pair_H_lda )

    dipole_mat = [ [ 0.0 for i in range( pair_H_lda ) ] for j in range( pair_H_lda ) ]

    # Filling in off-diagonal elements
    # number of states must differ = 1

    # <bra_first_state| <bra_second_state| V |ket_first_state> |ket_second_state>
    #print("Building Pair-wise Transition Dipole Matrix:")
    x_index = 0
    for i in range( first_lda ):
        bra_first_states  = sorted_index1[ i ]
        for j in range( second_lda ):
            bra_second_states = sorted_index2[ j ]
            y_index = 0
            for k in range( first_lda ):
                ket_first_states = sorted_index1[ k ]
                for l in range( second_lda ):
                    ket_second_states = sorted_index2[ l ]
                    dipole_mat[ x_index ][ y_index ] += \
                           pair_coupling_const *         \
                           dipole_integral_engine( 
                                                 bra_first_states     ,  
                                                 ket_first_states     , 
                                                 first_vectors        ,
                                                 bra_second_states    ,
                                                 ket_second_states    ,
                                                 second_vectors       ,
                                                 first_mass_list      ,
                                                 second_mass_list     ,  
                                                 first_rotated_k_list ,
                                                 second_rotated_k_list,
                                                 unit_length          ,  
                                                 unit_mass            , 
                                                 unit_hbar            )             
                    y_index += 1
            x_index += 1 

    prefactor = unit_hbar ** 2 / ( unit_mass * unit_length **2 )
    if prefactor != 1.0:
        for i in range( len(dipole_mat) ):
            for j in range( len(dipole_mat) ):
                dipole_mat[i][j] *= prefactor
        print("Energy Prefactor NOT Equal to ONE")
    else:
        pass
        #print("Energy Prefactor = 1.0")

    #print("Exiting Transition Dipole Matrix Module")
    #print("##########################################################")
    #print("\n\n")
    return dipole_mat 





