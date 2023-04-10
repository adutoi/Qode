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
import numpy as np
from qode.coupled_oscillators.c_mat_vec_prod_routine import mat_vec_prod_wrapper
#from qode.coupled_oscillators.c_mat_vec_prod_routine.omp_prod import mat_vec_prod_wrapper



def frag_num_to_coupling_mat_index( num1, num2, num_frag , coupling_mat_list):
    # We have a list of coupling matrices,
    # from Fragment num1 and num2, we can find the corresponding coupling matrix in the long list
    # This function returns the target coupling matrix.
    if num1 > num2: # A guard in case num1 > num2 ( the list is built in such pattern, V12 V13 V14.... V23 V24... Vn-1 n )
        tmp  = num1
        num1 = num2
        num2 = tmp

    index_in_list  = 0
    num_column = num_frag - 2
    for i in range(num1):
        index_in_list += num_column
        num_column    -= 1
    index_in_list = index_in_list + num2 - 1  # num2 starts at 1 but the index is still 0, therefore minus one is needed.
    return coupling_mat_list[ index_in_list ]


def get_total_dimension( num_states_list ):
    # num_states_list is a list of highest state numbers for each fragment, e.g.: [ N1, N2, ... Nn ] 
    # This function returns the Hamiltonian/Vector Dimension.
    total_ld = 1
    for item in num_states_list:
        total_ld *= item
    return total_ld


def get_jumping_distances( num_states_list ):
    # num_states_list is a list of highest state numbers for each fragment, e.g.: [ N1, N2, ... Nn ]
    # This function computes the jumping distances for given i-th fragment, then make a list.
    # For i-th fragment, the jumping distance = N_i+1 * N_i+2 * ... * N_n
    # It returns a list of jumping distances [ jump1, jump2 ... jump_n ]
    num_frag = len( num_states_list )
    jumps = []
    for i in range(num_frag):
        multiplier = 1
        for j in range(i+1,num_frag):
            multiplier *= num_states_list[j]
        jumps += [ multiplier ]
    return jumps 



def compute_loop_limits( num_states_list, position1 , position2 ):
    # num_states_list is a list of highest state numbers for each fragment, e.g.: [ N1, N2, ... Nn ] 
    # It returns the product of N_position1 * ... * N_position2, which computes the loop limits for L1, L2 and L3
    # position2 is included ( not same as python range function )
    multiplier = 1
    for i in range(position1, position2+1):
        multiplier *= num_states_list[i]
    return multiplier



def compute_loop_parameters( couple_num1, couple_num2, num_frag, jumps, num_states_list ):
    # couple_num1 == 0 : No Left Loop
    if couple_num1 == 0:
        left_loop_limit  =  1   # Set as one to allow other loops to run
        left_jump        =  0   # set as zero to eliminate the effect ( to be honest it's meaningless )
    else:
        left_loop_limit  =  compute_loop_limits( num_states_list, 0 , couple_num1 - 1 )
        left_jump        =  jumps[ couple_num1 - 1 ]
    
    # couple_num2 ==  num_frag - 1 : No Right Loop
    if couple_num2 == num_frag - 1:
        right_loop_limit  =  1       # Set as one to allow other loops to run
    else:
        right_loop_limit  =  compute_loop_limits( num_states_list, couple_num2 + 1, num_frag -1 )
    # right_jump is ALWAYS one.


    # couple_num1 and couple_num2 are adjacent : No Middle Loop
    if couple_num1 + 1 == couple_num2:
        mid_loop_limit    =  1        # set as one to allow the rest to run
        mid_jump          =  0        # set as zero to eliminate the effect ( to be honest it's meaningless )
    else:
        mid_loop_limit    =  compute_loop_limits( num_states_list , couple_num1 + 1, couple_num2 - 1 )
        mid_jump          =  jumps[ couple_num2 -1 ]

    return left_loop_limit, mid_loop_limit, right_loop_limit, left_jump, mid_jump






def diag_act_on_vec( eigval_list , orig_vec ):
    # eigval_list         is a Python List of Numpy 1D Arrays
    # orig_vec            is s Numpy Column Matrix   dim = (n, 1)
    total_ld = orig_vec.shape[0]
    num_frag  = len( eigval_list )
    ld_list  = [ item.shape[0] for item in eigval_list ]
    new_vec   = np.matrix( np.zeros((total_ld,1)) )  # Column Vector
    # Pre-Compute Big Jumps of each index
    jumps = get_jumping_distances( ld_list ) # A full list of least jumps if state number changes.

    for frag_i in range(num_frag):
        new_vec += mat_vec_prod_wrapper.diag_act_on_vec(
                                                        eigval_list[ frag_i ],
                                                        ld_list[ frag_i ],    # number of states
                                                        jumps[ frag_i ],
                                                        orig_vec
                                                        )

    
    return new_vec










def coupling_act_on_vec( eigval_list , coupling_mat_list , orig_vec ):
    # eigval_list         is a Python List of Numpy 1D Arrays
    # coupling_mat_list   is a Python List of Numpy 2D Matrices
    # orig_vec            is s Numpy Column Matrix   dim = (n, 1)
    num_frag = len( eigval_list )
    num_states_list = [ item.shape[0] for item in eigval_list ]
    total_ld = get_total_dimension( num_states_list )
    
    #print("HIGHEST STATES  =",num_states_list)
    #print("TOTAL DIMENSION =", total_ld)


    # Get an emtpy new vector for holding results
    # Same dimension identity (n,1) as input orig_vec
    new_vec = np.matrix( np.zeros(( total_ld , 1 )) )  

    # Pre-Compute the Gaps of each fragment state index
    jumps = get_jumping_distances( num_states_list )
    #print("JUMPS =",jumps)


    # Looping over pairs of fragments. Indices are named neq_num1 and neq_num2
    for couple_num1 in range( num_frag ):
        for couple_num2 in range( couple_num1 + 1, num_frag ):
            # Pre-Compute All Parameters for loops:
            # Ma =  couple1_jump
            # Mb =  couple2_jump
            # M1 =  left_jump  
            # M2 =  mid_jump  
            # I1 =  left_loop_limit
            # I2 =  mid_loop_limit 
            # I3 =  right_loop_limit
            
            target_mat         = frag_num_to_coupling_mat_index( couple_num1, couple_num2, num_frag, coupling_mat_list )
            
            
            couple1_jump = jumps[ couple_num1 ]
            couple2_jump = jumps[ couple_num2 ]
            
            left_loop_limit, mid_loop_limit, right_loop_limit, left_jump, mid_jump  =    \
                       compute_loop_parameters( couple_num1, couple_num2, num_frag, jumps, num_states_list )

            right_jump = 1  # Always one, since the smallest jump is just one.


            new_vec += mat_vec_prod_wrapper.coupling_act_on_vec(
                                                                num_states_list[ couple_num1 ],
                                                                num_states_list[ couple_num2 ],
                                                                jumps[ couple_num1 ],
                                                                jumps[ couple_num2 ],
                                                                target_mat,
                                                                left_loop_limit,
                                                                mid_loop_limit,
                                                                right_loop_limit,
                                                                left_jump,
                                                                mid_jump,
                                                                right_jump,
                                                                orig_vec,
                                                                )



    return new_vec






def mat_vec_product_main( eigval_list , coupling_mat_list , orig_vec ):
    # eigval_list         is a Python List of Numpy 1D Arrays
    # coupling_mat_list   is a Python List of Numpy 2D Matrices
    # orig_vec            is s Numpy Column Matrix   dim = (n, 1)
    
    # For now, I will just do a conversion inside this main function.
    # eigval_list = [ np.array( item ) for item in eigen_list ]
    # coupling_mat_list = [ np.matrix( item ) for item in dipole_mat_list ]
    
    # Compute each vector, return a sum
    diag_new_vec     = diag_act_on_vec( eigval_list , orig_vec ) 
    off_diag_new_vec = coupling_act_on_vec( eigval_list , coupling_mat_list , orig_vec )
    
    return diag_new_vec + off_diag_new_vec








