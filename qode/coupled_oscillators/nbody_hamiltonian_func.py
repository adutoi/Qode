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
################## nbody Hamiltonian Individual Element Generator  ####################

# INPUT:   A list of Fragment eigen energies ( 2D Array )
#          A list of two-body coupling matrices  ( list of 2D arrays )
#          x and y position in the total Hamiltonian ( never built explicitly )     
#
# OUTPUT:  n-body Hamiltonian value at location (x, y)
#
#
from . import nbody_position_to_frag_index
import numpy as np


def get_delta_flags( x_index_list, y_index_list ):
    # x_index_list, y_index_list are Fragment index numbers 
    # return zeroflag = 1  this element is zero due to zero delta function(s)
    #                      IF NOT DIAGONAL: diag_flag = False 
    #                          n body,  (n-2) equal indecies. two +/- one indecies NOT satisfied.
    #                          
    #        zeroflag = 0  1) Diagonal is non-zero  return num_equ = n, diagflag = True
    #                      2) n body,  (n-2) equal indecies. two +/- 1  indecies SATISFIED.
    #                                 return two locations of +/- 1 delta, diagflag = False
    #
    lda = len( x_index_list )
    non_zero_delta = 0
    non_zero_location = []

    for i in range(lda):
        diff = x_index_list[i] - y_index_list[i]
        if diff == 0:
            pass 
        else:
            non_zero_delta += 1
            non_zero_location += [ i ]

    if non_zero_delta == 0:
        zeroflag = False 
        diagflag = True
    elif non_zero_delta == 2:
        zeroflag = False
        diagflag = False
    else:
        zeroflag = True
        diagflag = True

    return zeroflag, diagflag, non_zero_location 


def frag_num_to_coupling_mat_index( num1, num2, num_frag ):
    if num1 > num2:
        tmp  = num1
        num1 = num2
        num2 = tmp

    index_sum  = 0
    num_column = num_frag - 2
    for i in range(num1):
        index_sum += num_column
        num_column -= 1
    index_sum = index_sum + num2 - 1  # num2 starts at 1 but the index is still 0, therefore minus one is needed.
    return index_sum




def get_hamiltonian_value( eigen_list , coupling_mat_list , x_position , y_position ):
    # CAUTION: x_position denotes LEADING DIMENSION.
    # eigen_list          is a Python List of Numpy 1D Arrays
    # coupling_mat_list   is a Python List of Numpy 2D Matrices
    size_list    = [ item.shape[0] for item in eigen_list ]
    x_index_list = nbody_position_to_frag_index.position_to_index( x_position , size_list ) 
    y_index_list = nbody_position_to_frag_index.position_to_index( y_position , size_list )

    zeroflag , diagflag, non_zero_location = get_delta_flags( x_index_list , y_index_list )
    if zeroflag:
        return 0.0
    else:
        if diagflag:
            h_sum = 0.0
            for i in range( len( size_list ) ):
                h_sum += eigen_list[i][ x_index_list[i] ] 
                # eigen_list ==> [ [h11, h12, ... , h1m] , [h21, h22, ... , h2n] , ... , [hn1, hn2, ... , hnp] ]
                # [i][ x_index_list[i] ] retrieves current fragment i and its current index number x_index_list[i]  
                # h_ij ==> eigen_list[i][ x_index_list[i] ] 
            return h_sum
        else:
            n1,n2 = non_zero_location # n1 and n2 are the n1-th and n2-th fragment

            shift = len( eigen_list[ n2 ] )
            # < i1 |< j1 | V_{n1,n2} | i2 >| j2 > 
            i1_n1 = x_index_list[ n1 ]  # < i1 |
            i2_n1 = y_index_list[ n1 ]  # | i2 >
            j1_n2 = x_index_list[ n2 ]  # < j1 | 
            j2_n2 = y_index_list[ n2 ]  # | j2 > 
            
         
            # Find V_{n1,n2} to retrieve < i1 |< j1 | V_{n1,n2} | i2 >| j2 > 
            # coupling_mat_list = [ V12, V13, ... , V1n, V23 ,..., V2n , V34 ... V3n , ... ... V(n-2)(n-1), 
            #                       V(n-2)(n) , V(n-1)(n) ]
            couple_mat_index = frag_num_to_coupling_mat_index( n1, n2, len( size_list ) )            
            target_mat = coupling_mat_list[ couple_mat_index ] 

            return target_mat[i1_n1 * shift + j1_n2 , i2_n1 * shift + j2_n2]



def rec_obj_build_full_hamiltonian( rec_obj ):
    eigen_list   = [ np.array(item) for item in rec_obj.get_sorted_eigen_energy() ]
    coupling_mat_list = [ np.array(item) for item in rec_obj.get_pair_dipole_mat() ]
    size_list    = [ item.shape[0] for item in eigen_list ]
    
    ld = 1
    for i in size_list:
        ld *= i
    
    H = np.matrix( np.zeros((ld , ld)) )
    
    x = 0
    while x < ld:
        y = 0
        while y < ld:
            H[x,y] = get_hamiltonian_value( eigen_list , coupling_mat_list , x , y )
            y += 1
        x += 1

    return H



if __name__ == "__main__":
    n = 10
    num_frag = 0
                  
    for i in range(n):
        for j in range(i+1,n):
            index_in_list = frag_num_to_coupling_mat_index( i, j, n )
            print(i,j,"(",index_in_list,")",end='\t')  
            num_frag += 1
        print(" ")
    print("num_frag =",num_frag)




