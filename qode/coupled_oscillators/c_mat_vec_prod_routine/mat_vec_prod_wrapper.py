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
from ctypes import c_int, c_double, c_void_p, cdll, pointer
import numpy as np



def diag_act_on_vec(
                    eigval_list,     # Numpy 1D array : EigenValue list for One Fragment (1D-Array).
                    num_state,       # integer type : Number of states for this Fragment.
                    jumps,           # integer type : Least jump if current fragment state-number plus one
                    orig_vec         # Numpy Column Matrix
                    ):
    #
    # This is a wrapper for function "diag_act_on_vec" from "libMatVecProd.so" library.
    # Its purpose it to convert Python/Numpy Data Types to C Data Types and Pass them to C function.
    #
    # This C function has such input interface:
    #
    # int diag_act_on_vec(
    #                     int      total_ld   /* Leading Dimension of the vector ( H Matrix as well ) */
    #                     double  *eigval_list     ,  /* ONE of the eigenvalue list (1D array) */
    #                     int       num_state       ,  /* number of states of current fragment  */
    #                     int       jumps           ,  /* Least jump if current fragment state-number plus one */
    #                     double  *orig_vec        ,
    #                     double  *new_vec )
    #
    
    # Convert Numpy to C types:
    #
    total_ld    = orig_vec.shape[0]
    total_ld_c  = c_int( total_ld )  # Convert LD to C-type int
    num_state_c = c_int( num_state )
    jumps_c     = c_int( jumps )
    
    
    eigval_list_c  = c_void_p( eigval_list.ctypes.data )  # a pointer to numpy array casted to python ctypes ( you can pass this to C function directly )
    orig_vec_c     = c_void_p( orig_vec.ctypes.data )   # Stored Continuously. In C function, treated still Continuously. Numpy Column Matrix is OK to pass.     
    new_vec        = np.matrix( np.zeros((total_ld,1)) )
    new_vec_c      = c_void_p( new_vec.ctypes.data )

    
    
    
    # Calling C function diag_act_on_vec()
    #
    lib = "./libMatVecProd.so"
    product = cdll.LoadLibrary(lib)
    
    product.diag_act_on_vec(
                            total_ld_c,
                            eigval_list_c,
                            num_state_c,
                            jumps_c,
                            orig_vec_c,
                            new_vec_c
                            )
    

    
    return new_vec 



def coupling_act_on_vec(
                        num_state1,
                        num_state2,
                        num1_jump,
                        num2_jump,
                        coupling_mat,  # Only this is an 2D-Matrix, the rest are of integer type
                        left_limit,
                        mid_limit,
                        right_limit,
                        left_jump,
                        mid_jump,
                        right_jump,
                        orig_vec       # Only this is an Column Numpy Matrix, the rest are of integer type
                        ):
    #
    # This is a wrapper for function "coupling_act_on_vec" from "libMatVecProd.so" library.
    # Its purpose it to convert Python/Numpy Data Types to C Data Types and Pass them to C function.
    #
    # This C function has such input interface:
    #
    #
    # int coupling_act_on_vec(
    #                        int num_state1,
    #                        int num_state2,         /* shift in coupling mat = num_state2 */
    #                        int num1_jump,
    #                        int num2_jump,
    #                        int total_ld,           /* Leading Dimension of the vector ( H Matrix as well ) */
    #                        double *coupling_mat,   /* coupling matrix for this pair                        */
    #                        int left_limit,         /* loop limits */
    #                        int mid_limit,          /* loop limits */
    #                        int right_limit,        /* loop limits */
    #                        int left_jump,
    #                        int mid_jump,
    #                        int right_jump,
    #                        double *orig_vec,
    #                        double *new_vec
    #                        )
    #
    #
    
    # Convert Numpy to C types:
    #
    total_ld    = orig_vec.shape[0]
    
    num_state1_c    =  c_int( num_state1 )
    num_state2_c    =  c_int( num_state2 )
    num1_jump_c     =  c_int( num1_jump )
    num2_jump_c     =  c_int( num2_jump )
    total_ld_c      =  c_int( total_ld )
    
    left_limit_c    =  c_int( left_limit )
    mid_limit_c     =  c_int( mid_limit )
    right_limit_c   =  c_int( right_limit )
    
    left_jump_c     =  c_int( left_jump )
    mid_jump_c      =  c_int( mid_jump )
    right_jump_c    =  c_int( right_jump )
    

    orig_vec_c     = c_void_p( orig_vec.ctypes.data )   # Stored Continuously. In C function, treated still Continuously. Numpy Column Matrix is OK to pass.
    new_vec        = np.matrix( np.zeros((total_ld,1)) )
    new_vec_c      = c_void_p( new_vec.ctypes.data )
    
    coupling_mat_c = c_void_p( coupling_mat.ctypes.data )
    
    # Calling C function coupling_act_on_vec()
    #
    lib = "./libMatVecProd.so"
    product = cdll.LoadLibrary(lib)
    product.coupling_act_on_vec(
                                num_state1_c,
                                num_state2_c,
                                num1_jump_c,
                                num2_jump_c,
                                total_ld_c,
                                coupling_mat_c,
                                left_limit_c,
                                mid_limit_c,
                                right_limit_c,
                                left_jump_c,
                                mid_jump_c,
                                right_jump_c,
                                orig_vec_c,
                                new_vec_c
                                )

    
    
    return new_vec 


