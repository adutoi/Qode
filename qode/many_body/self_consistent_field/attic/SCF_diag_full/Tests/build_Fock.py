#    (C) Copyright 2016 Yuhong Liu
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
import sum_repulsion
import numpy as np
import repul_transform
from matrix_operation import printline,printnpline,writenpline,build_size_mat
from matrix_operation import build_same_np_mat,empty_same_np_mat

######### SPIN NOT INCLUDED, ONLY ALPHA MATRIX COMPUTED  #############
######################################################################




#########################  Main Function ###############################
# Build F(-1) = T_on + N_on as the initial guess (simple and reasonable)
# on: OrthoNormal basis
# Diagonalize F0 (Optimize orbitals the first time) to get U1

# First Fock matrix: F0 = T0 + N0 + V(HF)_pq
# Transform F0 to F1 and then attach summed V(HF)_pq parts:
# V(HF)_pq = Sum_i <p(1)|<i(0)|[V(1-P10)]|q(1)>|i(0)>
#   = Sum_i <p(1)|<i(0)|V|q(1)>|i(0)> - <p(1)|<i(0)|V|i(1)>|q(0)>


def Fock_main(T_on,N_on,V_on,num_e,multi):
    # Initial Guess: F(-1) = T_on + N_on
    F_n = T_on + N_on
    E0,U0 = np.linalg.eig(F_n)
    
    # Transform into MOLECULAR ORBITAL basis
    T0 = U0.H * T_on * U0
    N0 = U0.H * N_on * U0
    if num_e > 1:
        V0 = repul_transform.transform_V(V_on,U0)
    else:
        V0 = V_on
        
    # BUILD FIRST FOCK
    if num_e > 1:    
        vHF = sum_repulsion.sum_repul(V0,num_e,multi)
        F0 = T0 + N0 + vHF
    else:
        F0 = T0 + N0
    #writenpline(V1,"transformed_V.txt")
    #print('F0:')
    #printnpline(F0)
    
    return F0,T0,N0,V0,E0
