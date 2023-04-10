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
import numpy as np
from qode.SCF_Code.matrix_operation import build_size_mat


def idx_M(size,step_num,idx1,idx2,idx3,idx4,idx):
    if step_num == 1:
        x = idx1*size + idx2
        y = idx3*size + idx
    elif step_num == 2:
        x = idx1*size + idx2
        y = idx*size + idx4
    elif step_num == 3:
        x = idx1*size + idx
        y = idx3*size + idx4
    elif step_num == 4:
        x = idx*size + idx2
        y = idx3*size + idx4
        
    return x,y

def idx_U(size,step_num,idx1,idx2,idx3,idx4,idx):
    x = idx
    if step_num == 1:
        y = idx4
    elif step_num == 2:
        y = idx3
    elif step_num == 3:
        y = idx2
    elif step_num == 4:
        y = idx1
        
    return x,y    

def trans_loop(M,U,step_num):
    size = len(U)
    N = build_size_mat(size*size)
    N = np.matrix(N)
    
    for idx1 in range(size):
        for idx2 in range(size):
            for idx3 in range(size):
                for idx4 in range(size):
                    for idx in range(size):
                        N[idx1*size+idx2,idx3*size+idx4] += M[idx_M(size,step_num,idx1,idx2,idx3,idx4,idx)]\
                                                               * U[idx_U(size,step_num,idx1,idx2,idx3,idx4,idx)]

    return N

def transform_V(M,U):
    # No Spin Considered Here! ONLY SPATIAL INTEGRALS.
    'INPUT:Unitary M, Repulsion M; OUTPUT:transformed Rep M'
    #Loop over i,j,k,l and in each one assume p is the l for big mat, and p is p for Upl.
    #So the loop is i,j,k,p -- p,l
    #Then sum
    # m,n,o,p -- p,l = m,n,o,l
    # m,n,o,l -- o,k = m,n,k,l
    # m,n,k,l -- n,j = m,j,k,l
    # m,j,k,l -- m,i = i,j,k,l
  
    size = len(U)
    M = np.matrix(M)
    
    A1 = trans_loop(M,U,1)
    A2 = trans_loop(A1,U,2)
    A3 = trans_loop(A2,U,3)
    A4 = trans_loop(A3,U,4)                                      
    return A4

