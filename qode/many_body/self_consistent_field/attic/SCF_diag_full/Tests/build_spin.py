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
from matrix_operation import build_size_mat,printline,printnpline,writenpline
import numpy as np

# This code builds big matrices with spins taken into considerations.
# Format of spins: (T and N)
#  alpha alpha, alpha beta
#  beta  alpha, beta  beta

def build_loop(M,N):
    'M:small mat; N:new big mat'
    size = len(M)
    size_n = len(N)
    
    for i in range(size_n):
        for j in range(size_n):
            if (i<= size - 1 and j<= size - 1):
                N[i,j]=M[i,j]
            elif (i>= size and j>= size):
                N[i,j]=M[i-size,j-size]
                
    return N

def build_rep_loop(M,N):
    'M:small mat; N:new big mat'
    size = len(M)
    size_n = len(N)
    
    for i in range(size_n):
        for j in range(size_n):
            if (i< size and j < size):
                N[i,j]=M[i,j]
            elif (size <= i < 2*size  and size <= j < 2*size):
                N[i,j]=M[i - size,j - size]
            elif (2*size <= i < 3*size  and 2*size <= j < 3*size):
                N[i,j]=M[i - 2*size,j - 2*size]
            elif (i >= 3*size and j >= 3*size):
                N[i,j]=M[i - 3*size,j - 3*size]
    
    return N

def spin_mat(M):
    size = len(M)
    N = build_size_mat(size*2)
    N = np.matrix(N)
    N = build_loop(M,N)
    return N

def spin_rep_mat(M):
    size = len(M)
    N = build_size_mat(size*4)
    N = np.matrix(N)
    N = build_rep_loop(M,N)
    return N

#================= Main Function ============================

def spin_main(S_cnt,Uni_t,T_t,N_t,V_t):
    # If you diagonalize S_sp, and you NORMALIZE it!
    # You will get the same Unitary Matrix as the doubled U_sp.
    
    S = np.matrix(S_cnt)
    
    S_sp = spin_mat(S)
    U_sp = spin_mat(Uni_t)
    T_sp = spin_mat(T_t)
    N_sp = spin_mat(N_t)
    V_sp = spin_rep_mat(V_t)
    
    
    #writenpline(V_t,'V-transed.txt')
    #writenpline(V_sp,'V-spin.txt')
    
    return S_sp,U_sp,T_sp,N_sp,V_sp









