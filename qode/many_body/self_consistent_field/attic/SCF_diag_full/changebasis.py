#    (C) Copyright 2016 Yuhong Liu
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
# This code diagonalize the overlap matrix
import numpy as np
from qode.SCF_diag_full import repul_transform
from math import sqrt
from qode.SCF_diag_full.matrix_operation import printline,build_same_np_mat,build_size_mat,printnpline

def diagonalize(M):
    'diagonalize a matrix, returns E and Unitary Matrix'
    # M:Python Matrix; E,U: numpy Matrices
    M_np = np.matrix(M)
    E,U = np.linalg.eigh(M_np)
    return E,U



def normalize(E,U):
    'INPUT:Eigen values, Unitary Matrix; OUTPUT:Normalized U'
    # E,U: numpy Matrices
    size = len(U)
    for i in range(size):
        for j in range(size):
            U[j,i] = U[j,i]/sqrt(E[i])
    return U



def test_orthonormal(U_on,S_c):
    Test = U_on.H * S_c * U_on
    return Test


    
def transform(M_c,U_on):
    Uc = U_on.H # Hermitian Conjugate
    M_c = np.matrix(M_c) 
    M_on = Uc * M_c * U_on
    return M_on




    
####################### Change of Basis Main Fuction #########################

def chbasis_main(S_c,T_c,N_c,V_c,num_e):
    E,U = diagonalize(S_c)
    U_on = normalize(E,U)

    T_on = transform(T_c,U_on)
    N_on = transform(N_c,U_on)
    if num_e > 1:
        V_on = repul_transform.transform_V(V_c,U_on)
    else:
        V_on = V_c
    
    #Test = test_orthonormal(U_on,S_c)
    #printnpline(Test)
    
    return E,U_on,T_on,N_on,V_on
