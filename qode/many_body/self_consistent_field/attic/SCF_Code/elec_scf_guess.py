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
import numpy as np
from qode.SCF_Code import build_spin
from qode.SCF_Code.matrix_operation import printnpline
from qode.SCF_Code import eigensort


def cp_col(M1, col1, M2, ol2):
    'copy column 2 in M2 to column 1 in M1'
    for i in range(len(M1)):
        M1[i,col1] = M2[i,col2]

def guess_main(T_sp, N_sp, alpha_e, beta_e):
    # Initial Guess: F(-1) = T_sp + N_sp
    # diagonalize F(-1) to get E0, U0
    n = len(T_sp)
    E,U = np.linalg.eig(T_sp + N_sp)
    Ia,Ib = eigensort.sort_main(E)
    
    guess = np.matrix( np.zeros((n, alpha_e + beta_e)) )
    for i in range(alpha_e):
        for j in range(n):
            guess[j,i] = U[j,Ia[i]]
    for i in range(beta_e):
        for j in range(n):
            guess[j,alpha_e + i] = U[j,Ib[i]]


    return guess, E, U


