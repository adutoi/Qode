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
from math import sqrt
from qode.SCF_Code.matrix_operation import build_size_mat,printnpline



def sum_clb_exchng(p,q,V1,alpha_e,beta_e,n):
    "Sum_i <p(1)|<i(0)|V|q(1)>|i(0)> - <p(1)|<i(0)|V|i(1)>|q(0)>"
    sum_rep = 0.0
    
    
    for j in range(alpha_e):
        sum_rep += ( V1[p*n + j, q*n + j] - V1[p*n + j, j*n + q] )
    if (beta_e > 0):
        for j in range(n//2, n//2 + beta_e):
            sum_rep += ( V1[p*n + j, q*n + j] - V1[p*n + j, j*n + q] )
            
    return sum_rep


##################### MAIN FUNCTION  ####################################
##############  ONLY SUM OVER ALPHA MATRIX  #############################

def sum_repul(V1,alpha_e,beta_e):

    size2 = len(V1)
    n = int(sqrt(size2)) # a+b

    "Build vHF repulsion matrix"
    vHF = np.zeros((n,n))
    vHF = np.matrix(vHF)

    for p in range(n):
        for q in range(n):
            vHF[p,q] = sum_clb_exchng(p,q,V1,alpha_e,beta_e,n)

    return vHF
