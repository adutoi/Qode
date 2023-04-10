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
from qode.SCF_Code import repul_transform
from qode.SCF_Code.matrix_operation import cp_arbit_np_mat

def HF_ENERGY(U1,T1,N1,V1, alpha_e,beta_e,NRE):
    "EHF = Sum_i [<i1|T|i1> + <i1|N|i1> ] + 1/2 * Sum_i,j <i1|<j0|V(1-P10)|i1>|j0>"
    #Initiate looping indicies
    n = len(T1)//2
    U = cp_arbit_np_mat(U1)
    #First transform all matrices into latest basis, then calculate E
    T = cp_arbit_np_mat(T1)
    N = cp_arbit_np_mat(N1)
    
    
    T = U.H * T * U
    N = U.H * N * U
    if (alpha_e + beta_e) > 1:
        V = cp_arbit_np_mat(V1)
        V = repul_transform.transform_V(V,U)
    
    "Kinetic Energy, Nuclear Attraction E., Repulsion E."
    KE = 0.0
    NAE = 0.0
    ERE = 0.0
    # Kinetic and Attraction Summing Loops ========
    
    for i in range(alpha_e):
        KE +=  T[i,i]
        NAE += N[i,i]
    
    if beta_e > 0:
        for i in range(n,n + beta_e):
            KE +=  T[i,i]
            NAE += N[i,i]


    # Repulsion Summing Loops ======================
    # alpha-alpha repulsion
    if (alpha_e + beta_e > 1):
        for i in range(alpha_e):
            for j in range(alpha_e):
                ERE += ( V[i*2*n+j,i*2*n+j]-V[i*2*n+j,j*2*n+i] )

        if (beta_e >0):
            # alpha-beta repulsion
            for i in range(alpha_e):
                for j in range(n, n + beta_e):                    
                    ERE += ( V[i*2*n+j,i*2*n+j]-V[i*2*n+j,j*2*n+i] )
            # beta-alpha repulsion 
            for j in range(n, n + beta_e):
                for i in range(alpha_e):
                    ERE += ( V[j*2*n+i,j*2*n+i]-V[j*2*n+i,i*2*n+j] )
            # beta-beta repulsion
            for i in range(n, n + beta_e):
                for j in range(n, n + beta_e):
                    ERE += ( V[i*2*n+j,i*2*n+j]-V[i*2*n+j,j*2*n+i] )

    ERE = ERE/2.0
    
    EHF = KE + NAE + NRE + ERE 
    return EHF,KE,NAE,ERE

