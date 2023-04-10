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
# All the contraction coefficients are multiplied in this code.
from math import pi
from qode.SCF_diag_full.matrix_loop import mat_loop,rep_loop
from qode.SCF_diag_full.matrix_operation import mat_add,mat_multiply,build_mat,build_mat_rep,printline,build_size_mat
from qode.SCF_diag_full import primitive


# The coefficients are held in a same size matrix and then multiplied correspondingly
# to the primitive matrix.
# Then undergo contraction.



######################### SHARED ROUTINE  #########################################
# mat_loop, rep_loop, mat_add, mat_multiply, build_mat, build_mat_rep

def NFactor(alist):
    'Normalization Factor'
    return pow(2.0*alist[3]/pi,0.75)
####################################################################################




# Two types of loops:
# Loop over every element without referring to its index
# Loop over every element by index
# Here the first one is chosen b/c of simplicity




####### CONTRACTION LOOP FOR OVERLAP, KINETIC, AND NUCLEAR ATTRACTION ##############

def nor_cont_coeff(list1,list2):
    'Normal Contraction'
    return NFactor(list1)* NFactor(list2)* list1[4] * list2[4]


def cont_nor_loop(M,ng):
    #Build matrix M to hold numbers.
    N = build_mat(ng)
    #The Loop:
    x = 0
    for m in range(len(ng)):
        for i in range(ng[m]):
            y = 0
            for n in range(len(ng)):
                for j in range(ng[n]):
                    N[m][n] += M[x][y]
                    y += 1
            x += 1
            
    return N

  
def contract_normal(M,ng,input_list):
    'M: Matrix in primitive gaussians; ng: No. of gaussians'
    # INPUT: M and ng; OUTPUT: Contracted Matrix in terms of orbitals
    C = mat_loop(input_list,nor_cont_coeff)
    M = mat_multiply(M,C)
    N = cont_nor_loop(M,ng)
    
    #print('==============')
    #printline(M)
    #print('==============')
    #printline(N)

    return N

###################################################################################




################ CONTRACTION LOOP FOR ELECTRONIC REPULSION #########################

def rep_cont_coeff(p,q,r,s):
    'Repulsion Matrix Contraction'
    return NFactor(p)* NFactor(q)* p[4] * q[4]*NFactor(r)* NFactor(s)* r[4] * s[4]


def cont_rep_loop(M,ng):
    #N = build_mat a different size one
    #The Loop:
    
    size = len(ng)
    N = build_size_mat(size*size)

    x = 0
    for p in range(size):
        for i in range(ng[p]):        
            for q in range(size):
                for j in range(ng[q]):
                    y = 0
                    for r in range(size):
                        for k in range(ng[r]):
                            for s in range(size):
                                for l in range(ng[s]):
                                    N[p*size + q][r*size + s] += M[x][y]
                                    y += 1
                    x += 1

    return N



def contract_repulsion(M,ng,input_list):
    'M: REPULSION Matrix in primitive gaussians; ng: No. of gaussians'
    # INPUT: M and ng; OUTPUT: Contracted Matrix in terms of orbitals
    C = rep_loop(input_list,rep_cont_coeff)
    M = mat_multiply(M,C)
    N = cont_rep_loop(M,ng)

    return N

####################################################################################





####################### Contraction Main Function  #############################

def cont_main(input_list,ngaussians,nulist,num_e,A_num_list):
    S,T,N,V = primitive.pri_main(input_list,nulist,num_e,A_num_list)

    
    ng = ngaussians
    S_c = contract_normal(S,ng,input_list)
    T_c = contract_normal(T,ng,input_list)
    N_c = contract_normal(N,ng,input_list)

    if num_e > 1:
        V_c = contract_repulsion(V,ng,input_list)
    else:
        V_c = V


    return S_c,T_c,N_c,V_c
