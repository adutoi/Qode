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
# CONTRACTION Coefficients not included in this code.
import sys
import numpy as np
from math import exp,pi,sqrt,erf
from qode.SCF_diag_full.matrix_loop import mat_loop,rep_loop
from qode.SCF_diag_full.matrix_operation import mat_add,mat_multiply,build_mat,printline
'''
_list_ is a datatype defined as:
[ [x, y, z, exponential, contraction coefficient],[...],[...],[...],... ]
'''

####################### Shared Routine: ####################################
# mat_loop, rep_loop, mat_multiply(M,C), mat_add



####################### Shared Functions ###################################

def distance2(t1,t2):   
    # INPUT: _list_ elements t1 and t2
    # OUTPUT: First three x,y,z coordinates to give distance squared
    R2 = pow(t1[0]-t2[0],2) + pow(t1[1]-t2[1],2) + pow(t1[2]-t2[2],2)
    return R2/(0.52917721092*0.52917721092)

def centerf(p,q,x,y):
    return (p*x+q*y)/(p+q)


def center_coor(l1,l2):
    p = l1[3]
    q = l2[3]
    x = centerf(p,q,l1[0],l2[0])
    y = centerf(p,q,l1[1],l2[1])
    z = centerf(p,q,l1[2],l2[2])
    return [x,y,z,0.0,0.0]


def F0(T):
    # For atoms with only s orbitals (H,He)
    if T <= 15.0:
        t = 1.0
        s = t*exp(-T)
        n = 1       
        while n<=100:
            t *= 2.0*T/(2.0*n+1.0)
            s += t*exp(-T)
            n += 1
        #prt += str(round(T,2))+'\t'
        #correct = erf(sqrt(T))/(2*sqrt(T/pi))
        #prt += str(fabs((s-correct)/correct))+'\n'
        return s
    else:
        print('T out of range!')
        sys.exit(1)
##############################################################################

################## OVERLAP INTEGRAL  #########################################
def S_pri_int(list1,list2):
    '''
    INPUT:
        _list_ elements: list1, list2
    OUTPUT:
        OVERLAP INTEGRAL
    '''
        
    # Compute OVERLAP Matrix.
    #
    # Eqn: INT_r c_i * N_i*exp[-p(r-P)^2] * c_j * N_j * exp[-q(r-Q)^2] * dr
    # = c_i * c_j * N_i * N_j * (pi/a)^(3/2) * exp[-miu*R^2]
    #
    # where N_i = (2*p/pi)^(3/4) and N_j = (2*q/pi)^(3/4)
    # where a = p + q, R =|P - Q|, miu = pq/(p + q)

    p = list1[3]
    q = list2[3]
    a = p + q
    mu = p * q / a
    R2 = distance2(list1,list2)
    return pow(pi/a,1.5)*exp(-mu*R2)

def overlap_mat(input_list):
    'INPUT: _list_, OUTPUT: OVERLAP MATRIX'
    S = mat_loop(input_list,S_pri_int)
    return S

##############################################################################


################## KINETIC INTEGRAL  #########################################
def K_pri_int(list1,list2):
    '''
    INPUT:
        _list_ elements: list1, list2
    OUTPUT:
        KINETIC ENERGY INTEGRAL
    '''
    # Compute T Matrix.
    #
    # Eqn: INT_r c_i * N_i*exp[-p(r-P)^2] T c_j * N_j * exp[-q(r-Q)^2] * dr
    # = c_i * c_j * N_i * N_j * (1/2)*mu*(6-4muR^2)*(pi/a)^(3/2) * exp[-mu*R^2]
    #
    # where N_i = (2*p/pi)^(3/4) and N_j = (2*q/pi)^(3/4)
    # where a = p + q, R =|P - Q|, mu = pq/(p + q)

    p = list1[3]
    q = list2[3]
    a = p + q
    mu = p * q / a
    R2 = distance2(list1,list2)
    return 0.5*mu*(6.0-4.0*mu*R2)*pow(pi/a,1.5)*exp(-mu*R2)

def kin_mat(input_list):       
    '''
    INPUT: _list_: (x, y, z, exponential, contraction coefficient)
    OUTPUT: Kinetic Energy Matrix '''
    T = mat_loop(input_list,K_pri_int)
    return T

##############################################################################


########### NUCLEAR ATTRACTION INTEGRAL  #####################################


def C_N_dis(list1,list2,nucleus):
    N = nucleus
    C = center_coor(list1,list2)
    d_A = sqrt(pow(C[0]-N[0],2) + pow(C[1]-N[1],2) + pow(C[2]-N[2],2))
    return d_A/0.52917721092

def U_pri_coeff(list1,list2):
    p = list1[3]
    q = list2[3]
    a = p + q

    d2 = distance2(list1,list2)
    return exp(-p*q*d2/a)*pow(pi/a,1.5)


def U_pri_int(list1,list2,nucleus,ncharge):
    # Return integrals for each pair of gaussians
    # and nucleus.
    p = list1[3]
    q = list2[3]
    a = p + q
    mu =  a
    R = C_N_dis(list1,list2,nucleus)
    T = mu*R*R
    
    if T > 15.0:
        return -1.0*ncharge* erf(sqrt(mu)*R)/R
    else:
        return -2.0* ncharge * sqrt(mu/pi)*F0(T)

   
def att_loop(input_list,operator,nucleus,ncharge):
    #Compute integral for each nucleus.
    M = []  
    for i in input_list:
        row = []
        for j in input_list:
            row += [operator(i,j,nucleus,ncharge)]
        M += [row]
    return M   


def nucl_loop(input_list,nulist,A_num_list):
    #Loop over nuclei, accumulate matrix for each nucleus.
    M = build_mat(input_list)
    #print(len(_nulist_))
    for i in range(len(nulist)):
        T = att_loop(input_list,U_pri_int,nulist[i],A_num_list[i])
        #print("DEBUG:")
        #printline(T)
        #print('=============================================')
        M = mat_add(M,T)
        #printline(M)
    #print("FINAL M:")
    #printline(M)
    
    return M
    

def att_mat(input_list,nulist,A_num_list):
    # mat_loop() to compute coefficients
    # sum over nuclei for the rest part of integrals
    # multiply coefficients and integrals correspondingly
    C = mat_loop(input_list,U_pri_coeff)
    M = nucl_loop(input_list,nulist,A_num_list)
    M = mat_multiply(M,C)
    return M

##############################################################################


########## ELECTRONIC REPULSTION INTEGRAL  ###################################    
    
def V_pri_combine(list1,list2):
    p = list1[3]
    q = list2[3]
    a = p + q
    mu = p * q / a
    C = center_coor(list1,list2)
    R2 = distance2(list1,list2)   
    # To return a and coefficient, attach them to C[3] and C[4].
    C[3] = a
    C[4] = exp(-mu*R2)*pow(pi/a,1.5)
    return C


def V_pri_coeff(p,q,r,s):
    C1 = V_pri_combine(p,r)
    C2 = V_pri_combine(q,s)
    return C1[4]*C2[4]


def V_pri_int(p,q,r,s):
    # list1 and list2 are combined new lists
    # in the format of [x,y,z,exp,cont]
    # integral goes this way:
    # Loop: <p(1)|<q(2)|V|r(1)>|q(2)>
    # but combine Psi_p(1)Psi_r(1)V Psi_q(2)Psi_s(2)
    C1 = V_pri_combine(p,r)
    C2 = V_pri_combine(q,s)
    
    a1 = C1[3]
    a2 = C2[3]
    a = a1 + a2
    mu = a1 * a2 / a
    R = sqrt(distance2(C1,C2))
    T = mu*R*R
    
    if T > 15.0:
        return erf(sqrt(mu)*R)/R
    else:
        return 2.0*sqrt(mu/pi)*F0(T)  

def rep_mat(_list_):
    C = rep_loop(_list_,V_pri_coeff)
    M = rep_loop(_list_,V_pri_int)
    M = mat_multiply(M,C)
    return M




###############################################################################



############################# MAIN FUNCTION ###################################

def pri_main(input_list,nulist,num_e,A_num_list):
    # INPUT: _list_: (x, y, z, exponential, contraction coefficient)
    # OUTPUT:
    #   Overlap,
    #   Kinetic,
    #   Attraction,
    #   Repulsion Matrices
    #
    # They are ALL IN PRIMITIVE BASIS!
    # 

    
    S = overlap_mat(input_list)
    T = kin_mat(input_list)
    N = att_mat(input_list,nulist,A_num_list)
    if num_e > 1:
        V = rep_mat(input_list)
    else:
        V = []

    return S,T,N,V


