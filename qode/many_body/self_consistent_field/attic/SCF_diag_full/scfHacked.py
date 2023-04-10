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
from qode.SCF_diag_full import repul_transform
from qode.SCF_diag_full import sum_repulsion
from qode.SCF_diag_full import eigensort
from qode.SCF_diag_full import build_spin
from math import sqrt
from time import ctime
from qode.SCF_diag_full.matrix_operation import printline,printnpline,writenpline,build_size_mat
from qode.SCF_diag_full.matrix_operation import build_same_np_mat,empty_same_np_mat

np.set_printoptions(precision=6,linewidth=270,threshold=np.nan)


################## SHARED ROUTINE  ##################################
def ab_spin(num_e,multi):
    # Return numbers of alpha and beta electrons.
        
    if ((num_e > 1) and (multi > 0)): 
        single_e = multi - 1
        if (num_e - single_e) < 0 or (num_e - single_e)%2 != 0:
            print("Invalid Charge/Electron Configuration")
            sys.exit(1)
        beta_e = (num_e - single_e) // 2
        alpha_e = beta_e + single_e
        
    elif (num_e == 1) and (multi == 2):
        alpha_e = 1
        beta_e = 0

    else:
        print("Invalid Charge/Electron Configuration")
        sys.exit(1)
        
    return alpha_e,beta_e

######################################################################



###################  HF ENERGY CODE  ##################################

def HF_energy(E,T,N,V, alpha_e,beta_e,NRE):
    "EHF = Sum_i [<i1|T|i1> + <i1|N|i1> ] + 1/2 * Sum_i,j <i1|<j0|V(1-P10)|i1>|j0>"
    #Initiate looping indicies
    n = len(T)//2
    #print(n)
    Ia,Ib = eigensort.sort_main(E)
    #print(Ia,Ib)
    # Make an alpha and a beta electrons summing-index array

    
    "Kinetic Energy, Nuclear Attraction E., Repulsion E."
    KE = 0.0
    NAE = 0.0
    ERE = 0.0
    # Kinetic and Attraction Summing Loops ========
    
    for k1 in range(alpha_e):
        i = Ia[k1]
        #print(i)
        KE +=  T[i,i]
        NAE += N[i,i]
    
    if beta_e > 0:
        for k2 in range(beta_e):
            j = Ib[k2]
            #print(j)
            KE +=  T[j,j]
            NAE += N[j,j]


    # Repulsion Summing Loops ======================
    if (alpha_e + beta_e > 1):
        for k1 in range(alpha_e):
            i = Ia[k1]
            for k2 in range(alpha_e):
                j = Ia[k2]
                ERE += ( V[i*2*n+j,i*2*n+j]-V[i*2*n+j,j*2*n+i] )

        if (beta_e >0):
            for k1 in range(alpha_e):
                i = Ia[k1]
                for k2 in range(beta_e):
                    j = Ib[k2]
                    ERE += ( V[i*2*n+j,i*2*n+j]-V[i*2*n+j,j*2*n+i] )
    
            for k2 in range(beta_e):
                j = Ib[k2]
                for k1 in range(alpha_e):
                    i = Ia[k1]
                    ERE += ( V[j*2*n+i,j*2*n+i]-V[j*2*n+i,i*2*n+j] )
                
            for k1 in range(beta_e):
                i = Ib[k1]
                for k2 in range(beta_e):
                    j = Ib[k2]
                    ERE += ( V[i*2*n+j,i*2*n+j]-V[i*2*n+j,j*2*n+i] )
                    

        
        # ONLY LOOP OVER ALPHA
    ###############################################################
    # HAVE TO ADDRESS THE PROBLEM OF ALPHA != BETA (Multiplets)
    # Keep Both Alpha Repulsion and Beta Repulsion added with N, T
    #################################################################
    ERE = ERE/2.0
    
    EHF = KE + NAE + NRE + ERE 
    return EHF,KE,NAE,ERE
######################################################################




############# CONVERGENCE CHECK  #####################################

def chk_converge(EHF1,EHF2,threshold):
    
    ###### Default assmues SCF is Done. ########
    loop = False
    if abs(EHF1 - EHF2) >= threshold:
        loop = True
    return loop


#######################################################################

########################################################################
def sorted_eigenval_rebuild_spin(E):
    # print("EIGEN SHAPE = ",E.shape[0])
    num_spatial_orb = E.shape[0] // 2
    E_alpha = np.array( [ E[i*2] for i in range(num_spatial_orb) ] )
    E_rebuilt = np.concatenate((E_alpha,E_alpha), axis=1)
    # print("NEW EIGEN VALUES =", E_rebuilt)
    return E_rebuilt


def sorted_vector_rebuild_spin(U):
    # Since numpy.linalg.eigh is used to diagonalize all Fock matrices,
    # Vector spins are mixed but sorted by energy. Therefore, this function
    # takes vectors and rebuild them in an alpha-beta pattern. Inside alpha 
    # or beta blocks, the vectors are sorted by energy.
    # 
    # print(U)
    num_spatial_orb = U.shape[0] // 2
    #
    dummpy_column = np.zeros((2*num_spatial_orb,1))  
    #
    for i in range(num_spatial_orb):
        if U[0,2*i] != 0.0:
            dummpy_column = np.concatenate((dummpy_column, U[ :,[2*i] ]), axis=1)
        else:
            dummpy_column = np.concatenate((dummpy_column, U[ :,[2*i+1] ]), axis=1)
    alpha_block = dummpy_column[:,1:]
    # print(alpha_block)
    # print("BETA TOP")
    # print(alpha_block[num_spatial_orb:,:])
    # print("BETA BOTTOM")
    # print(alpha_block[:num_spatial_orb:,:])
    beta_block  = np.concatenate((alpha_block[num_spatial_orb:,:], alpha_block[:num_spatial_orb:,:]), axis=0)
    # print("Beta Block")
    # print(beta_block)
    U_rebuilt = np.concatenate((alpha_block, beta_block), axis=1)
    # print("rebuilding DONE!!")
    # print(U_rebuilt)
    return U_rebuilt


############## SCF CYCLE STEPS ########################################
def scf_cycle(E1,U1,T0,N0,V0,alpha_e,beta_e):
    T1 = U1.H * T0 * U1
    N1 = U1.H * N0 * U1
    if (alpha_e + beta_e) > 1:
        V1 = repul_transform.transform_V(V0,U1)
        vHF = sum_repulsion.sum_repul(E1,V1,alpha_e,beta_e)
        #printnpline(vHF)
        F1 = T1 + N1 + vHF 
   
    else:
        F1 = T1 + N1
        V1 = V0


    E2, U2 = np.linalg.eigh(F1)
    E2 = sorted_eigenval_rebuild_spin(E2)
    U2 = sorted_vector_rebuild_spin(U2)
    
    return E2,U2,F1,T1,N1,V1




########################################################################
def printE(EHF,KE,NAE,NRE,ERE):
    
    print('Kinetic   E =',KE)
    print('Nucl Attr.E =',NAE)
    print('Nucl RepulE =',NRE)
    print('Elec.RepulE =',ERE)
    print("EHF =",EHF,'\n')



####################### Main Function #################################

    
def scf_main(T_on, N_on, V_on, num_e, multi, threshold,NuRepulE):

    NRE = NuRepulE
    T_sp,N_sp,V_sp = build_spin.spin_main(T_on, N_on, V_on)
    
    print("Sorry!  SCF code hacked for testing ... short-circuited to just return input matrices")
    return T_sp, N_sp, V_sp

    alpha_e,beta_e = ab_spin(num_e,multi)
    print('# of alpha e-:',alpha_e)
    print('# of beta e-:',beta_e)
    # Initial Guess: F(-1) = T_on + N_on
    E0,U0 = np.linalg.eigh(T_sp + N_sp)
    E0 = sorted_eigenval_rebuild_spin(E0)
    U0 = sorted_vector_rebuild_spin(U0)

    print("Fock Guess E =", E0)
    print("Guess U")
    print(U0)
    
    count = 1
    print("\nSCF starts at",ctime())
    print("\nCycle 1:")
    E1, U1, F0, T0, N0, V0 = scf_cycle(E0,U0,T_sp,N_sp,V_sp,alpha_e,beta_e)
    EHF0,KE0,NAE0,ERE0 = HF_energy(E0, T0, N0, V0, alpha_e,beta_e,NRE)
    
    printE(EHF0,KE0,NAE0,NRE,ERE0)
    
    
    scf_loop = True
    while scf_loop:
        count += 1
        print("Cycle",count,":")
        E2,U2,F1,T1,N1,V1 = scf_cycle(E1,U1,T0,N0,V0,alpha_e,beta_e)
        EHF1,KE1,NAE1,ERE1 = HF_energy(E1, T1, N1, V1, alpha_e,beta_e,NRE)
        printE(EHF1,KE1,NAE1,NRE,ERE1)
        scf_loop = chk_converge(EHF0, EHF1, threshold)
        EHF0 = EHF1
        KE0 = KE1
        NAE0 = NAE1
        ERE0 = ERE1
        E1 = E2
        U1 = U2
        T0 = T1
        N0 = N1
        V0 = V1
        #U = U * Ui
        #print("SCF in progress")
        #printnpline(F0)
        #print("EHF =",EHF1)
    
    #printnpline(U)
    print("Convergence criteria met at Cycle",count)
    #printnpline(F1)
    print('=====================================================')
    print("Final Energies:")
    printE(EHF0,KE0,NAE0,NRE,ERE0)
    print("Total HF ENERGY =",EHF0)
    print(' ')
    #E,U = get_unitary(F1)
    #print("E:\n",E)
    #printnpline(U)
    
    return EHF0,KE0,NAE0,NRE,ERE0,E1,U1,T0,N0,V0,alpha_e,beta_e
