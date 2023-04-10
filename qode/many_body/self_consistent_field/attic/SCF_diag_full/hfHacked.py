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
from qode.SCF_diag_full.matrix_operation import printline,printnpline,writenpline
from qode.SCF_diag_full import molecule
from qode.SCF_diag_full import nuclear_list
from qode.SCF_diag_full import nucl_repul
from qode.SCF_diag_full import contracted
from qode.SCF_diag_full import changebasis
from qode.SCF_diag_full import Fock_test
from qode.SCF_diag_full import sum_repulsion
from qode.SCF_diag_full import scfHacked as scf
from qode.SCF_diag_full import mp2


np.set_printoptions(precision=3,linewidth=270,threshold=np.nan)


def main(InputFile,threshold=1.0e-8):

    print("INPUT FILE NAME:",InputFile)
    print("Computation Thresthold:",threshold)
    
    input_list, ngaussians, orbitype, atoms, num_e, multi, A_num_list,correlation = molecule.mol_main(InputFile)
    nulist =  nuclear_list.numain(InputFile)
    
    print('nuclei list and charge:')
    print(nulist)
    print(A_num_list)
    
    print("Num Gaussian =", ngaussians)

    NuRepulE = nucl_repul.nurepul_main(nulist,A_num_list)
    S_c, T_c, N_c, V_c = contracted.cont_main(input_list, ngaussians, nulist,num_e,A_num_list)

    # Before this line, all matrices in Python Arrays.
    # After this line, they are all in Numpy Arrays due to diagonalization.

    E, U_on, T_on, N_on, V_on = changebasis.chbasis_main(S_c, T_c, N_c, V_c, num_e)

    #Fock_test.fock_main(U_on,T_on,N_on,V_on,input_list,ngaussians)

    # Before this line, Alpha Block == Beta Block
    # After this line, they may not be equal ( multiplets ).
    return scf.scf_main(T_on, N_on, V_on, num_e, multi, threshold, NuRepulE)				# sorry about this nasty hack ... just want access to these for some testing
    EHF,KE,NAE,NRE,ERE,E,U,T,N,V,alpha_e,beta_e = scf.scf_main(T_on, N_on, V_on, num_e, multi, threshold, NuRepulE)
    
    # np.savetxt("E.txt",E)
    # np.savetxt("V.txt",V)

    if correlation == 'mp2':
        EMP2 = mp2.mp2_main(EHF,E,V,alpha_e,beta_e)
    
    #input("Press Any Key to Continue...")
    # np.set_printoptions(precision=3,threshold=np.nan,linewidth=200)
    # print("Eigen Values =",E)
    # print("U matrix =")
    # print(U)
    # print("V matrix =")
    # print(V)
    # print("eigenvals =",E)

    # Get Fock Matrix elec-elec part vHF

    vHF = sum_repulsion.sum_repul(E,V,alpha_e,beta_e)

    # Make Fock matrix as a diagonal matrix using the numbers in E
    F = np.matrix( np.zeros((E.shape[0], E.shape[0])) )
    for i in range(E.shape[0]):
        F[i,i] = E[i]
        

    num_spin_orb = 2 * len(ngaussians)
    return EHF,KE,NAE,NRE,ERE, F, T, N, V, vHF, alpha_e, beta_e, num_spin_orb



if __name__=='__main__':
    input_file = 'he.in'
    main(input_file,threshold)













    
