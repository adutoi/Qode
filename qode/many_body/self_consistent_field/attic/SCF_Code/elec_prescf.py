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
from qode.SCF_Code.matrix_operation import printline,printnpline
from qode.SCF_Code import molecule
from qode.SCF_Code import nuclear_list
from qode.SCF_Code import nucl_repul
from qode.SCF_Code import contracted
from qode.SCF_Code import changebasis
from qode.SCF_Code import abspin
from qode.SCF_Code import build_spin
from qode.SCF_Code import elec_scf_guess






def main(InputFile):
    print("Pre-SCF process for electrons")
    print("INPUT FILE NAME:",InputFile)
    
    input_list, ngaussians, orbitype, atoms, num_e, multi, A_num_list, correlation = molecule.mol_main(InputFile)
    nulist =  nuclear_list.numain(InputFile)
    
    
    print('nuclei list and charge:')
    print(nulist)
    print(A_num_list) #Atomic Number List
    
    NuRepulE = nucl_repul.nurepul_main(nulist,A_num_list)
    S_c, T_c, N_c, V_c = contracted.cont_main(input_list, ngaussians, nulist,num_e,A_num_list)

    # Before this line, all matrices in Python Arrays.
    # After this line, they are all in Numpy Arrays due to diagonalization.

    E, U_on, T_on, N_on, V_on = changebasis.chbasis_main(S_c, T_c, N_c, V_c, num_e)
    
    alpha_e, beta_e = abspin.ab_spin(num_e,multi)
    
    T_sp,N_sp,V_sp = build_spin.spin_main(T_on, N_on, V_on)
    
    print("T spin matrix dimension =", len(T_sp) )
    print("N spin matrix dimension =", len(N_sp) )
    print("V spin matrix dimension =", len(V_sp) )
    
    guess, E0, U0 = elec_scf_guess.guess_main(T_sp, N_sp, alpha_e, beta_e)

    return guess, U0, T_sp, N_sp, V_sp, alpha_e, beta_e, NuRepulE


if __name__=='__main__':
    input_file='h2.in'
    main(input_file)
