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
from qode.SCF_Code.matrix_operation import printnpline
from qode.SCF_Code import elec_prescf
from qode.SCF_Code import elec_fock_builder
from qode.SCF_Code import generic_scf
from qode.SCF_Code import elec_HF_energy
from qode.SCF_Code import elec_mat2fock

def copy_alpha_vec_to_beta(U):
	dimension = U.shape[0]
	num_spatial_orb = dimension // 2
	for i in range(num_spatial_orb):
		for j in range(num_spatial_orb):
			U[i+num_spatial_orb,j+num_spatial_orb] = U[i,j]
	



def scf_main( file, thresh = 1.0e-10 ):
	np.set_printoptions(threshold=np.nan,linewidth=260,precision=3)

	guess, U0, T_sp, N_sp, V_sp, alpha_e, beta_e, NRE = elec_prescf.main( file )
	print("Alpha e-:",alpha_e,"Beta e-:", beta_e)
	print("Guess Vectors:")
	printnpline(guess)

	fock = elec_fock_builder.fock_builder(T_sp, N_sp, V_sp, alpha_e, beta_e, NRE)

	orb = generic_scf.SCF(fock,guess,thresh)

	U = elec_fock_builder.build_basis(orb, alpha_e, beta_e)
	print("UHF U Matrix:")
	print(U)


	copy_alpha_vec_to_beta(U)
	print("U")
	print(U)


	F_scf, T_scf, N_scf, V_scf = elec_mat2fock.build_Fock(U, T_sp, N_sp, V_sp, alpha_e, beta_e)

	print("FOCK MATRIX:")
	print(F_scf)
	# print(V_scf)
	

	EHF,KE,NAE,ERE = elec_HF_energy.HF_ENERGY(U, T_sp, N_sp, V_sp, alpha_e, beta_e, NRE)
	print("HF ENERGY =",EHF)
	# print("KE =",KE)
	# print("NAE =",NAE)
	# print("ERE =",ERE)

	# print("SPATIAL ORB DIM =", T_scf.shape[0] // 2)
	num_spatial_orb = T_sp.shape[0] // 2

	return EHF, F_scf, V_scf, alpha_e, beta_e, num_spatial_orb
