#    (C) Copyright 2018 Anthony D. Dutoi
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
import copy
import time
import numpy as np



# DIIS solver Imports
from diis_solver import diis_solver
from ccsd_omega0 import cisd_vec_product, omega0_act_on_ref
import cisd_amplitude





#------------------------------ CCSD ENERGY FUNCTION ------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------


# def ccsd_energy(t_amp_obj , H, T, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
# 	#
# 	#
# 	print("Solving CCSD Energy...")
# 	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
# 		generate_elec_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
# 	#
# 	#
# 	occ_orbs  = alpha_orbs + beta_orbs
# 	virt_orbs = alpha_virt_orbs + beta_virt_orbs
# 	#
# 	dimension = generate_elec_config.get_num_reference() + generate_elec_config.get_num_singles(num_occ_orb,num_virt_orb) +\
# 	            generate_elec_config.get_num_doubles(num_occ_orb,num_virt_orb) 
# 	#
# 	#
# 	omega0_list = omega0_ops(H, T)
# 	flat_omega0_list = flatten_operator_list(omega0_list)
# 	del omega0_list
# 	# flat_omega0_list.print_info()
# 	#
# 	#
# 	hf_state    = state.state(t_amp_obj.get_configs())
# 	hf_state.coeffs[0] = 1.0
# 	#
# 	#
# 	return_state = state.state(t_amp_obj.get_configs())
# 	#
# 	# for op_list_obj in omega0_list:
# 	# 	# Top level operator list element
# 	# 	print("-----------------------")

# 	operators = flat_omega0_list.get_operators()
# 	coeffs    = flat_omega0_list.get_coeffs()

# 	input_list = []
# 	for i in range(len(coeffs)):
# 		input_list += [[ operators[i], coeffs[i], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat ]]
	
# 	# # print("NUM TERMS =", len(omega0_list))
# 	omega0_pool = Pool( len(coeffs) )
# 	list_return_states = omega0_pool.map(omega0_wrapper, input_list)
# 	omega0_pool.terminate()
# 	#
# 	#
# 	# Combine all states.
# 	for i in range(len(list_return_states)):
# 		return_state = state.add_state( return_state, list_return_states[i] )


# 	return state.state_dot(hf_state, return_state), return_state


def fock_get_diag_elemnt(fock_mat):
	size = fock_mat.shape[0]
	fock_orb_energy = []
	for i in range(size):
		fock_orb_energy += [ copy.deepcopy(fock_mat[i,i]) ]
	return fock_orb_energy


def compute_omega0_vector(ccsd_amp_obj, t_amp_obj, F_mat, V_mat, num_threads):
	computed_vec = omega0_act_on_ref(ccsd_amp_obj, t_amp_obj, F_mat, V_mat, num_threads)
	omega0_size  = ccsd_amp_obj.get_vec_dimension()
	omega0_vec   = np.zeros((omega0_size,1))

	num_alpha_elec  = ccsd_amp_obj.get_num_alpha_elec()
	num_beta_elec   = ccsd_amp_obj.get_num_beta_elec()
	num_spatial_orb = ccsd_amp_obj.get_num_spatial_orb()
	num_spin_orb    = num_spatial_orb * 2
	num_occ_orb     = num_alpha_elec + num_beta_elec
	num_vrt_orb     = num_spin_orb   -  num_occ_orb

	omega0_idx = 0
	mu_bra_state = copy.deepcopy( ccsd_amp_obj )
	mu_bra_state.clean_all_amplitude()
	mu_bra_state.update_ref_amplitude(1.0)
	omega0_vec[omega0_idx,0] = cisd_vec_product(mu_bra_state, computed_vec) 

	omega0_idx += 1
	for i in range(num_occ_orb):
		for a in range(num_vrt_orb):
			mu_bra_state = copy.deepcopy( ccsd_amp_obj )
			mu_bra_state.clean_all_amplitude()
			temp_single_mat = copy.deepcopy(mu_bra_state.get_single_amplitude())
			temp_single_mat[i,a] = 1.0
			mu_bra_state.update_single_amplitude( copy.deepcopy(temp_single_mat) )
			omega0_vec[omega0_idx,0] = cisd_vec_product(mu_bra_state, computed_vec)
			omega0_idx += 1

	for i in range(num_occ_orb):
		for j in range(i+1, num_occ_orb):
			for a in range(num_vrt_orb):
				for b in range(a+1, num_vrt_orb):
					mu_bra_state = copy.deepcopy( ccsd_amp_obj )
					mu_bra_state.clean_all_amplitude()
					temp_double_mat = copy.deepcopy(mu_bra_state.get_double_amplitude())
					# print(temp_double_mat.shape)
					# print(i,j,a,b,'==>',i*num_vrt_orb+a, j*num_vrt_orb+b)
					temp_double_mat[i*num_vrt_orb+a, j*num_vrt_orb+b] = 1.0
					# print(temp_double_mat)
					mu_bra_state.update_double_amplitude( copy.deepcopy(temp_double_mat) )
					omega0_vec[omega0_idx,0] = cisd_vec_product(mu_bra_state, computed_vec)
					omega0_idx += 1
	return omega0_vec




# def get_excitation_idx(num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
# 	excitation = []
# 	for i in range(num_occ_orb):
# 		orb_i = occ_orbs[i]
# 		for a in range(num_virt_orb):
# 			orb_a = virt_orbs[a]
# 			excitation += [[orb_a,orb_i]]
# 	for i in range(num_occ_orb):
# 		orb_i = occ_orbs[i]
# 		for j in range(i+1, num_occ_orb):
# 			orb_j = occ_orbs[j]
# 			for a in range(num_virt_orb):
# 				orb_a = virt_orbs[a]
# 				for b in range(a+1, num_virt_orb):
# 					orb_b = virt_orbs[b]
# 					excitation += [[orb_a, orb_i, orb_b, orb_j]]
# 	return excitation




# THIS IS BAD B/C IT IMPLICITLY USES INTERNAL DATA STRUCTURE, LET'S BUILD THIS ON THE FLIGHT.
# def compute_orb_energy_diff(fock_orb_energy, excitation_idx):
# 	energy_diff_list = []
# 	for item in excitation_idx:
# 		if len(item) == 2:
# 			energy_diff_list += [ fock_orb_energy[item[0]] - fock_orb_energy[item[1]] ]
# 		elif len(item) == 4:
# 			energy_diff_list += [ fock_orb_energy[item[0]] - fock_orb_energy[item[1]] + fock_orb_energy[item[2]] - fock_orb_energy[item[3]] ]
# 	return energy_diff_list

def compute_orb_energy_diff(fock_orb_energy, occ_orb_idx, vrt_orb_idx):
	# fock_orb_energy is a 1D array.
	num_occ_orb = len(occ_orb_idx)
	num_vrt_orb = len(vrt_orb_idx)
	inverted_fock_diff = []
	for i in range(num_occ_orb):
		for a in range(num_vrt_orb):
			inverted_fock_diff += [ fock_orb_energy[vrt_orb_idx[a]] - fock_orb_energy[occ_orb_idx[i]] ]
	for i in range(num_occ_orb):
		for j in range(i+1, num_occ_orb):
			for a in range(num_vrt_orb):
				for b in range(a+1, num_vrt_orb):
					inverted_fock_diff += [ fock_orb_energy[vrt_orb_idx[a]] - fock_orb_energy[occ_orb_idx[i]] + fock_orb_energy[vrt_orb_idx[b]] - fock_orb_energy[occ_orb_idx[j]] ]
	dimension = len(inverted_fock_diff)
	return np.array(inverted_fock_diff).reshape((dimension,1))

def compute_delta_t_vec(omega0_vector, orb_energy_diff_list):
	size    = omega0_vector.shape[0] - 1  # Long vec
	dt_amps = np.matrix(np.zeros((size,1)))  # short vec
	for i in range(size):
		dt_amps[i,0] = -1.0/orb_energy_diff_list[i] * omega0_vector[i+1,0]  # Need plus one
	return dt_amps




#-------------------------------------------- MAIN FUNCTION -------------------------------------------------
#------------------------------------------------------------------------------------------------------------



def ccsd_main(F_mat, V_mat, fock_orb_energy, nuclear_repul_energy, num_alpha_elec, num_beta_elec, num_spatial_orb, num_max_subspace, num_itr_start_diis):
	if num_max_subspace < 1:
		print("DIIS Subspace Must Be Greater or Equal to One")
		raise ValueError
	if num_itr_start_diis < 1:
		print("DIIS Iteration Must At Least Start at Cycle One")
		raise ValueError
	#
	#
	ccsd_diis_subspace = num_max_subspace
	print('==================================================================================')
	print('CCSD Solver Starts ...')
	print("DIIS Max Subspcae Size =", ccsd_diis_subspace)
	print("DIIS Will Be Invoked at Cycle ", num_itr_start_diis)
	#
	#
	#
	#
	num_spin_orb    = num_spatial_orb * 2
	num_occ_orb     = num_alpha_elec+num_beta_elec
	num_vrt_orb     = 2 * num_spatial_orb - num_alpha_elec - num_beta_elec
	#
	occ_orb_idx = [ i for i in range(num_alpha_elec)] + [ i+num_spatial_orb for i in range(num_beta_elec)]
	vrt_orb_idx = [ i+num_alpha_elec for i in range(num_spatial_orb - num_alpha_elec) ] + [ i+num_spatial_orb+num_beta_elec for i in range(num_spatial_orb - num_beta_elec) ]
	#
	#
	# 
	#  CCSD CYCLES START HERE: 
	cyc_count = 1

	print('==================================================================================')
	print("CCSD CYCLE", cyc_count)
	
	ccsd_amp_obj = cisd_amplitude.cisd_amplitude(num_alpha_elec, num_beta_elec, num_spatial_orb)
	

	t_amp_obj = cisd_amplitude.cisd_amplitude(num_alpha_elec, num_beta_elec, num_spatial_orb)
	t_amp_obj.update_ref_amplitude(0.0)

	t_ccsd_start = time.time()

	print("Start building Vector Omega(0)\n")
	# t0_start = time.time()
	omega0_vector = compute_omega0_vector(ccsd_amp_obj, t_amp_obj, F_mat, V_mat, num_threads=20)
	# t0_end   = time.time()
	# np.savetxt("OMEGA0_time.txt", np.array([[t0_end - t0_start]]) )
	print("Omega(0) Vector Built.")
	print(omega0_vector)
	#
	#
	# print("Solving Vector dt")
	#
	orb_energy_diff_list = compute_orb_energy_diff(fock_orb_energy, occ_orb_idx, vrt_orb_idx)
	dt_amps = compute_delta_t_vec(omega0_vector, orb_energy_diff_list)
	#
	#
	#
	#
	if num_itr_start_diis == 1:
		print("Quasi-Newton PT using DIIS ...")
		t_amps = np.matrix( np.zeros((ccsd_amp_obj.get_vec_dimension()-1, 1)) )
		t_amps_list = None
		diis_error_vecs = None
		t_amps_long, t_amps_list, diis_error_vecs = diis_solver(ccsd_diis_subspace, t_amps, dt_amps, t_amps_list, diis_error_vecs)
		t_amp_obj.update_amp_from_long_vec( t_amps_long )
	else:
		print("Quasi-Newton PT Solver")
		dt_amps_long = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
		t_amp_obj.update_amp_from_long_vec( dt_amps_long )

	# t_amp_obj.print_info()
	print("New t Amplitudes Done.")
	omega0_vector = compute_omega0_vector(ccsd_amp_obj, t_amp_obj, F_mat, V_mat, num_threads=20)
	E_CCSD = nuclear_repul_energy + omega0_vector[0,0]
	print('--------------------------------------------------------------')
	print("E_CCSD =", E_CCSD)
	print('--------------------------------------------------------------')

	# sys.exit(0)

	cyc_count += 1

	
	converged = False
	while not converged:
		print('==================================================================================')
		print("CCSD CYCLE", cyc_count)
		# print("Start Solving Omega_0")
		# print("Omega_0 Done.")
		#
		#
		dt_amps = compute_delta_t_vec(omega0_vector, orb_energy_diff_list)
		#
		#
		if cyc_count < num_itr_start_diis:
			print("Quasi-Newton PT Solver")
			dt_amps_long = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
			t_amp_obj.update_amp_from_long_vec( t_amp_obj.amp_obj_to_long_vec() + dt_amps_long )
		else:
			print("Quasi-Newton PT using DIIS ...")
			#
			#
			if cyc_count == num_itr_start_diis: 
				# Check if t_amps_list and diis_error_vecs they exist
				t_amps_list     = None
				diis_error_vecs = None
				t_amps_long, t_amps_list, diis_error_vecs = diis_solver(ccsd_diis_subspace, t_amps_long[1:,:], dt_amps, t_amps_list, diis_error_vecs)
				t_amp_obj.update_amp_from_long_vec( t_amps_long )
			
			else:
				t_amps_long, t_amps_list, diis_error_vecs = diis_solver(ccsd_diis_subspace, t_amps_long[1:,:], dt_amps, t_amps_list, diis_error_vecs)
				t_amp_obj.update_amp_from_long_vec( t_amps_long )
				
		#
		#
		print("New t Amplitudes Done.")
		omega0_vector = compute_omega0_vector(ccsd_amp_obj, t_amp_obj, F_mat, V_mat, num_threads=20)
		print(omega0_vector)
		E_CCSD_new = nuclear_repul_energy + omega0_vector[0,0]
		print('--------------------------------------------------------------')
		print("E_CCSD =", E_CCSD_new )
		print('--------------------------------------------------------------')
		#
		#
		E_diff = abs(E_CCSD - E_CCSD_new) /abs(E_CCSD) 
		print("E_diff =", E_diff)
		if E_diff < 1e-10:
			print("CCSD Converged!")
			converged = True
		else:
			E_CCSD = E_CCSD_new
			cyc_count += 1
	
	print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")	
	print("CCSD converged at Cycle", cyc_count)
	print("CCSD Energy =", E_CCSD_new, 'Hartree')
	print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")	



	t_ccsd_end = time.time()

	np.savetxt("CCSD_TOTAL_TIME.txt", np.array([[ t_ccsd_end - t_ccsd_start ]]) )

	# np.save("ccsd_vec.npy", t_amp_obj.get_coeffs())
	f = open("ccsd_energy.txt",'w')
	f.write("CCSD converged at Cycle %d\n" %(cyc_count) )
	f.write("CCSD Energy = %.16f Hartree" %(E_CCSD_new)  )
	f.close()
	return E_CCSD_new	







if __name__ == "__main__":
	#
	# Helium Atom:
	# SCF energy in the final basis set = -2.8551604262  Hartree (Qchem4.4)
	# MP2 energy                        = -2.86814954    Hartree (Qchem4.4)
	# CCSD total energy                 = -2.87178917    Hartree (Qchem4.4)
	# QCHEM 4.4 DIIS_SUBSPACE_SIZE = 15 by default
	#
	# main(num_max_subspace, num_itr_start_diis)
	# iteration number starts at one!
	#
	# inputs = [ '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', '1.9']
	# for first_digit in range(2,8):
	# 	for second_digit in [0, 25, 50, 75]:
	# 		inputs += [str(first_digit) + '.' + str(second_digit)]
	# inputs +=[ '8.0']

	# print(inputs)
	# inputs = ['1.0']

	# for item in inputs:
	# 	t1 = time.time()
	# 	input_file = item + '.in'
	# 	E_0, KE, NAE, NRE, ERE, F_mat, T_mat, N_mat, V_mat, vHF, num_alpha_elec, num_beta_elec, num_spin_orb = scf_routine.main(input_file, 1e-10)
	# 	E_hf_orb = fock_get_diag_elemnt(F_mat)
	# 	h_mat = T_mat + N_mat
	# 	E_CCSD = main(h_mat, V_mat, E_hf_orb, NRE, num_alpha_elec, num_beta_elec, num_spin_orb, num_max_subspace=5, num_itr_start_diis=1)
	# 	f = open(item + '.out', 'w')
	# 	f.write('E_CCSD = ' + str(E_CCSD))
	# 	f.close()
	# 	t2 = time.time()
	# 	print("CCSD TOTAL TIME =", t2 - t1, "seconds")
	# input_file = '/home/yhliu/Qode-clone/Applications/ccsd/1.4.in'
	# E_0, KE, NAE, NRE, ERE, F_mat, T_mat, N_mat, V_mat, vHF, num_alpha_elec, num_beta_elec, num_spin_orb = scf_routine.main(input_file, 1e-10)
	# # E_hf_orb = fock_get_diag_elemnt(F_mat)
	# h_mat = T_mat + N_mat
	# np.save('1.4.h_mat.npy', h_mat)
	# np.save('1.4.V_mat.npy', V_mat)
	nuclear_repul_energy = 0.0
	F_mat  = np.load('F_mat.npy')
	V_mat  = np.load('V_mat.npy')
	fock_orb_energy = fock_get_diag_elemnt(F_mat)
	num_alpha_elec  = 1
	num_beta_elec   = 1
	num_spatial_orb = 4
	# num_alpha_elec  = 2
	# num_beta_elec   = 2



	E_CCSD =ccsd_main(F_mat, V_mat, fock_orb_energy, nuclear_repul_energy, num_alpha_elec, num_beta_elec, num_spatial_orb, num_max_subspace=5, num_itr_start_diis=1)

	# f = open('/home/yhliu/Qode-clone/Applications/ccsd/' + item + '.out', 'w')
	# f.write('E_CCSD = ' + str(E_CCSD))
	# f.close()















