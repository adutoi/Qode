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
from ccsd_code_normal_ordered import *  # All the h_ref, h_single and h_double ...
import numpy as np
import copy
from cisd_amplitude import add_and_update_cisd_amp_obj
from multiprocessing import Pool



def h_wrapper(input_args):
	return input_args[0](*input_args[1:])


def omega0_act_on_ref(ccsd_amp_obj, t_amp_obj, F_mat, V_mat, num_threads):
	result_amp_obj = copy.deepcopy(ccsd_amp_obj)
	result_amp_obj.clean_all_amplitude()

	list_of_h = [h0,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,\
	             h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,\
	             h21,h22,h23,h24,h25,h26,h27,h28,h29,h30,\
	             h31,h32,h33,h34,h35,h36,h37,h38,h39,h40,\
	             h41,h42,h43,h44,h45,h46,h47,h48,h49,h50,\
	             h51,h52,h53,h54,h55,h56,h57,h58,h59,h60,\
	             h61,h62,h63,h64,h65,h66,h67,h68,h69,h70,\
	             h71,h72,h73,h74,h75,h76,h77,h78,h79,h80,h81]

	list_of_inputs = [[list_of_h[i], ccsd_amp_obj, t_amp_obj, F_mat, V_mat] for i in range(len(list_of_h))]

	myPool = Pool(num_threads)
	return_list = myPool.map(h_wrapper, list_of_inputs)
	myPool.terminate()
	
	for item in return_list:
		result_amp_obj = add_and_update_cisd_amp_obj( result_amp_obj, item )



	return result_amp_obj





def cisd_vec_product(cisd_amp_obj1, cisd_amp_obj2):
	product  = cisd_amp_obj1.get_ref_amplitude() * cisd_amp_obj2.get_ref_amplitude()
	ld1, ld2 = cisd_amp_obj1.get_single_amplitude().shape
	single_amp_mat1 = cisd_amp_obj1.get_single_amplitude()
	single_amp_mat2 = cisd_amp_obj2.get_single_amplitude()
	for i in range(ld1):
		for j in range(ld2):
			product += single_amp_mat1[i,j] * single_amp_mat2[i,j]
	double_amp_mat1 = cisd_amp_obj1.get_double_amplitude()
	double_amp_mat2 = cisd_amp_obj2.get_double_amplitude()
	ld1, ld2 = double_amp_mat1.shape
	
	for i in range(ld1):
		for j in range(ld2):
			product += double_amp_mat1[i,j] * double_amp_mat2[i,j]
	return product


if __name__ == "__main__":
	import cisd_amplitude
	np.set_printoptions(precision=3,linewidth=270,threshold=np.nan)
	num_ae = 1
	num_be = 1
	# num_ae = 2
	# num_be = 2
	num_spat = 4
	num_occ = num_ae + num_be
	num_vrt = 2 * num_spat - num_occ
	ccsd_vec = cisd_amplitude.cisd_amplitude(num_ae, num_be, num_spat)
	ccsd_vec.update_ref_amplitude(1.0)
	ccsd_vec.print_info()
	# t is empty
	t_amp_obj = cisd_amplitude.cisd_amplitude(num_ae, num_be, num_spat)
	F_mat = np.load('F_mat.npy')
	V_mat = np.load('V_mat.npy')
	# E_0   = -2.85516042615

	# print(F_mat)
	# print(V_mat)

	first_state = omega0_act_on_ref(ccsd_vec, t_amp_obj, F_mat, V_mat, num_threads=4)
	print(cisd_vec_product(ccsd_vec, first_state))




# def cisd_main(cisd_amp_obj, F_mat, V_mat, num_occ, num_vrt):
# 	list_of_CISD_states  = []
# 	list_computed_states = []
# 	# H |R>
# 	ref_amp_obj = copy.deepcopy(cisd_amp_obj)
# 	list_of_CISD_states += [copy.deepcopy(ref_amp_obj)]
# 	# list_computed_states += [copy.deepcopy( H(ref_amp_obj, E_0, F_mat, V_mat) )]

# 	# All Singles
# 	single_amp_obj = copy.deepcopy(cisd_amp_obj)
# 	single_amp_obj.clean_all_amplitude()
# 	for i in range(num_occ):
# 		for j in range(num_vrt):
# 			single_amp_mat = np.zeros( single_amp_obj.get_single_amplitude().shape )
# 			single_amp_mat[i,j] = 1.0
# 			single_amp_obj.update_single_amplitude( single_amp_mat )
# 			# single_amp_obj.print_info()
# 			list_of_CISD_states += [copy.deepcopy( single_amp_obj )]

# 	# All doubles
# 	double_amp_obj = copy.deepcopy(cisd_amp_obj)
# 	double_amp_obj.clean_all_amplitude()
# 	for m in range(num_occ):
# 		for e in range(num_vrt):
# 			for n in range(m+1, num_occ):
# 				for f in range(e+1, num_vrt):
# 					double_amp_mat = np.zeros( double_amp_obj.get_double_amplitude().shape )
# 					double_amp_mat[m*num_vrt+e, n*num_vrt+f] = 1.0
# 					double_amp_obj.update_double_amplitude( double_amp_mat )
# 					# print( "C[%d,%d,%d,%d] => C[%d,%d]" %(m,n,e,f, m*num_vrt+e, n*num_vrt+f ))
# 					# double_amp_obj.print_info()
# 					list_of_CISD_states += [copy.deepcopy( double_amp_obj )]

# 	# print(list_of_CISD_states)
# 	# print(len(list_of_CISD_states))
# 	for cisd_state in list_of_CISD_states:
# 		list_computed_states += [copy.deepcopy( H(cisd_state, F_mat, V_mat) )]

# 	# print(list_computed_states)
# 	ld = len(list_computed_states)
# 	h_mat = np.matrix( np.zeros((ld,ld)) )
# 	for i in range(ld):
# 		for j in range(ld):
# 			h_mat[i,j] = cisd_vec_product(list_of_CISD_states[i], list_computed_states[j])
# 	print(h_mat)

# 	E,V = np.linalg.eigh(h_mat)

# 	cisd_energy = E[0]
# 	print("CISD ENERGY =", cisd_energy)

# 	return cisd_energy













