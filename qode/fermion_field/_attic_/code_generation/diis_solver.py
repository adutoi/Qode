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
import numpy as np


def diis_solver(diis_subspace, t_amps, dt_amps, t_amps_list, diis_error_vecs):
	# diis_subspace = integer, others are NUMPY column matrices or list of NUMPY column matrices,
	# which can be shared by electronic and oscillator code.
	#
	if t_amps_list == None and diis_error_vecs == None:
		new_t_amps_list     = [ copy.deepcopy(t_amps) ]
		new_diis_error_vecs = [ copy.deepcopy(dt_amps) ]
		num_subspace        = 1
	else:
		if len(t_amps_list) != len(diis_error_vecs):
			print("Data Corruption: Data Not Aligned!")
			raise AssertionError
		if len(diis_error_vecs) < diis_subspace:
			new_t_amps_list     = copy.deepcopy(t_amps_list)     + [ copy.deepcopy(t_amps)  ]  # List + one element
			new_diis_error_vecs = copy.deepcopy(diis_error_vecs) + [ copy.deepcopy(dt_amps) ]
		else:
			new_t_amps_list     = copy.deepcopy(t_amps_list[1:])     + [ copy.deepcopy(t_amps) ] 
			new_diis_error_vecs = copy.deepcopy(diis_error_vecs[1:]) + [ copy.deepcopy(dt_amps) ]
		num_subspace = len(new_diis_error_vecs)

		
	# print("DIIS ERROR VECS =")
	# print(new_diis_error_vecs)

	# print("t_amps list=")
	# print(new_t_amps_list)

	 
	print("NUM SUBSPACE =", num_subspace)
	error_mat = -1.0 * np.matrix( np.ones(( num_subspace+1, num_subspace+1 )) )
	error_mat[-1,-1] = 0.0    # Last element set to zero
	for i in range(num_subspace):
		for j in range(num_subspace):
			error_mat[i,j] = new_diis_error_vecs[i].T * new_diis_error_vecs[j]

	# print("ERROR MATRIX")
	# print(error_mat)
	print(error_mat)

	error_mat_inv = np.linalg.inv(error_mat)
	# print("INVERTED ERROR MATRIX")
	# print(error_mat_inv)

	right_hand_vec = np.matrix( np.zeros((num_subspace+1,1)) )
	right_hand_vec[-1,0]  = -1.0

	# print("Right Hand Side Vec =")
	# print(right_hand_vec)


	weights = (error_mat_inv * right_hand_vec)[:-1,:]  # DO NOT INCLUDE THE LAST ELEMENT WHICH IS LAMBDA

	print("WEIGHTS")
	print(weights)
	
	# Try to normalize the weights but NOT SURE if it is necessary.
	# total_weight = 0.0
	# for i in range(num_subspace):
	# 	total_weight += weights[i,0]
	# norm_weights = copy.deepcopy(weights)
	# norm_weights = norm_weights / total_weight

	# print("normalized weights")
	# print(norm_weights)
	

	t_amps_short = np.matrix( np.zeros(( dt_amps.shape[0], 1)) )

	for i in range(num_subspace):
		t_amps_short += weights[i,0] * ( new_t_amps_list[i] + new_diis_error_vecs[i] ) 

	# print("SHORT t amplitudes =")
	# print(t_amps_short)

	t_amps_long  = np.concatenate( (np.zeros((1,1)), t_amps_short), axis=0)
	# print("LONG t amplitudes =")
	# print(t_amps_long)

	return t_amps_long, new_t_amps_list, new_diis_error_vecs











