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
#
# This code only works for TWO Harmonic Oscillator Molecules (dirty code)
# The code generates all the FCI configurations.
#
import copy
import numpy as np
from . import config
from . import state

from qode.coupled_oscillators import oscillator_system as osys
from qode.coupled_oscillators import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution



def generate_empty_config(list_of_num_state):
	vec_size = 0
	for size in list_of_num_state:
		vec_size += size
	return [False for i in range(vec_size)]



def generate_ref(list_of_num_state):
	config = generate_empty_config(list_of_num_state)
	idx = 0
	for num_state in list_of_num_state:
		# idx = num_state
		config[idx] = True
		idx += num_state

	return config



def generate_singles(list_of_num_state):
	num_mol = len(list_of_num_state)
	configs = []
	idx_shift = 0
	for i in range(num_mol):
		# Loop over "molecules"
		for j in range(list_of_num_state[i] - 1):
			# Loop over virtual states.
			this_config = generate_ref(list_of_num_state)
			this_config[ idx_shift     ] = False
			this_config[ idx_shift+j+1 ] = True
			configs += [ this_config ]
		idx_shift += list_of_num_state[i]
	return configs



def generate_doubles(list_of_num_state):
	num_mol = len(list_of_num_state)
	configs = []
	first_idx_shift  = 0
	for i in range(num_mol):               # Mol 1
		# Compute the second index shift
		#
		second_idx_shift = 0
		for num in list_of_num_state[:i+1]:
			second_idx_shift += num
		#
		#
		for j in range(i+1, num_mol):      # Mol 2
			for k in range(list_of_num_state[i] - 1):      # Virtual States in Mol 1
				for l in range(list_of_num_state[j] - 1):  # Virtual States in Mol 2
					this_config = generate_ref(list_of_num_state)
					this_config[ first_idx_shift    ] = False
					this_config[ second_idx_shift   ] = False
					this_config[ first_idx_shift+k+1  ] = True
					this_config[ second_idx_shift+l+1 ] = True
					# print(this_config, first_idx_shift, second_idx_shift, '==>', first_idx_shift+k+1, second_idx_shift+l+1)
					configs += [ this_config ]
			second_idx_shift += list_of_num_state[j]
		first_idx_shift += list_of_num_state[i]
	return configs



def generate_main(list_of_num_state):
	# list_of_num_state is e.g. [3,4,4,6,3] 
	# Each nnumber denotes total numbers of states for each molecular system
	#
	# This 'configs' holds all the states.
	ref_config     = generate_ref(list_of_num_state)
	# print(ref_config)
	print("NUM REFERENCE = 1")
	single_configs = generate_singles(list_of_num_state)
	# for line in single_configs:
	# 	print(line)
	print("NUM SINGLES   =",len(single_configs))
	double_configs = generate_doubles(list_of_num_state)
	# for line in double_configs:
	# 	print(line)
	print("NUM DOUBLES   =",len(double_configs))
	# 
	configs = [ref_config] + single_configs + double_configs
	#
	#
	return configs




def get_lowest_analytical(mol_obj):
    required_states = 5
    k_eff, states_configure, E_analytical, vectors = \
                analytical_solution.analytical_main(
                                                    mol_obj.list_of_masses(),
                                                    mol_obj.get_k_mat(),
                                                    required_states,
                                                    1.0, 1.0, 1.0
                                                    )
    return sorted( E_analytical )[0]



def get_rec_obj( mol_obj, rec_num_states):
	rec_obj = hamiltonian_class.RecSys_HO_hamiltonian(
                                                      mol_obj.get_k_mat(),
                                                      mol_obj.top_level_groups(),
                                                      mol_obj.top_level_coupling_mat(),
                                                      rec_num_states
                                                      )
	return rec_obj







class diag_rec_h0_mat(object):
	def __init__(self, list_of_eigvals):
		self.list_of_eigvals = copy.deepcopy(list_of_eigvals)
		self.h_diag = []
		for eigvals in self.list_of_eigvals:
			self.h_diag += eigvals

	def __call__(self, first_idx, second_idx):
		if first_idx != second_idx:
			return 0.0
		else:
			return self.h_diag[first_idx]

	def get_dimension(self):
		return len(self.h_diag)


class diag_rec_v_matrix(object):
	def __init__(self, rec_num_states, list_of_dipole_mat):
		self.rec_num_states     = copy.deepcopy(rec_num_states) 
		self.list_of_dipole_mat = copy.deepcopy(list_of_dipole_mat)
		self.V_mat_like_list    = self.get_trivial_2d_list()

	def __call__(self, mol1_num, mol2_num, mol1_idx1, mol2_idx1, mol1_idx2, mol2_idx2 ):
		return self.V_mat_like_list[mol1_num][mol2_num][mol1_idx1 * self.rec_num_states[mol2_num] + mol2_idx1][mol1_idx2 * self.rec_num_states[mol2_num] + mol2_idx2]

	def get_trivial_2d_list(self):
		num_mol = len(self.rec_num_states)
		trivial_2d_list = []
		list_idx = 0
		for i in range(num_mol):
			tmp_list = []
			for j in range(num_mol):
				if i < j:
					tmp_list +=  [ copy.deepcopy(self.list_of_dipole_mat[list_idx]) ]
					list_idx +=  1
				else:
					tmp_list += [ [] ]
			trivial_2d_list += [copy.deepcopy(tmp_list)]
		return trivial_2d_list




















