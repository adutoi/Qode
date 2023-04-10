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
from state import operator, summation, matrix_term, vector_state, term, state
from state import concatenate, compute_state
from vector_states import get_reference_state, get_single_state, get_double_state


#
#  Helper Function Only Maintain Its Inner Indentations,
#  Script-wise Indentation is handled by the computed_ci_state_to_code ( the "main" function )
#


#
# This dictionary code is very dirty but OK for now since this is not inside the automatic derivation
index_field_dict = { 'a':'vrt_orb_idx',  'b':'vrt_orb_idx',  'c':'vrt_orb_idx',  'd':'vrt_orb_idx',  'e':'vrt_orb_idx', 'f':'vrt_orb_idx', \
                     'i':'occ_orb_idx',  'j':'occ_orb_idx',  'k':'occ_orb_idx',  'l':'occ_orb_idx',  'm':'occ_orb_idx', 'n':'occ_orb_idx', \
                     'p':'occ_orb_idx',  'q':'occ_orb_idx',  'r':'occ_orb_idx',  's':'occ_orb_idx',  'g':'vrt_orb_idx', 'h':'vrt_orb_idx', \
                     'A1':'vrt_orb_idx', 'A2':'vrt_orb_idx', 'A3':'vrt_orb_idx', 'A4':'vrt_orb_idx', 'A5':'vrt_orb_idx', \
                     'B1':'vrt_orb_idx', 'B2':'vrt_orb_idx', 'B3':'vrt_orb_idx', 'B4':'vrt_orb_idx', 'B5':'vrt_orb_idx', \
                     'C1':'vrt_orb_idx', 'C2':'vrt_orb_idx', 'C3':'vrt_orb_idx', 'C4':'vrt_orb_idx', 'C5':'vrt_orb_idx', \
                     'I1':'occ_orb_idx', 'I2':'occ_orb_idx', 'I3':'occ_orb_idx', 'I4':'occ_orb_idx', 'I5':'occ_orb_idx', \
                     'J1':'occ_orb_idx', 'J2':'occ_orb_idx', 'J3':'occ_orb_idx', 'J4':'occ_orb_idx', 'J5':'occ_orb_idx', \
                     'K1':'occ_orb_idx', 'K2':'occ_orb_idx', 'K3':'occ_orb_idx', 'K4':'occ_orb_idx', 'K5':'occ_orb_idx', }


def pre_conditioner():
	return [ "num_alpha_elec  = cisd_amp_obj.get_num_alpha_elec()", 
	         "num_beta_elec   = cisd_amp_obj.get_num_beta_elec()",
             "num_spatial_orb = cisd_amp_obj.get_num_spatial_orb()", 
             "num_spin_orb    = num_spatial_orb * 2", 
             "num_occ_orb     = num_alpha_elec+num_beta_elec",
             "num_vrt_orb     = 2 * num_spatial_orb - num_alpha_elec - num_beta_elec",
             "occ_orb_idx = [ i for i in range(num_alpha_elec)] + [ i+num_spatial_orb for i in range(num_beta_elec)]",
             "vrt_orb_idx = [ i+num_alpha_elec for i in range(num_spatial_orb - num_alpha_elec) ] + [ i+num_spatial_orb+num_beta_elec for i in range(num_spatial_orb - num_beta_elec) ]" ]
   




def matrix_terms_to_initialzed_code(mat_term):
	this_symbol = mat_term.get_symbol()
	if this_symbol == 'C':
		# 'C' means CISD amplitude matrices
		indices = mat_term.get_index()
		if len(indices) == 0:
			return ["ref_amp = cisd_amp_obj.get_ref_amplitude()"]
		elif len(indices) == 2:
			return ["single_amp_mat = cisd_amp_obj.get_single_amplitude()"]
		elif len(indices) == 4:
			return ["double_amp_mat = cisd_amp_obj.get_double_amplitude()"]
		# elif len(indices) == 6:
		# 	return ["triple_amp_mat = cisd_amp_obj.get_triple_amplitude()"]
		# elif len(indices) == 8:
		# 	return ["quadruple_amp_mat = cisd_amp_obj.get_quadruple_amplitude()"]
		else:
			print("Quintuples Are Not Stored in this Scope, ERROR!")
			raise ValueError
	elif this_symbol == 't':
		# t amplitudes
		indices = mat_term.get_index()
		if len(indices) == 2:
			return ["t_single_amp_mat = t_amp_obj.get_single_amplitude()"]
		elif len(indices) == 4:
			return ["t_double_amp_mat = t_amp_obj.get_double_amplitude()"]
		else:
			raise ValueError
	else:
		return None



def vec_state_to_initialized_code(vec_state):
	return_buf = []
	occ_idx = vec_state.get_occ_idx()
	vrt_idx = vec_state.get_vrt_idx()
	#
	# Initialize the target amplitude container.
	#
	if occ_idx == None:
		return_buf += ["new_ref_amp = 0.0"]
		# Ref is done, waiting for contraction loops
	elif len(occ_idx) == 1:
		return_buf += ["new_single_amp_mat = np.zeros(cisd_amp_obj.get_single_amplitude().shape)"]
	elif len(occ_idx) == 2:
		return_buf += ["new_double_amp_mat = np.zeros(cisd_amp_obj.get_double_amplitude().shape)"]
	else:
		print("Triples Are Not Stored in this Scope, ERROR!")
		raise TypeError
	return return_buf


def sum_obj_to_code(sum_obj):
	# Take care of the upper_limit
	# if sum_obj.get_upper_limit() != None:
	# 	upper_limit = sum_obj.get_upper_limit()
	# else:
	if sum_obj.get_field() == 'occupied':
		if sum_obj.get_upper_limit() == None:
			upper_limit = "num_occ_orb"
		else:
			upper_limit = sum_obj.get_upper_limit()
	else:
		if sum_obj.get_upper_limit() == None:
			upper_limit = "num_vrt_orb"
		else:
			upper_limit = sum_obj.get_upper_limit()

	# Take care of the lower_limit
	if sum_obj.get_lower_limit() != None:
		lower_limit = sum_obj.get_lower_limit() + '+1, '
	else:
		lower_limit = ''
	return 'for ' + sum_obj.get_index() + ' in range(' + lower_limit + upper_limit + '):'


def list_of_mat_terms_to_line(mat_terms):
	return_buf = []
	for item in mat_terms:
		term_symbol = item.get_symbol()
		if term_symbol == 'F':
			indices = item.get_index()    # Dirty Code, find occ/vrt from index itself using a dictionary defined at the beginning.
			return_buf += ['F_mat[%s,%s]' %(index_field_dict[indices[0]] + '[' + indices[0] + ']', index_field_dict[indices[1]] + '[' + indices[1] + ']')  ] 

		elif term_symbol == 'V':
			indices = item.get_index() 
			return_buf += ['( -V_mat[%s * num_spin_orb + %s ,%s * num_spin_orb + %s] +  V_mat[%s * num_spin_orb + %s ,%s * num_spin_orb + %s] )' \
			                                                                             %( index_field_dict[indices[0]] + '[' + indices[0] + ']', \
				                                                                            index_field_dict[indices[1]] + '[' + indices[1] + ']', \
				                                                                            index_field_dict[indices[2]] + '[' + indices[2] + ']', \
				                                                                            index_field_dict[indices[3]] + '[' + indices[3] + ']', \
				                                                                            index_field_dict[indices[0]] + '[' + indices[0] + ']', \
				                                                                            index_field_dict[indices[1]] + '[' + indices[1] + ']', \
				                                                                            index_field_dict[indices[3]] + '[' + indices[3] + ']', \
				                                                                            index_field_dict[indices[2]] + '[' + indices[2] + ']' ) ]

		elif term_symbol == 'C':
			indices = item.get_index()
			if len(indices) == 0:
				return_buf += ['ref_amp']
			elif len(indices) == 2:
				return_buf += ['single_amp_mat[%s,%s]' %(indices[0], indices[1])]
			elif len(indices) == 4:
				return_buf += ['double_amp_mat[%s * num_vrt_orb + %s, %s * num_vrt_orb + %s]' %(indices[0], indices[2], indices[1], indices[3]) ]
			# elif len(indices) == 6:
			# 	pass # DON'T KNOW WHAT TO DO YET
			# 	# return_buf += ['quadruple_amp_mat[']
			# elif len(indices) == 8:
			# 	pass
		elif term_symbol == 't':
			indices = item.get_index()
			if len(indices) == 2:
				return_buf += ['t_single_amp_mat[%s,%s]' %(indices[0], indices[1])]
			elif len(indices) == 4:
				return_buf += ['t_double_amp_mat[%s * num_vrt_orb + %s, %s * num_vrt_orb + %s]' %(indices[0], indices[2], indices[1], indices[3]) ]
			else:
				raise ValueError
	#
	return_line = ' '
	for i in range(len(return_buf) - 1):
		return_line += return_buf[i] + ' * '
	return_line += return_buf[ -1 ] 
	#
	#
	return return_line


def generate_contraction_line(this_term, vec_state):
	mat_terms = this_term.get_mat_term()
	coeff     = this_term.get_coeff()
	occ_idx   = vec_state.get_occ_idx()
	vrt_idx   = vec_state.get_vrt_idx()
	#
	return_line = ''
	# This condition determines which matrix is holding all the numbers.
	if occ_idx == None:
		return_line += 'new_ref_amp += ' + str(coeff) + ' * ' + list_of_mat_terms_to_line(mat_terms) 
	elif len(occ_idx) == 1:
		return_line += 'new_single_amp_mat[%s,%s] += ' %(occ_idx[0], vrt_idx[0]) + str(coeff) + ' * ' + list_of_mat_terms_to_line(mat_terms) 
	elif len(occ_idx) == 2:
		return_line += 'new_double_amp_mat[%s * num_vrt_orb + %s, %s * num_vrt_orb + %s] += ' %(occ_idx[0],vrt_idx[0],occ_idx[1],vrt_idx[1]) +\
		                str(coeff) + ' * ' + list_of_mat_terms_to_line(mat_terms) 

	else:
		print("Triples Are Not Stored in this Scope, ERROR!")
		raise TypeError
	return return_line


# def generate_suffix(vec_state):
# 	return_buf  = []
# 	return_buf += ['new_ccsd_amp_obj = copy.deepcopy( ccsd_amp_obj )']
# 	return_buf += ['new_ccsd_amp_obj.clean_all_amplitude()']

# 	occ_idx   = vec_state.get_occ_idx()
# 	if occ_idx == None:
# 		return_buf += ['new_ccsd_amp_obj.update_ref_amplitude( new_ref_amp )']
# 	elif len(occ_idx) == 1:
# 		return_buf += ['new_ccsd_amp_obj.update_single_amplitude( new_single_amp_mat )']
# 	elif len(occ_idx) == 2:
# 		return_buf += ['new_ccsd_amp_obj.update_double_amplitude( new_double_amp_mat )']
# 	else:
# 		print("Triples Are Not Stored in this Scope, ERROR!")
# 		raise TypeError

# 	return_buf += ['return new_ccsd_amp_obj']
# 	return return_buf


def generate_suffix(vec_state):
	return_buf  = []
	return_buf += ['new_cisd_amp_obj = copy.deepcopy( cisd_amp_obj )']
	return_buf += ['new_cisd_amp_obj.clean_all_amplitude()']

	occ_idx   = vec_state.get_occ_idx()
	if occ_idx == None:
		return_buf += ['new_cisd_amp_obj.update_ref_amplitude( new_ref_amp )']
	elif len(occ_idx) == 1:
		return_buf += ['new_cisd_amp_obj.update_single_amplitude( new_single_amp_mat )']
	elif len(occ_idx) == 2:
		return_buf += ['new_cisd_amp_obj.update_double_amplitude( new_double_amp_mat )']
	else:
		print("Triples Are Not Stored in this Scope, ERROR!")
		raise TypeError

	return_buf += ['return new_cisd_amp_obj']
	return return_buf


def add_lines_to_buffer(return_buf, lines, prefix='\t', suffix='\n'):
	# lines is a list of strings
	# Add '\t' to the front and '\n' to the end
	for line in lines:
		return_buf += prefix + line + suffix
	return return_buf


def computed_cc_state_to_code( state_obj, func_name ):
	'''Take a computed CI state, generate a multiline long string of Python code'''
	return_buf =  "def " + func_name + "( cisd_amp_obj, t_amp_obj , F_mat , V_mat ):\n"
	lines = pre_conditioner()
	return_buf = add_lines_to_buffer(return_buf, lines)
	terms = state_obj.get_list_of_terms()
	if len(terms) > 0:
		# Non-Zero returns, keep going ...
		first_mat_terms = terms[0].get_mat_term()
		lines = []
		for each_term in first_mat_terms:
			return_term = matrix_terms_to_initialzed_code(each_term)
			if return_term != None:
				lines += return_term
		return_buf = add_lines_to_buffer(return_buf, lines)

		first_vec_term = terms[0].get_vec_state()
		return_buf = add_lines_to_buffer(return_buf, vec_state_to_initialized_code( first_vec_term ) )

		# 2) Read summations and translate to code.
		for each_term in terms:
			list_of_sum = each_term.get_sum_objs()
			lines = []
			for each_sum in list_of_sum:
				lines += [sum_obj_to_code(each_sum)]

			counts = 0
			for line in lines:
				return_buf = add_lines_to_buffer( return_buf, [line], '\t' * (counts+1), '\n')
				counts += 1
			#
			# 3) Read vector_state, coeff, matrix terms, accumulate them into new_states ( contraction line )
			return_buf  = add_lines_to_buffer(return_buf, [ generate_contraction_line(each_term, each_term.get_vec_state()) ], '\t' * (counts+1) )

		return_buf = add_lines_to_buffer(return_buf, generate_suffix(first_vec_term) )
	print(return_buf)
	return return_buf




if __name__ == "__main__":
	import pickle
	ham_states = pickle.load( open("non_zero_ccsd_states.p","rb") )
	
	# for item in ham_states:
	# 	item.print_info()


	output_buf = ''
	for i in range(len(ham_states)):
		ham_states[i].print_info()
		output_buf += computed_cc_state_to_code( ham_states[i], 'h' + str(i) ) + '\n'
	print(output_buf)
	output_buf = 'import copy\nimport numpy as np\n' + output_buf
	f = open('ccsd_code_normal_ordered.py','w')
	f.write(output_buf)
	f.close()


