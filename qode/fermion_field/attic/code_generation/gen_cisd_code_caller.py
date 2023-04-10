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
from vector_states import get_reference_state, get_single_state, get_double_state
from state import concatenate, compute_state
from gen_cisd_code import computed_ci_state_to_code


from h01 import h1
from h02 import h2
from h03 import h3
from h04 import h4
from h05 import h5
from h06 import h6
from h07 import h7
from h08 import h8
from h09 import h9
from h10 import h10
from h11 import h11
from h12 import h12
from h13 import h13
# from h14 import h14



def get_ham_ref(ham, num_i):
	# reference state
	# print("============================================")
	ref_state = get_reference_state()
	terms = []
	for item in ham:
		terms += [concatenate(item, ref_state)]
	ham_ref = state( terms )
	ham_ref_return = compute_state(ham_ref)
	# print("h|R> =", ham_ref_return.return_info() )
	# cont_ref = ham_ref_return.get_contraction_relation()
	# print("Equations:")
	# for item in cont_ref:
	# 	item.print_info()
	returned_code = computed_ci_state_to_code(ham_ref_return, 'h' + str(num_i+1) + '_ref')
	# print("============================================")
	# print(returned_code)
	# print("============================================")
	return returned_code

def get_ham_single(ham, num_i):
	# single state
	# print("============================================")
	single_state = get_single_state()
	terms = []
	for item in ham:
		terms += [concatenate(item, single_state)]
	ham_single = state( terms )
	ham_single_return = compute_state(ham_single)
	# print("h|S> =", ham_single_return.return_info() )
	# cont_single = ham_single_return.get_contraction_relation()
	# print("Equations:")
	# for item in cont_single:
	# 	item.print_info()
	# print("============================================")
	returned_code = computed_ci_state_to_code(ham_single_return, 'h' + str(num_i+1) + '_single')
	# print(returned_code)
	# print("============================================")
	return returned_code

def get_ham_double(ham, num_i):
	# double state
	# print("============================================")
	double_state = get_double_state()
	terms = []
	for item in ham:
		terms += [concatenate(item, double_state)]
	ham_double  = state( terms )
	ham_double_return = compute_state(ham_double)
	# print("h|D> =", ham_double_return.return_info() )
	# cont_double = ham_double_return.get_contraction_relation()
	# print("Equations:")
	# for item in cont_double:
	# 	item.print_info()
	# print("============================================")
	returned_code = computed_ci_state_to_code(ham_double_return, 'h' + str(num_i+1) + '_double')
	# print(returned_code)
	# print("============================================")
	return returned_code

# def get_ham_triple(ham, num_i):
# 	# double state
# 	# print("============================================")
# 	triple_state = get_triple_state()
# 	terms = []
# 	for item in ham:
# 		terms += [concatenate(item, triple_state)]
# 	ham_triple = state( terms )
# 	ham_triple_return = compute_state(ham_triple)
# 	# print("h|D> =", ham_double_return.return_info() )
# 	# cont_double = ham_double_return.get_contraction_relation()
# 	# print("Equations:")
# 	# for item in cont_double:
# 	# 	item.print_info()
# 	# print("============================================")
# 	returned_code = computed_ci_state_to_code(ham_triple_return, 'h' + str(num_i+1) + '_triple')
# 	# print(returned_code)
# 	# print("============================================")
# 	return returned_code


# def get_ham_quadruple(ham, num_i):
# 	# double state
# 	# print("============================================")
# 	quadruple_state = get_quadruple_state()
# 	terms = []
# 	for item in ham:
# 		terms += [concatenate(item, quadruple_state)]
# 	ham_quadruple  = state( terms )
# 	ham_quadruple_return = compute_state(ham_quadruple)
# 	# print("h|D> =", ham_double_return.return_info() )
# 	# cont_double = ham_double_return.get_contraction_relation()
# 	# print("Equations:")
# 	# for item in cont_double:
# 	# 	item.print_info()
# 	# print("============================================")
# 	returned_code = computed_ci_state_to_code(ham_quadruple_return, 'h' + str(num_i+1) + '_quadruple')
# 	# print(returned_code)
# 	# print("============================================")
# 	return returned_code

if __name__ == "__main__":
	# reference state
	list_of_h = [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13] 


	output_buf = "import numpy as np\nimport copy\n"
	i = 0
	for item in list_of_h:
		# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		# print(" ")
		print(item.__name__)
		h_ref      = get_ham_ref( item(), i )
		h_single   = get_ham_single( item(), i )
		h_double   = get_ham_double( item(), i )
		# h_triple   = get_ham_triple( item(), i )
		# h_quaduple = get_ham_quadruple( item(), i )

		if h_ref != None:
			output_buf += '\n' + h_ref
		if h_single != None:
			output_buf += '\n' + h_single
		if h_double != None:
			output_buf += '\n' + h_double
		# if h_triple != None:
		# 	output_buf += '\n' + h_triple
		# if h_quaduple != None:
		# 	output_buf += '\n' + h_quaduple
		i += 1
		# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		

	# print(output_buf)
	f = open("normal_order_h_cisd_py_generated.py", 'w')
	f.write(output_buf)
	f.close()


















