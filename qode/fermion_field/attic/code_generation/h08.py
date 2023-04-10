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
# import state

from state import operator, summation, matrix_term, vector_state, term, state




def h8():
	'''H_8 = \sum_{i j k l} 1/4 V_{i j k l} i+ j+ k l '''
	i_op = operator('create',      'occupied', 'i')
	j_op = operator('create',      'occupied', 'j')
	k_op = operator('annihilate',  'occupied', 'k')
	l_op = operator('annihilate',  'occupied', 'l')

	sum_i = summation('i', 'occupied') 
	sum_j = summation('j', 'occupied') 
	sum_k = summation('k', 'occupied') 
	sum_l = summation('l', 'occupied') 
	
	mat_term = matrix_term( 'V', [ 'i', 'j', 'k', 'l' ] )
	#
	#
	list_of_sum   = [ sum_i, sum_j, sum_k, sum_l ]
	list_of_delta = []
	op_string     = [i_op, j_op, k_op, l_op]
	mat_term_obj  = [ mat_term ]
	vec_term      = None
	coeff         = 0.25
	#
	#
	h_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	return [h_term]


if __name__ == "__main__":
	from vector_states import get_double_state
	from state import concatenate, compute_state

	# double state
	# print("============================================")
	double_state = get_double_state()
	ham = h8()
	terms = []
	for item in ham:
		terms += [concatenate(item, double_state)]
	ham_double  = state( terms )
	ham_double_return = compute_state(ham_double)
	print("h|D> =", ham_double_return.return_info() )
	cont_double = ham_double_return.get_contraction_relation()
	print("Equations:")
	for item in cont_double:
		item.print_info()
	# print("============================================")
	# returned_code = computed_ci_state_to_code(ham_double_return, 'h' + str(num_i+1) + '_double')
	# print(returned_code)
	# print("============================================")
	
