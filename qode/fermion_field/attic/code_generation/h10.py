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



def h10():
	# a+ i+ b c
	a_op = operator('create',     'virtual',  'a')
	i_op = operator('create',     'occupied', 'i')
	b_op = operator('annihilate', 'virtual',  'b')
	c_op = operator('annihilate', 'virtual',  'c')

	sum_i = summation('i', 'occupied') 
	sum_a = summation('a', 'virtual') 
	sum_b = summation('b', 'virtual')
	sum_c = summation('c', 'virtual')

	mat_term = matrix_term( 'V', [ 'a', 'i', 'b', 'c' ] )
	#
	#
	list_of_sum   = [ sum_i, sum_a, sum_b, sum_c ]
	list_of_delta = []
	op_string     = [a_op, i_op, b_op, c_op]
	mat_term_obj  = [ mat_term ]
	vec_term      = None
	coeff         = 0.5
	#
	#
	h_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	return [h_term]