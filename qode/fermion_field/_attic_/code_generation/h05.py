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



def h5():
	'''H_5 = '''
	a_op = operator('create',      'virtual', 'a')
	b_op = operator('create',      'virtual', 'b')
	c_op = operator('annihilate',  'virtual', 'c')
	d_op = operator('annihilate',  'virtual', 'd')
	sum_a = summation('a', 'virtual') 
	sum_b = summation('b', 'virtual')  
	sum_c = summation('c', 'virtual') 
	sum_d = summation('d', 'virtual') 
	
	mat_term = matrix_term( 'V', [ 'a','b','c','d' ] )
	#
	#
	list_of_sum   = [ sum_a, sum_b, sum_c, sum_d ]
	list_of_delta = []
	op_string     = [a_op, b_op, c_op, d_op]
	mat_term_obj  = [ mat_term ]
	vec_term      = None
	coeff = 0.25
	#
	#
	h_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	return [h_term]

