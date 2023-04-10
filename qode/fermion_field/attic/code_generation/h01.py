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



def h1():
	'''H_1 =  '''
	a_op = operator('create',     'virtual',  'a')
	i_op = operator('annihilate', 'occupied', 'i')
	sum_i = summation('i', 'occupied') 
	sum_a = summation('a', 'virtual')
	sum_p = summation('p', 'occupied')
	mat_term1 = matrix_term( 'F', [ 'a', 'i' ] )
	mat_term2 = matrix_term( 'V', [ 'a', 'p', 'i','p'])
	#
	#
	list_of_sum1  = [ sum_a, sum_i ]
	list_of_sum2  = [ sum_a, sum_i, sum_p ]
	list_of_delta = []
	op_string     = [a_op, i_op]
	mat_term_obj1 = [ mat_term1 ]
	mat_term_obj2 = [ mat_term2 ]
	vec_term      = None
	coeff = 1.0
	#
	#
	h_term1 = term( list_of_sum1, list_of_delta, op_string, mat_term_obj1, vec_term, coeff )
	h_term2 = term( list_of_sum2, list_of_delta, op_string, mat_term_obj2, vec_term, coeff )
	return [h_term1, h_term2]


