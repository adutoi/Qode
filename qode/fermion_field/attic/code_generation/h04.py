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
from state import operator, summation, matrix_term, term, state, vector_state



def h4():
	i_op = operator('create',      'occupied', 'i')
	j_op = operator('annihilate',  'occupied', 'j')

	sum_i = summation('i', 'occupied') 
	sum_j = summation('j', 'occupied') 
	sum_p = summation('p', 'occupied') 


	mat_term1 = matrix_term( 'F', [ 'i', 'j' ] )
	mat_term2 = matrix_term( 'V', [ 'i', 'p', 'j', 'p' ] )
	#
	#
	list_of_sum1  = [ sum_i, sum_j ]
	list_of_sum2  = [ sum_i, sum_j, sum_p ]

	list_of_delta = []
	op_string     = [i_op, j_op]

	mat_term_obj1 = [ mat_term1 ]
	mat_term_obj2 = [ mat_term2 ]

	vec_term = None
	coeff    =  1.0
	#
	#
	h_term1 = term( list_of_sum1, list_of_delta, op_string, mat_term_obj1, vec_term, coeff )
	h_term2 = term( list_of_sum2, list_of_delta, op_string, mat_term_obj2, vec_term, coeff )
	return [h_term1, h_term2]





