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
from state import operator, summation, matrix_term, term, state, vector_state


def T1(index_occ, index_vrt):
	'''T1 = \sum_{A I} t_I^A a+_A a_I'''
	A_op = operator('create',      'virtual',  index_vrt)
	I_op = operator('annihilate',  'occupied', index_occ)
	sum_A = summation(index_vrt, 'virtual') 
	sum_I = summation(index_occ, 'occupied')
	mat_term = matrix_term( 't', [ index_occ, index_vrt ] )
	#
	#
	list_of_sum   = [ sum_A, sum_I ]
	list_of_delta = []
	op_string     = [A_op, I_op]
	mat_term_obj  = [ mat_term ]
	vec_term      = None
	coeff = 1.0
	#
	#
	T1_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	return T1_term


def T2(index_occ1, index_occ2, index_vrt1, index_vrt2):
	'''T2 = \sum_{A<B I<J} t_{IJ AB} a+_A a_I a+_B a_J = 1/4 \sum_{A B I J} t_{IJ AB} a+_A a_I a+_B a_J'''
	A_op = operator('create',      'virtual',  index_vrt1)
	B_op = operator('create',      'virtual',  index_vrt2)

	I_op = operator('annihilate',  'occupied', index_occ1)
	J_op = operator('annihilate',  'occupied', index_occ2)

	sum_A = summation(index_vrt1, 'virtual')
	sum_B = summation(index_vrt2, 'virtual',  lower_limit=index_vrt1)

	sum_I = summation(index_occ1, 'occupied')
	sum_J = summation(index_occ2, 'occupied', lower_limit=index_occ1)

	mat_term = matrix_term( 't', [ index_occ1, index_occ2, index_vrt1, index_vrt2 ] )
	#
	#
	list_of_sum   = [ sum_A, sum_B, sum_I, sum_J ]
	list_of_delta = []
	op_string     = [A_op, B_op, I_op, J_op]  
	mat_term_obj  = [ mat_term ]
	vec_term      = None
	coeff = 1.0
	#
	#
	T2_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	return T2_term


def get_T(num_of_T_needed):
	'''T = T1 + T2'''
	# Max num of T needed is 4 (Above 4 times T, term goes to ZERO)
	# for T1, index C (virtual) and K (occupied) are used. 
	# for T2, index A,B (virtual) and I,J (occupied) are used.
	# They are capitalized to avoid conflicts with H operator parts.
	#
	list_of_T1 = []
	list_of_T2 = []
	for i in range(1,num_of_T_needed+1):
		T1_vrt_idx = 'C' + str(i)
		T1_occ_idx = 'K' + str(i)
		list_of_T1 += [ copy.deepcopy( T1(T1_occ_idx,T1_vrt_idx) ) ]
		#
		T2_vrt_idx1 = 'A' + str(i)
		T2_vrt_idx2 = 'B' + str(i)
		T2_occ_idx1 = 'I' + str(i)
		T2_occ_idx2 = 'J' + str(i)
		list_of_T2 += [ copy.deepcopy( T2(T2_occ_idx1, T2_occ_idx2, T2_vrt_idx1, T2_vrt_idx2) ) ]

	return list_of_T1, list_of_T2



if __name__ == '__main__':
	for i in range(1,5):
		list_of_T1, list_of_T2 = get_T(i)
		#
		#
		for item in list_of_T1:
			item.print_info()
		for item in list_of_T2:
			item.print_info()




