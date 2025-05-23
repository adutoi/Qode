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




def get_reference_state():
	coeff_ref      =  matrix_term( 'C', [] )
	list_of_sum    =  []
	list_of_delta  =  [] 
	op_string      =  []
	mat_term_obj   =  [ coeff_ref ]
	vec_term       =  vector_state()
	coeff = 1.0
	#
	#
	ref_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	#
	return ref_term



def get_single_state():
	sum_m = summation('m', 'occupied')
	sum_e = summation('e', 'virtual' )
	coeff_single = matrix_term( 'C', [ 'm', 'e' ] )
	#
	list_of_sum   = [ sum_m, sum_e ]
	list_of_delta = []
	op_string     = []
	mat_term_obj  = [ coeff_single ]
	vec_term      = vector_state(['m'], ['e'])
	coeff = 1.0
	#
	#
	single_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	#
	return single_term


def get_double_state():
	sum_m = summation('m', 'occupied')
	sum_n = summation('n', 'occupied', lower_limit = 'm') 

	sum_e = summation('e', 'virtual' )
	sum_f = summation('f', 'virtual' , lower_limit = 'e')
	coeff_double = matrix_term( 'C', [ 'm', 'n', 'e', 'f' ] )
	#
	list_of_sum   = [ sum_m, sum_n, sum_e, sum_f ]
	list_of_delta = []
	op_string     = []
	mat_term_obj  = [ coeff_double ]
	vec_term      = vector_state(['m','n'], ['e','f'])
	coeff = 1.0
	#
	#
	double_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	#
	return double_term


# Triple    is  \sum_{m<n<r, e<f<g}     C_{m,n,r}^{e,f,g}     | \Psi_{m,n,r}^{e,f,g} >
# Quadruple is  \sum_{m<n<r<s, e<f<g<h} C_{m,n,r,s}^{e,f,g,h} | \Psi_{m,n,r,s}^{e,f,g,h} >


# def get_triple_state():
# 	sum_m = summation('m', 'occupied') #, upper_limit = 'n')  # This upper_limit doesn't really matter, none of these limits matter at all.
# 	sum_n = summation('n', 'occupied', lower_limit = 'm')  # b/c these outer loops are generated somewhere else according to final vector_state
# 	sum_r = summation('r', 'occupied', lower_limit = 'n')

# 	sum_e = summation('e', 'virtual' ) #, upper_limit = 'f')
# 	sum_f = summation('f', 'virtual' , lower_limit = 'e')
# 	sum_g = summation('g', 'virtual' , lower_limit = 'f') 
# 	coeff_triple = matrix_term( 'C', [ 'm', 'n', 'r', 'e', 'f', 'g' ] )
# 	#
# 	list_of_sum = [ sum_m, sum_n, sum_r, sum_e, sum_f, sum_g ]
# 	list_of_delta = []
# 	op_string     = []
# 	mat_term_obj  = [ coeff_triple ]
# 	vec_term      = vector_state(['m','n','r'], ['e','f','g'])
# 	coeff         = 1.0
# 	#
# 	triple_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
# 	return triple_term


# def get_quadruple_state():
# 	sum_m = summation('m', 'occupied') #, upper_limit = 'n')  # This upper_limit doesn't really matter, none of these limits matter at all.
# 	sum_n = summation('n', 'occupied', lower_limit = 'm')  # b/c these outer loops are generated somewhere else according to final vector_state
# 	sum_r = summation('r', 'occupied', lower_limit = 'n')
# 	sum_s = summation('s', 'occupied', lower_limit = 'r')

# 	sum_e = summation('e', 'virtual' ) #, upper_limit = 'f')
# 	sum_f = summation('f', 'virtual' , lower_limit = 'e')
# 	sum_g = summation('g', 'virtual' , lower_limit = 'f') 
# 	sum_h = summation('h', 'virtual' , lower_limit = 'g') 
# 	coeff_triple = matrix_term( 'C', [ 'm', 'n', 'r', 's', 'e', 'f', 'g', 'h' ] )
# 	#
# 	list_of_sum = [ sum_m, sum_n, sum_r, sum_s, sum_e, sum_f, sum_g, sum_h ]
# 	list_of_delta = []
# 	op_string     = []
# 	mat_term_obj  = [ coeff_triple ]
# 	vec_term      = vector_state(['m','n','r','s'], ['e','f','g','h'])
# 	coeff         = 1.0
# 	#
# 	quadruple_term = term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
# 	return quadruple_term




# if __name__ == "__main__":
# 	r = get_reference_state()
# 	s = get_single_state()
# 	d = get_double_state()
# 	t = get_triple_state()
# 	q = get_quadruple_state()

# 	r.print_info()
# 	s.print_info()
# 	d.print_info()
# 	t.print_info()
# 	q.print_info()
	








