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


# def T():
# 	# trivial, a function name needed to permuate
# 	pass

# def H():
# 	# trivial, a function name needed to permuate
# 	pass

# def get_all_operator_strings():
# 	omega0_list      = operator_commutation.omega0_ops(H, T)
# 	flat_omega0_list = operator_commutation.flatten_operator_list(omega0_list)
# 	# flat_omega0_list.print_info()
# 	ops_char = []
# 	for op_string in flat_omega0_list.get_operators():
# 		# print(op_string)
# 		current_list = []
# 		for op in op_string:
# 			if op == H:
# 				current_list += ['H']
# 				# print('H',end='  ')
# 			else:
# 				current_list += ['T']
# 				# print('T',end='  ')
# 		ops_char += [copy.deepcopy(current_list)]
# 	# print(ops_char)
# 	coeffs = flat_omega0_list.get_coeffs()
# 	# print(coeffs)
# 	# for i in range(len(coeffs)):
# 	# 	print("%.4f\t%s" %(coeffs[i], ops_char[i]) ) 	
# 	return coeffs, ops_char

# def get_all_terms():
# 	coeffs, ops_char_list = get_all_operator_strings()
# 	for i in range(len(coeffs)):
# 		print("%.4f\t%s" %(coeffs[i], ops_char_list[i]) ) 
# 	whole_ham_list = get_all_ham_char(ops_char_list)

# if __name__ == "__main__":
# 	coeffs, ops_char_list = get_all_operator_strings()
# 	for i in range(len(coeffs)):
# 		print("%.6f\t%s" %(coeffs[i], ops_char_list[i]) )



#
#  Generate All Commutator Terms: (This is for convenience)
#
#
#  1	['H']
#
#  1	['H', 'T']
# -1	['T', 'H']
#
#  1/2	['H', 'T', 'T']
# -1	['T', 'H', 'T']
#  1/2	['T', 'T', 'H']
#
#  1/6	['H', 'T', 'T', 'T']
# -1/2	['T', 'H', 'T', 'T']
#  1/2	['T', 'T', 'H', 'T']
# -1/6	['T', 'T', 'T', 'H']
#
#  1/24	['H', 'T', 'T', 'T', 'T']
# -1/6	['T', 'H', 'T', 'T', 'T']
#  1/4	['T', 'T', 'H', 'T', 'T']
# -1/6	['T', 'T', 'T', 'H', 'T']
#  1/24	['T', 'T', 'T', 'T', 'H']
#
#

import copy
import operator_commutation
from state import concatenate, state, compute_state
from vector_states import get_reference_state
from t_operator import get_T

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


ham_dict = {'h1':h1(), 'h2':h2(), 'h3':h3(), 'h4':h4(), 'h5':h5(), 'h6':h6(), 'h7':h7(),\
                'h8':h8(), 'h9':h9(), 'h10':h10(), 'h11':h11(), 'h12':h12(), 'h13':h13()}


def bch_nth_order(order_n):
	if order_n == 0:
		return [1.0], [['H']]
	elif order_n == 1:
		return [1.0, -1.0], [['H', 'T'], ['T', 'H']]
	elif order_n == 2:
		return [1.0/2.0, -1.0, 1.0/2.0], [['H', 'T', 'T'], ['T', 'H', 'T'], ['T', 'T', 'H']]
	elif order_n == 3:
		return [1.0/6.0, -1.0/2.0, 1.0/2.0, -1/6.0], [['H', 'T', 'T', 'T'], ['T', 'H', 'T', 'T'], ['T', 'T', 'H', 'T'], ['T', 'T', 'T', 'H']]
	elif order_n == 4:
		return [1.0/24.0,-1.0/6.0,1.0/4.0,-1.0/6.0,1.0/24.0], [['H', 'T', 'T', 'T', 'T'], ['T', 'H', 'T', 'T', 'T'], ['T', 'T', 'H', 'T', 'T'], ['T', 'T', 'T', 'H', 'T'], ['T', 'T', 'T', 'T', 'H']]
	else:
		return [],[]


def get_H():
	# These are term class objects, not functions
	return ['h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','h13']


def permuate_op_list(list_term1, list_term2):
	new_list = []
	for each_term1 in list_term1:
		for each_term2 in list_term2:
			new_list += [ each_term1 + ' ' + each_term2 ]
	return new_list




def operators_to_terms( list_of_ops ):
	# print(list_of_ops)
	if len(list_of_ops) == 1:
		# Only H is possible
		if list_of_ops[0] == 'H':
			return get_H()
		else:
			print("Stupid Bug, only H is possible.")
			raise AssertionError
	elif len(list_of_ops) > 1:
		list_of_separate_op = []
		for op_char in list_of_ops:
			if op_char == 'H':
				list_of_separate_op += [ get_H() ]
			elif op_char == 'T':
				list_of_separate_op += [ ['T1', 'T2'] ]
			else:
				print("UNKNOWN OPERATOR TYPE USED, PLEASE DOUBLE CHECK.")
				raise AssertionError

		# print(list_of_separate_op)
		list_of_terms = permuate_op_list(list_of_separate_op[0], list_of_separate_op[1])
		for i in range(2, len(list_of_separate_op)):
			list_of_terms = permuate_op_list(list_of_terms, list_of_separate_op[i])
		return list_of_terms
	else:
		print("THIS IS NOT EVEN POSSIBLE")
		raise AssertionError


def str_to_terms(op_in_str):
	ops = op_in_str.split()
	num_T_needed = len(ops) - 1
	T1_list, T2_list = get_T(num_T_needed)
	t_count = 0

	if ops[0].startswith('h'):
		op_terms_list = copy.deepcopy(ham_dict[ops[0]])
	elif ops[0] == 'T1':
		op_terms_list = [ copy.deepcopy(T1_list[t_count]) ]
		t_count += 1
	elif ops[0] == 'T2':
		op_terms_list = [ copy.deepcopy(T2_list[t_count]) ]
		t_count += 1
	else:
		raise AssertionError


	for i in range(1, len(ops)):
		op = ops[i]
		if op.startswith('h'):
			next_op = copy.deepcopy(ham_dict[op])
		elif op == 'T1':
			next_op = [ copy.deepcopy(T1_list[t_count]) ]
			t_count += 1
		elif op == 'T2':
			next_op = [ copy.deepcopy(T2_list[t_count]) ]
			t_count += 1
		else:
			raise AssertionError		
		temp_terms = []
		for item_in_list in op_terms_list:
			for new_item in next_op:
				temp_terms += [ copy.deepcopy( concatenate(item_in_list, new_item) ) ]
		op_terms_list = copy.deepcopy(temp_terms)
	return op_terms_list



def compute_this_state(state_obj):
	return compute_state(state_obj)

def get_full_term_list():
	import time
	import pickle
	from multiprocessing import Pool
	#
	#
	all_state_list = []
	max_order = 4   # inclusive
	for i in range(max_order+1):
		coeffs, op_terms = bch_nth_order(i)
		print(coeffs, op_terms)
		for j in range(len(coeffs)):
			same_coeff_terms = operators_to_terms(op_terms[j])
			# for line in same_coeff_terms:
			# 	print(line)
			for this_string in same_coeff_terms:
				terms_in_list = str_to_terms(this_string)
				state_term_list = []
				for item in terms_in_list:
					each_new_term = concatenate(item, get_reference_state())
					each_new_term.update_coeff( each_new_term.get_coeff() * coeffs[j] )
					state_term_list +=[ copy.deepcopy(each_new_term) ]
				all_state_list += [copy.deepcopy(state(state_term_list))]
	# for each_term in all_state_list:
	# 	each_term.print_info()
	print("TOTAL %d TERMS" %(len(all_state_list)))
	# pickle.dump(all_state_list, open('all_ccsd_states.p','wb'))
	# computed_states = []
	# for each_term in all_state_list:
	# 	each_computed_state = compute_state(each_term)
	# 	each_computed_state.print_info()
	# 	computed_states += [copy.deepcopy(each_computed_state)]
	myPool = Pool(24)
	print("CCSD Equation Derivation Starts at",time.ctime())
	t1 = time.time()
	computed_states = myPool.map(compute_this_state, all_state_list)
	t2 = time.time()
	myPool.terminate()
	print("Derivation Done at", time.ctime())
	print("Derivation Duration = %d sec" %(t2-t1))

	surviving_terms = []

	for each_state in computed_states:
		if len(each_state.get_list_of_terms()) > 0:
			surviving_terms += [ copy.deepcopy(each_state) ]
	
	print("%d SURVIVING TERMS" %(len(surviving_terms)))
	pickle.dump(surviving_terms, open('non_zero_ccsd_states.p','wb'))



if __name__ == "__main__":
	# for i in range(5):
	# 	c,op = bch_nth_order(i)
	# 	print("++++++++++++++++++++++++++++++++++")
	# 	print(c, op)

	# for item in (get_H()):
	# 	print("++++++++++++++++++++++++++++++++++")
	# 	for each_term in item:
	# 		each_term.print_info()

	# for i in range(1,5):
	# 	t1, t2 = get_T(i)
	# 	print("++++++++++++++++++++++++++++++++++")
	# 	for each_term in t1:
	# 		each_term.print_info()
	# 	for each_term in t2:
	# 		each_term.print_info()
	get_full_term_list()
































