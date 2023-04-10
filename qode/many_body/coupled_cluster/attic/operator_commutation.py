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
import qode

class operator_list(object):
	def __init__(self, list_of_operator_list, list_of_coeffs):
		#
		# list_of_operators = [ op1, op2, op3, ... ... ]
		# list_of_coeffs    = [   1,  -1,  -1, ... ... ]
		#
		self.operators = copy.deepcopy( list_of_operator_list )
		self.coeffs    = copy.deepcopy( list_of_coeffs )
		# Sanity check: len of operators == len of coeffs or raise error
		if len(self.operators) != len(self.coeffs):
			print("num of operator list NOT EQUAL to num coeffs")
			raise AssertionError
		
		#
	# def coeff_change_sign(self, op_list_index):
		# self.coeffs[op_list_index]  = -self.coeffs[op_list_index] 
	def get_operators(self):
		return copy.deepcopy(self.operators)
	def get_coeffs(self):
		return copy.deepcopy( self.coeffs )
	def update_operators(self, new_operators):
		self.operators = copy.deepcopy(new_operators)
	def update_coeffs(self, new_coeffs):
		self.coeffs = copy.deepcopy(new_coeffs)
	def combine_coeffs(self):
		new_coeffs    = []
		new_operators = []
		size = len(self.coeffs)
		for i in range(size):
			if self.operators[i] not in new_operators:
				new_operators += [ copy.deepcopy(self.operators[i]) ]
				new_coeffs    += [ copy.deepcopy(self.coeffs[i]) ]
			else:
				j = 0
				checking = True
				new_size = len(new_coeffs)
				while j < new_size and checking:
					if new_operators[j] == self.operators[i]:
						new_coeffs[j] += self.coeffs[i]
						checking = False
					else:
						j += 1
				if j == new_size:
					print("Found in list but failed to add the coeff, horrible error!")
					raise AssertionError
		self.operators = copy.deepcopy(new_operators)
		self.coeffs    = copy.deepcopy(new_coeffs)

	def print_info(self):
		size = len(self.coeffs)
		for i in range(size):
			print("Coeff =",self.coeffs[i], "OPs =", self.operators[i])





def commute( op_list_obj1, op_list_obj2 ):
	#
	# [A, B] = AB - BA
	#  return AB - BA as one operator_list object with a top level list of two element [AB, -BA].
	#
	# print(" ")
	# print("Commute", op_list_obj1, op_list_obj2)
	# print("------------------------------------------------------------------------------------------------")
	#
	first_ops     =  op_list_obj1.get_operators()
	first_coeffs  =  op_list_obj1.get_coeffs()
	second_ops    =  op_list_obj2.get_operators()
	second_coeffs =  op_list_obj2.get_coeffs()
	#
	size1 = len(first_coeffs)
	size2 = len(second_coeffs)
	#
	new_list   = []
	new_coeffs = []
	#
	# AB part:
	for i in range(size1):
		for j in range(size2):
			# print(i,j)
			new_list   += [ first_ops[i] + second_ops[j] ]
			new_coeffs += [ first_coeffs[i] * second_coeffs[j] ]
	#
	# -BA part:
	for i in range(size2):
		for j in range(size1):
			# print(i,j)
			new_list   += [ second_ops[i] + first_ops[j] ]
			new_coeffs += [ -second_coeffs[i] * first_coeffs[j] ]
	#
	#	
	new_op_list_obj = operator_list( new_list, new_coeffs )
	new_op_list_obj.combine_coeffs()
	return new_op_list_obj





def scale_py_list(py_list, value):
	new_list = [ item * value for item in py_list ]
	return new_list





def omega0_ops(H,T):
	# \Omega_\mu^(0) (t^(n))   = < \mu | exp(-T^(n)) H exp(T^(n)) | HF >
	# exp(-T^(n)) H exp(T^(n)) = H + [H,T] + 1/2 [[H,T],T] + 1/6 [[[H,T],T],T] + 1/24 [[[[H,T],T],T],T]
	# Return all the permuated operators.
	#
	# Check if H, T are primitive functions.
	#
	if isinstance(H, operator_list):
		# print("H is an OP List")
		H_op = copy.deepcopy(H)
	else:
		# print("H is a primitive function")
		H_op = operator_list( [[H]], [1.0] )

	if isinstance(T, operator_list):
		# print("T is an OP List")
		T_op = copy.deepcopy(T)
	else:
		# print("T is a primitive function")
		T_op = operator_list( [[T]], [1.0] )
	# H_op.print_info()
	# T_op.print_info()
	#
	temp_op_list = H_op
	operators = [ copy.deepcopy(temp_op_list) ]

	# print(temp_op_list)
	# print(operators)
	# temp_op_list.print_info()
	for i  in range(4):
		temp_op_list = commute(temp_op_list, T_op)
		operators += [ copy.deepcopy(temp_op_list) ]

	operators[2].update_coeffs( scale_py_list( operators[2].get_coeffs(), 0.5      ) )
	operators[3].update_coeffs( scale_py_list( operators[3].get_coeffs(), 1.0/6.0  ) )	
	operators[4].update_coeffs( scale_py_list( operators[4].get_coeffs(), 1.0/24.0 ) )
	
	return operators





def flatten_operator_list(commutator_list):
	flat_op_list = []
	flat_coeffs  = []
	#
	#
	for operator_list_obj in commutator_list:
		operator = operator_list_obj.get_operators()
		coeffs   = operator_list_obj.get_coeffs()
		for i in range(len(coeffs)):
			flat_op_list += [ operator[i] ]
			flat_coeffs  += [ coeffs[i]   ]
	return operator_list(flat_op_list, flat_coeffs)






def HBar_term_on_HF(op_sequence, coefficient, occ_strings):
	# Each operator in current term, accumulate them into return_state
	#
	return_state = qode.fermion_field.state.state(occ_strings,0)
	for op in reversed(op_sequence):  return_state = op(return_state) 
	return_state.scale(coefficient)
	#
	return return_state

def HBar_on_HF(H, T, occ_strings, resources=qode.util.parallel.resources(1), textlog=print):
	#
	omega0_list = omega0_ops(H, T)
	flat_omega0_list = flatten_operator_list(omega0_list)
	del omega0_list
	operators = flat_omega0_list.get_operators()
	coeffs    = flat_omega0_list.get_coeffs()
	#
	input_list = [ (o, c, occ_strings) for o,c in zip(operators,coeffs) ]
	textlog("Running on {} cores.".format(resources.n_cores))
	list_return_states = qode.util.parallelize_task(HBar_term_on_HF,input_list,resources.n_cores)
	#
	# Combine all states.
	return_state = qode.fermion_field.state.state(occ_strings)
	for state_term in list_return_states:  return_state.increment(state_term)
	#
	return return_state
