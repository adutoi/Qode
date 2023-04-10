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


class operator(object):
	"""primitive creation and annihilation operators"""
	def __init__(self, op_name, op_acting_space, op_idx):
		# op_name = 'create' or 'annihilate'
		# op_acting_space = 'occupied' or 'virtual'
		if op_name == 'create' or op_name == 'annihilate':
			self.op_name = op_name
		else:
			raise TypeError
		if op_acting_space == 'occupied' or op_acting_space == 'virtual':
			self.op_acting_space = op_acting_space
		else:
			raise TypeError
		self.op_idx = op_idx
	def name(self):
		return self.op_name
	def space(self):
		return self.op_acting_space
	def op_index(self):
		return self.op_idx
	def update_index(self, old_index, new_index):
		if old_index == self.op_idx:
			self.op_idx = new_index
	# quick info
	def print_info(self):
		print(self.op_name, self.op_acting_space, self.op_idx)
	def return_info(self):
		if self.op_name == 'create':
			output = ' a+_'
		else:
			output = ' a_'
		return output + str(self.op_idx) + ' '



class delta_func(object):
	"""delta functions used in the derivation"""
	def __init__(self, idx1, idx2, orb1_space, orb2_space):
		self.idx1 = idx1
		self.idx2 = idx2
		self.orb1_space = orb1_space
		self.orb2_space = orb2_space
	def indices(self):
		return self.idx1, self.idx2
	def name(self):
		return 'delta'
	def space(self):
		return self.orb1_space, self.orb2_space
	def update_index(self, old_index, new_index):
		if old_index == self.idx1:
			self.idx1 = new_idx
		if old_index == self.idx2:
			self.idx2 = new_idx
	# quick info
	def print_info(self):
		print("delta_{%s,%s}" %(self.idx1, self.idx2))
	def return_info(self):
		return " delta_{%s,%s} " %(self.idx1, self.idx2)
		


class summation(object):
	"""Summation Obj Class"""
	def __init__(self, summing_idx, index_field, lower_limit=None, upper_limit=None):
		# summing_idx   is a char
		# lower_limit, upper_limit => char, to store info like: i<j or a>b ...
		self.summing_idx  =  summing_idx
		self.lower_limit  =  lower_limit 
		self.upper_limit  =  upper_limit
		#
		# index_field is either 'occupied' or 'virtual'
		if index_field == 'occupied' or index_field == 'virtual':
			self.index_field  =  index_field
		else:
			raise TypeError
	# Readers:
	def get_index(self):
		return self.summing_idx
	def get_field(self):
		return self.index_field
	def get_lower_limit(self):
		return self.lower_limit
	def get_upper_limit(self):
		return self.upper_limit
	# Writers:
	def update_index(self, old_index, new_index):
		if old_index == self.summing_idx:
			self.summing_idx = new_index
	def update_limits(self, old_index, new_index):
		if old_index == self.lower_limit:
			self.lower_limit = new_index
		elif old_index == self.upper_limit:
			self.upper_limit = new_index
	def update_field(self, new_field):
		if new_field == 'occupied' or new_field == 'virtual':
			self.index_field  =  new_field
		else:
			raise TypeError
	def update_lower_limit(self, new_limit):
		self.lower_limit = new_limit
	def update_upper_limit(self, new_limit):
		self.upper_limit = new_limit
	# comparisons
	def is_identical_to(self, other_sum_obj):
		if self.summing_idx == other_sum_obj.get_index() and self.index_field == other_sum_obj.get_field() and\
		self.lower_limit == other_sum_obj.get_lower_limit() and self.upper_limit == other_sum_obj.get_upper_limit():
			return True
		else:
			return False
	# info
	def return_info(self):
		if self.lower_limit == None:
			lower_limit_str = ''
		else:
			lower_limit_str = self.summing_idx + ' > ' + self.lower_limit

		if self.upper_limit == None:
			upper_limit_str = ''
		else:
			upper_limit_str = self.summing_idx + ' < ' + self.upper_limit
		
		if lower_limit_str != '':
			output_str = ', ' + lower_limit_str 
		else:
			output_str = ''
			
		if upper_limit_str != '':
			output_str += ', ' + upper_limit_str 
		else:
			pass
		return "\\sum_{ " + self.summing_idx + output_str + " } " 
	def print_info(self):
		print(self.return_info())





class matrix_term(object):
	"""where the matrix elements are represented:
	            mat_symbol  = char;
	            list_of_idx = [ char, char, ... ]"""
	def __init__(self, mat_symbol, list_of_idx):
		self.mat_symbol = mat_symbol
		if len(list_of_idx) % 2 == 0:
			self.list_of_idx = copy.deepcopy(list_of_idx)
		else:
			raise TypeError
	def update_index(self, original_idx, new_idx):
		# Updating Only One Index
		if original_idx in self.list_of_idx:
			for i in range(len(self.list_of_idx)):
				if self.list_of_idx[i] == original_idx:
					self.list_of_idx[i] = new_idx
	def get_symbol(self):
		return self.mat_symbol
	def get_index(self):
		return copy.deepcopy(self.list_of_idx)
	# quick info
	def return_info(self):
		return " " + self.mat_symbol + "_{ " + ", ".join(self.list_of_idx) + " } "
	def print_info(self):
		print( self.return_info() )


class vector_state(object):
	"""vector_state is a CISD... or CCSD... state"""
	def __init__(self, list_of_occ_idx=None, list_of_vrt_idx=None):
		self.list_of_occ_idx = copy.deepcopy(list_of_occ_idx)
		self.list_of_vrt_idx = copy.deepcopy(list_of_vrt_idx)
	def get_occ_idx(self):
		return copy.deepcopy(self.list_of_occ_idx)
	def get_vrt_idx(self):
		return copy.deepcopy(self.list_of_vrt_idx)
	def update_occ_idx(self, occ_idx):
		self.list_of_occ_idx = copy.deepcopy(occ_idx)
	def update_vrt_idx(self, vrt_idx):
		self.list_of_vrt_idx = copy.deepcopy(vrt_idx)
	def update_index(self, old_index, new_index):
		if self.list_of_occ_idx == None:
			pass
		else:
			if old_index in self.list_of_occ_idx:
				for i in range(len(self.list_of_occ_idx)):
					if self.list_of_occ_idx[i] == old_index:
						self.list_of_occ_idx[i] = new_index
			if old_index in self.list_of_vrt_idx:
				for i in range(len(self.list_of_vrt_idx)):
					if self.list_of_vrt_idx[i] == old_index:
						self.list_of_vrt_idx[i] = new_index
	def return_info(self):
		if self.list_of_occ_idx != None and self.list_of_vrt_idx != None:
			output = " | \Psi_{ "
			for item in self.list_of_occ_idx:
				output += item + ' '
			output += '}^{ '
			for item in self.list_of_vrt_idx:
				output += item + ' '
			output += '} >'
			return output
		else:
			return " |0> "
	def print_info(self):
		print(self.return_info())


class term(object):
	"""A class holding summation, delta and operator class objects; also has a coeff stored."""
	def __init__(self, list_of_sum, list_of_delta, op_string, list_of_mat_term, vector_state_term, coeff=1.0):
		self.list_of_sum       = copy.deepcopy(list_of_sum)
		self.list_of_delta     = copy.deepcopy(list_of_delta)
		self.op_string         = copy.deepcopy(op_string)
		self.coeff             = coeff
		self.vector_state_term = copy.deepcopy(vector_state_term)
		self.list_of_mat_term  = copy.deepcopy(list_of_mat_term)
	# Readers:
	def get_sum_objs(self):
		return copy.deepcopy(self.list_of_sum)
	def get_delta(self):
		return copy.deepcopy(self.list_of_delta)
	def get_op_string(self):
		return copy.deepcopy(self.op_string)
	def get_coeff(self):
		return copy.deepcopy(self.coeff)
	def get_vec_state(self):
		return copy.deepcopy(self.vector_state_term)
	def get_mat_term(self):
		return copy.deepcopy(self.list_of_mat_term)
	# Writers:
	def update_coeff(self, new_coeff):
		self.coeff = new_coeff
	def update_delta(self, new_delta_list):
		self.list_of_delta = copy.deepcopy(new_delta_list)
	def update_sum(self, new_sum_list):
		self.list_of_sum = copy.deepcopy(new_sum_list)
	def update_op_string(self, new_op_string):
		self.op_string = copy.deepcopy(new_op_string)
	def update_vec_state(self, vec_state):
		self.vector_state_term = copy.deepcopy(vec_state)
	def update_mat_term(self, new_mat_terms):
		self.list_of_mat_term = copy.deepcopy(new_mat_terms)
	# quick info (for debugging)
	def return_info(self):
		output_buf = str(self.coeff) + ' '
		if self.list_of_sum != None:
			for item in self.list_of_sum:
				output_buf += item.return_info() 
		for item in self.list_of_mat_term:
			output_buf += item.return_info()
		for item in self.list_of_delta:
			output_buf += item.return_info()
		for item in self.op_string:
			output_buf += item.return_info()
		if self.vector_state_term != None:
			output_buf += self.vector_state_term.return_info()
		return output_buf
	def print_info(self):
		print( self.return_info() )



def swap_op(op1, op2):
	# Rules:
	# a+_p a+_q  = - a+_q a+_p
	# a_p  a_q   = - a_q  a_p
	#
	# a+_p a_q  = \delta_{pq} - a_q  a+_p
	# a_q  a+_p = \delta_{qp} - a+_p a_q
	#
	if op1.name() == op2.name():
		# No Delta Function
		# Send back an empty list to hold the return orders
		empty_list = []
		return [-1.0], \
		       empty_list, \
		       [ copy.deepcopy(op2), copy.deepcopy(op1) ]
	else:
		# Add Delta Function
		return [1.0, -1.0], \
		       [ delta_func( op1.op_index(), op2.op_index(), op1.space(), op2.space() ) ], \
		       [ copy.deepcopy(op2), copy.deepcopy(op1) ] 



def swap_position(left_idx, term_obj):
	# Swap primitive operators, leave out delta.
	#
	op_string = term_obj.get_op_string()
	prefix    = copy.deepcopy(op_string[:left_idx])
	suffix    = copy.deepcopy(op_string[left_idx+2:])
	coeff, delta_func, swapped_op = swap_op(op_string[left_idx], op_string[left_idx+1])	
	#
	#
	if len(coeff) == 1:
		# No delta function
		new_op_string = [ copy.deepcopy( prefix + swapped_op + suffix ) ]
		new_term_obj  = copy.deepcopy(term_obj)
		new_term_obj.update_coeff( coeff[0] * term_obj.get_coeff() )
		new_term_obj.update_op_string( new_op_string[0] )
		new_terms = [new_term_obj]
	else:
		# There's a delta function
		new_op_string = [ copy.deepcopy( prefix + suffix ), copy.deepcopy( prefix + swapped_op + suffix ) ]
		new_terms = []
		#
		# First term contains delta
		new_term_obj1 = copy.deepcopy(term_obj)
		new_term_obj1.update_coeff( coeff[0] * term_obj.get_coeff() )
		new_term_obj1.update_op_string( new_op_string[0] )
		new_term_obj1.update_delta( term_obj.get_delta() + delta_func )
		#
		new_terms += [new_term_obj1]
		#
		# Second term has no delta
		new_term_obj2 = copy.deepcopy(term_obj)
		new_term_obj2.update_coeff( coeff[1] * term_obj.get_coeff() )
		new_term_obj2.update_op_string( new_op_string[1] )
		#
		new_terms += [new_term_obj2]
	#
	#
	return copy.deepcopy(new_terms)


def operator_wanted(current_operator):
	# Takes prmitive operator and evaluates it.
	if current_operator.name() == 'create' and current_operator.space() == 'occupied':
		return True
	elif current_operator.name() == 'annihilate' and current_operator.space() == 'virtual':
		return True
	else:
		return False



def term_is_zero(op_string):
	# Check the last operator in the string
	# if a |R> or m+ |R> found, send back message that this term is zero
	last_op = op_string[-1]
	if operator_wanted(last_op):
		return True
	else:
		return False



def find_swap_position(op_string):
	# return the left index for calling swap_position,
	# if index == -1, then move on to its left and look for another one,
	# then if nothing found, figure out a way to send a "Ordered" Message back.
	length = len(op_string)
	i = length - 1
	op_found = False
	while i >= 0 and not op_found:
		if operator_wanted( op_string[i] ):
			op_found = True
		else:
			i -= 1
	# if op_wanted not found in this string,
	# the return index = -1, which indicates op_wanted not in here,
	# go ahead and add this op_string to the final list.
	return i

def find_term_swap_position(term_obj):
	op_string = term_obj.get_op_string()
	return find_swap_position(op_string)

# def update_coeffs(current_coeff, returned_coeffs):
# 	coeffs = copy.deepcopy(returned_coeffs)
# 	# print("CHECK HERE COEFF =", returned_coeffs)
# 	for i in range(len(coeffs)):
# 		coeffs[i] = coeffs[i] * current_coeff
# 	return coeffs


# def print_op_string(op_string):
# 	output_buf = ""
# 	for item in op_string:
# 		output_buf += item.return_info()
# 	print(output_buf)



class state(object):
	"""state object contains term class and is running all the re-ordering processes"""
	def __init__(self, list_of_term_obj, state_collapsed=True):
		self.list_of_term_obj = copy.deepcopy(list_of_term_obj)
		self.state_collapsed  = state_collapsed
		# self.op_string = op_string
		# self.state = 'reference'
		
	def get_list_of_terms(self):
		return copy.deepcopy(self.list_of_term_obj)


	def reorder(self):
		# Move Annihilate Virtual to the right side
		#
		# print("Re-ordering This Term:")
		# for item in self.list_of_term_obj:
		# 	item.print_info()

		final_list   = []
		# final_coeffs = []
		ordered = False
		current_list   = copy.deepcopy(self.list_of_term_obj) 
		# current_coeffs = [ self.coeff ]
		
		while not ordered:
			# print("Ordering...")
			new_list   = []
			# new_coeffs = []
			for i in range(len(current_list)):
				op_string = current_list[i].get_op_string()
				op_wanted_idx = find_swap_position(op_string)
				# print("op_wanted_idx =", op_wanted_idx) 
				#
				#
				if op_wanted_idx == -1:
					# This term is non-zero, send it to the final_list
					final_list   += [ copy.deepcopy(current_list[i]) ]
					# print("Non-Zero Term Found: ", current_list[i].return_info())
				else:
			 		# Re-Order is Needed
			 		if term_is_zero(op_string):
			 			# Do nothing
			 			# print("Removing Zero Term: ", current_list[i].return_info())
			 			pass
			 		else:
			 			# SWAP HERE
			 			new_terms   = swap_position(op_wanted_idx, current_list[i])
			 			new_list   += new_terms

			current_list = copy.deepcopy(new_list)
			# print("==================================")
			# print("CURRENT_LIST")
			# for item in current_list:
			# 	item.print_info()
			# print("==================================")
			# print("FINAL_LIST")
			# for item in final_list:
			# 	item.print_info()
			# print("==================================")

			if len(current_list) == 0:
				# All term are ordered and zero terms are removed,
				# when the new_list is empty
				ordered = True
		# print("Ordered Operator Strings:")
		# for item in final_list:
		# 	item.print_info()
		return state( final_list, False )


	# def compare_sum_range(self, sum_obj1, sum_obj2):
	# 	if sum_obj1.get_index() == sum_obj2.get_index() and \
	# 	   sum_obj1.get_field() == sum_obj2.get_field():
	# 		# do ...
	# 		pass
	# 	else:
	# 		pass

	def take_sum_subset(self, sums_to_condense):
		if len(sums_to_condense) > 2:
			print("INTERNAL DESIGN ERROR, HAVE TO RETHINK ABOUT THIS! KILLED IN take_sum_subset FUNCTION")
			raise RuntimeError
		lower_limit = sums_to_condense[0].get_lower_limit() or sums_to_condense[1].get_lower_limit()
		# upper_limit = sums_to_condense[0].get_upper_limit() or sums_to_condense[1].get_upper_limit()
		sum_idx   = sums_to_condense[0].get_index()
		idx_field = sums_to_condense[0].get_field()
		return summation(sum_idx, idx_field, lower_limit)

	def remove_term_redundant_sum(self, list_of_sum):
		# First find all existing indices and make them a list
		# Then retrieve all corresponding sums having the same index,
		# send them to 'take_sum_subset' to combine the summation conditions.
		sum_index = []
		for each_sum in list_of_sum:
			if each_sum.get_index() not in sum_index:
				sum_index += [ each_sum.get_index() ]
		# print("THIS SUM INDICES =",sum_index)
		new_sum_list = []
		for each_index in sum_index:
			sums_to_condense = [ each_sum for each_sum in list_of_sum if each_sum.get_index() == each_index ]
			# print(sums_to_condense)
			if len(sums_to_condense) > 1:
				# print("Condensing...")
				new_sum = self.take_sum_subset(sums_to_condense)
			else:
				# print("Only One Sum is Present...")
				new_sum = copy.deepcopy(sums_to_condense[0])
			new_sum_list += [ new_sum ]
		return new_sum_list

	def remove_redundant_sum(self):
		# ONLY USE THIS IMMEDIATELY AFTER REMOVING DELTA FUNCTIONS
		# 
		for each_term in self.list_of_term_obj:
			sums = each_term.get_sum_objs()
			if len(sums) >0:
				new_sums = self.remove_term_redundant_sum(sums)
				each_term.update_sum(new_sums)

	def condense_delta(self, term_obj):
		deltas = term_obj.get_delta()
		new_index = max( deltas[0].indices() )
		old_index = min( deltas[0].indices() )
		# The alphabetically lower char is chosen
		# 1) summations and limits
		sums = term_obj.get_sum_objs()
		new_sums = []
		for each_term in sums:
			each_term.update_index(old_index,new_index)
			each_term.update_limits(old_index,new_index)
			new_sums += [ copy.deepcopy(each_term) ]
		# 2) Operators
		op_string = term_obj.get_op_string()
		new_op_string = []
		for each_term in op_string:
			each_term.update_index(old_index,new_index)
			new_op_string += [ copy.deepcopy(each_term) ]
		# print("NEW OP STRING", new_op_string)
		# 3) other delta functions
		new_deltas = []
		for each_term in deltas[1:]:
			each_term.update_index(old_index,new_index)
			new_deltas += [copy.deepcopy(each_term)]
		# 4) Matrix (like) terms
		mat_terms = term_obj.get_mat_term()
		new_mat_terms = []
		for each_term in mat_terms:
			each_term.update_index(old_index, new_index)
			new_mat_terms += [ copy.deepcopy(each_term) ]
		# 5) coefficient is not changing for the whole term
		#
		# Update them into a new term to iterate over.
		new_term = copy.deepcopy(term_obj)
		new_term.update_sum(new_sums)
		new_term.update_delta(new_deltas)
		new_term.update_op_string(new_op_string)
		new_term.update_mat_term(new_mat_terms)
		return new_term

	def remove_delta(self):
		new_list =  []
		# Loop over all the terms in this state.
		for item in self.list_of_term_obj:
			if len(item.get_delta()) == 0:
				# No delta in this term
				new_list += [copy.deepcopy(item)]
			else:
				new_term = self.condense_delta(item)
				while len(new_term.get_delta()) > 0:
					new_term = self.condense_delta(new_term)
				new_list += [copy.deepcopy(new_term)]
		return state(new_list, self.state_collapsed)


	def sort_alphabetic(self, op_string):
		new_coeff = 1.0
		temp_string = copy.deepcopy(op_string)
		size = len(temp_string)
		for i in range(size-1):  # sort thru n-1 times
			for j in range(size-i-1): # bubble sort, each time loop over n-i-1 pairs
				if temp_string[j].op_index() > temp_string[j+1].op_index():
					temp_string = copy.deepcopy(temp_string[:j])   + [copy.deepcopy(temp_string[j+1])] +\
					              [copy.deepcopy(temp_string[j])]  + copy.deepcopy(temp_string[j+2:])
					new_coeff *= -1.0
		return temp_string, new_coeff


	def sort_term_op_string(self, term_obj):
		op_string = term_obj.get_op_string()
		new_coeff = 1.0
		size = len(op_string) // 2  # half the size is needed since No. of a+ == No. of a
		# 1) Separation of creation/ annihilation operators
		check_idx = 0
		temp_string = copy.deepcopy(op_string)
		while check_idx < size :
			# Create on the left, annihilate on the right
			# There is ONLY a+ on virtual and a on occupied left! Other zero terms are all removed.
			# Check a from position zero to size-1, move a to the right-most position if found.
			if temp_string[check_idx].name() == 'annihilate':
				temp_string = copy.deepcopy( temp_string[:check_idx] ) + copy.deepcopy(temp_string[check_idx+1:]) +\
				              [copy.deepcopy( temp_string[check_idx] )]
				new_coeff  *= pow(-1.0, 2*size - 1 - check_idx) 
			else:
				check_idx += 1

		# 2) Sorting alphabetically within each group
		#
		# a) Sort creations on virtual
		creation_string = copy.deepcopy( temp_string[:size] )
		creation_string, temp_coeff = self.sort_alphabetic( creation_string )
		new_coeff  *= temp_coeff



		# b) Sort annihilations on occupied
		annihilation_string = copy.deepcopy( temp_string[size:] )
		annihilation_string, temp_coeff = self.sort_alphabetic( annihilation_string )
		new_coeff  *= temp_coeff

		new_op_string = copy.deepcopy(creation_string + annihilation_string)
		new_term_obj  = copy.deepcopy(term_obj)
		new_term_obj.update_op_string(new_op_string)
		new_term_obj.update_coeff( new_term_obj.get_coeff() * new_coeff )
		return new_term_obj


	def sort_state_op_string(self):
		# In here, triples or higher are not removed, but they will be in CISD/CCSD calculations.
		# Sort all the operators using sort_op_string function
		#
		new_terms = []
		for item in self.list_of_term_obj:
			new_terms += [ self.sort_term_op_string(item) ]
		
		return state(new_terms, self.state_collapsed)


	def collapse_operators(self):
		# This function does this: e+ f+ m n |0> ==> |\Psi_{mn}^{ef} >
		if self.state_collapsed:
			new_state = state( self.list_of_term_obj )
		else:
			# ... collapse it into a vector_state object
			new_state = self.sort_state_op_string()
			# new_state.print_info()
			for each_term in new_state.list_of_term_obj:
				if len(each_term.get_op_string()) > 1:
					list_of_index = []
					for op in each_term.get_op_string():
						list_of_index += [ op.op_index() ]
					each_term.update_op_string([])
					each_term.update_vec_state( vector_state( list_of_index[len(list_of_index)//2:], list_of_index[:len(list_of_index)//2] ) )
			new_state.state_collapsed = True
		return new_state

	def populate_operators(self):
		if self.state_collapsed:
			# ... recover string of operators from the collpased matrix_term obj
			for each_term in self.list_of_term_obj:
				vec_state = each_term.get_vec_state()
				list_of_occ_idx = vec_state.get_occ_idx()
				list_of_vrt_idx = vec_state.get_vrt_idx()
				if list_of_occ_idx != None and list_of_vrt_idx != None:
					new_op_string  = [ operator('create','virtual',      item) for item in list_of_vrt_idx ]
					new_op_string += [ operator('annihilate','occupied', item) for item in list_of_occ_idx ]
				else:
					new_op_string = []
				each_term.update_op_string( copy.deepcopy(each_term.get_op_string()) + new_op_string )
				each_term.update_vec_state( vector_state())
			self.state_collpased = False
		else:
			pass


	def remove_higher_order_term(self):
		# At this final stage, any Triples or higher will be removed
		# print("=============== Removing Higher-Order Terms ====================")
		new_list_of_terms = []
		if self.state_collapsed:
			# Length of list_of_vrt_idx must be less than or equal to 2 
			# (same to the list_of_occ_idx but only one is needed to be checked)
			# discard anything else.
			for each_term in self.list_of_term_obj:
				vec_state = each_term.get_vec_state()
				list_of_vrt_idx = vec_state.get_vrt_idx()
				if list_of_vrt_idx != None:
					if len(list_of_vrt_idx) <= 2:
						new_list_of_terms += [copy.deepcopy(each_term)]
						# print("Found ref/single/double Term:", each_term.return_info())
					else:
						# print("Removing   high-order   Term:", each_term.return_info())
						pass
				else:
					# print("Found ref/single/double Term:", each_term.return_info())
					new_list_of_terms += [copy.deepcopy(each_term)]
		else:
			# length of op_string must be less than or equal to 4
			# discard anything else.
			for each_term in self.list_of_term_obj:
				op_string = each_term.get_op_string()
				if len(op_string) <= 4:
					new_list_of_terms += [copy.deepcopy(each_term)]
					# print("Found ref/single/double Term:", each_term.return_info())
				else:
					# print("Removing   high-order   Term:", each_term.return_info())
					pass
		#
		#
		return state( new_list_of_terms, self.state_collapsed )

	def split_occ_sum_index(self, current_term):
		term_obj = copy.deepcopy(current_term)
		list_of_sum = term_obj.get_sum_objs()
		vec_state   = term_obj.get_vec_state()
		# 1) check occupied indices
		list_of_occ_idx = vec_state.get_occ_idx()
		#
		for i in range(len(list_of_sum)):
			if list_of_sum[i].get_index() == list_of_occ_idx[0]:
				sum_obj1 = list_of_sum[i] # A shallow copy
			elif list_of_sum[i].get_index() == list_of_occ_idx[1]:
				sum_obj2 = list_of_sum[i]

		# None set, a fresh new one.
		if sum_obj1.get_lower_limit() == None and sum_obj2.get_lower_limit() == None:
			sum_obj2.update_lower_limit( sum_obj1.get_index() )
			term_obj.update_sum(list_of_sum)
			term1 = copy.deepcopy(term_obj)
			# now swap
			sum_obj2.update_lower_limit( None )
			sum_obj1.update_lower_limit( sum_obj2.get_index() )
			term_obj.update_sum(list_of_sum)
			new_occ_idx = [ list_of_occ_idx[1], list_of_occ_idx[0] ]
			vec_state.update_occ_idx(new_occ_idx)
			term_obj.update_vec_state(vec_state)
			term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
			term2 = copy.deepcopy( term_obj )
			new_terms = [ term1, term2 ]

		elif sum_obj1.get_lower_limit() == None and sum_obj2.get_lower_limit() != None:
			# One side need to use the upper_limit
			if sum_obj2.get_lower_limit() == sum_obj1.get_index():
				new_terms = [copy.deepcopy(term_obj)]
			else:
				# set sum_obj1's upper_limit and lower_limit
				sum_obj1.update_upper_limit( sum_obj2.get_index() )
				term_obj.update_sum(list_of_sum)
				term1 = copy.deepcopy(term_obj)
				#
				sum_obj1.update_upper_limit( None )
				sum_obj1.update_lower_limit( sum_obj2.get_index() )
				term_obj.update_sum(list_of_sum)
				new_occ_idx = [ list_of_occ_idx[1], list_of_occ_idx[0] ]
				vec_state.update_occ_idx(new_occ_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				term2 = copy.deepcopy( term_obj )
				new_terms = [ term1, term2 ]
			
		elif sum_obj1.get_lower_limit() != None and sum_obj2.get_lower_limit() == None:
			if sum_obj1.get_lower_limit() == sum_obj2.get_index():
				# easy one term problem
				new_occ_idx = [ list_of_occ_idx[1], list_of_occ_idx[0] ]
				vec_state.update_occ_idx(new_occ_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				new_terms = [copy.deepcopy(term_obj)]
			else:
				# Two term problem, set sum_obj2's lower_limit and upper_limit
				sum_obj2.update_lower_limit( sum_obj1.get_index() )
				term_obj.update_sum(list_of_sum)
				term1 = copy.deepcopy(term_obj)
				# next term
				sum_obj2.update_lower_limit( None )
				sum_obj2.update_upper_limit( sum_obj1.get_index() )
				term_obj.update_sum(list_of_sum)
				new_occ_idx = [ list_of_occ_idx[1], list_of_occ_idx[0] ]
				vec_state.update_occ_idx(new_occ_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				term2 = copy.deepcopy( term_obj )
				new_terms = [ term1, term2 ]

		elif sum_obj1.get_lower_limit() != None and sum_obj2.get_lower_limit() != None: 
			if sum_obj2.get_lower_limit() == sum_obj1.get_index():
				new_terms = [copy.deepcopy(term_obj)]

			elif sum_obj1.get_lower_limit() == sum_obj2.get_index():
				new_occ_idx = [ list_of_occ_idx[1], list_of_occ_idx[0] ]
				vec_state.update_occ_idx(new_occ_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				new_terms = [copy.deepcopy(term_obj)]
			else:
				sum_obj1.update_upper_limit( sum_obj2.get_index() )
				term_obj.update_sum(list_of_sum)
				term1 = copy.deepcopy(term_obj)
				sum_obj1.update_upper_limit( None )
				sum_obj2.update_upper_limit( sum_obj1.get_index() )
				term_obj.update_sum(list_of_sum)
				new_occ_idx = [ list_of_occ_idx[1], list_of_occ_idx[0] ]
				vec_state.update_occ_idx(new_occ_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				term2 = copy.deepcopy( term_obj )
				new_terms = [term1, term2]
		else:
			print("THIS IS IMPOSSIBLE TO BE RAN.")
			raise AssertionError
		return new_terms
				

	def split_vrt_sum_index(self, current_term):
		term_obj = copy.deepcopy(current_term)
		list_of_sum = term_obj.get_sum_objs()
		vec_state   = term_obj.get_vec_state()
		# check virtual indices
		list_of_vrt_idx = vec_state.get_vrt_idx()
		#
		for i in range(len(list_of_sum)):
			if list_of_sum[i].get_index() == list_of_vrt_idx[0]:				
				sum_obj1 = list_of_sum[i] # A shallow copy
			elif list_of_sum[i].get_index() == list_of_vrt_idx[1]:
				sum_obj2 = list_of_sum[i]
		
		# Lastest way to doing it...
		if sum_obj1.get_lower_limit() == None and sum_obj2.get_lower_limit() == None:
			sum_obj2.update_lower_limit( sum_obj1.get_index() )
			term_obj.update_sum(list_of_sum)
			term1 = copy.deepcopy(term_obj)
			# now swap
			sum_obj2.update_lower_limit( None )
			sum_obj1.update_lower_limit( sum_obj2.get_index() )
			term_obj.update_sum(list_of_sum)
			new_vrt_idx = [ list_of_vrt_idx[1], list_of_vrt_idx[0] ]
			vec_state.update_vrt_idx(new_vrt_idx)
			term_obj.update_vec_state(vec_state)
			term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
			term2 = copy.deepcopy( term_obj )
			new_terms = [ term1, term2 ]

		elif sum_obj1.get_lower_limit() == None and sum_obj2.get_lower_limit() != None:
			# One side need to use the upper_limit
			if sum_obj2.get_lower_limit() == sum_obj1.get_index():
				new_terms = [copy.deepcopy(term_obj)]
			else:
				# set sum_obj1's upper_limit and lower_limit
				sum_obj1.update_upper_limit( sum_obj2.get_index() )
				term_obj.update_sum(list_of_sum)
				term1 = copy.deepcopy(term_obj)
				#
				sum_obj1.update_upper_limit( None )
				sum_obj1.update_lower_limit( sum_obj2.get_index() )
				term_obj.update_sum(list_of_sum)
				new_vrt_idx = [ list_of_vrt_idx[1], list_of_vrt_idx[0] ]
				vec_state.update_vrt_idx(new_vrt_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				term2 = copy.deepcopy( term_obj )
				new_terms = [ term1, term2 ]
			
		elif sum_obj1.get_lower_limit() != None and sum_obj2.get_lower_limit() == None:
			if sum_obj1.get_lower_limit() == sum_obj2.get_index():
				# easy one term problem
				new_vrt_idx = [ list_of_vrt_idx[1], list_of_vrt_idx[0] ]
				vec_state.update_vrt_idx(new_vrt_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				new_terms = [copy.deepcopy(term_obj)]
			else:
				# Two term problem, set sum_obj2's lower_limit and upper_limit
				sum_obj2.update_lower_limit( sum_obj1.get_index() )
				term_obj.update_sum(list_of_sum)
				term1 = copy.deepcopy(term_obj)
				# next term
				sum_obj2.update_lower_limit( None )
				sum_obj2.update_upper_limit( sum_obj1.get_index() )
				term_obj.update_sum(list_of_sum)
				new_vrt_idx = [ list_of_vrt_idx[1], list_of_vrt_idx[0] ]
				vec_state.update_vrt_idx(new_vrt_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				term2 = copy.deepcopy( term_obj )
				new_terms = [ term1, term2 ]

		elif sum_obj1.get_lower_limit() != None and sum_obj2.get_lower_limit() != None: 
			if sum_obj2.get_lower_limit() == sum_obj1.get_index():
				new_terms = [copy.deepcopy(term_obj)]

			elif sum_obj1.get_lower_limit() == sum_obj2.get_index():
				new_vrt_idx = [ list_of_vrt_idx[1], list_of_vrt_idx[0] ]
				vec_state.update_vrt_idx(new_vrt_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				new_terms = [copy.deepcopy(term_obj)]
			else:
				sum_obj1.update_upper_limit( sum_obj2.get_index() )
				term_obj.update_sum(list_of_sum)
				term1 = copy.deepcopy(term_obj)
				sum_obj1.update_upper_limit( None )
				sum_obj2.update_upper_limit( sum_obj1.get_index() )
				term_obj.update_sum(list_of_sum)
				new_vrt_idx = [ list_of_vrt_idx[1], list_of_vrt_idx[0] ]
				vec_state.update_vrt_idx(new_vrt_idx)
				term_obj.update_vec_state(vec_state)
				term_obj.update_coeff( term_obj.get_coeff() * -1.0 )
				term2 = copy.deepcopy( term_obj )
				new_terms = [term1, term2]
		else:
			print("THIS IS IMPOSSIBLE TO BE RAN.")
			raise AssertionError
		return new_terms



	def split_doubles_sum_index(self):
		new_list_of_terms = []
		if self.state_collapsed:  # A collapsed state is REQUIRED.
			for each_term in self.list_of_term_obj:
				# This is the major checkpoint
				vec_state   = each_term.get_vec_state()
				if vec_state.get_vrt_idx() == None:
					new_list_of_terms += [copy.deepcopy(each_term)]
				elif len(vec_state.get_vrt_idx()) == 2:
					# check both index if they already have got lower_limit
					# The one has the lower_limit (resulting from removing the delta) should be swapped 
					# to the first position (where lower_limit is not worried about)
					# since the separation of summation range will create new lower_limit for the second index,
					# The second index must be free of lower_limit by the time of separation, thus a swap is needed.
					#
					new_terms = self.split_occ_sum_index(each_term)
					if len(new_terms) == 1:
						final_terms = self.split_vrt_sum_index(new_terms[0])
					elif len(new_terms) == 2:
						final_terms  = self.split_vrt_sum_index( new_terms[0] )
						final_terms += self.split_vrt_sum_index( new_terms[1] )
					else:
						print("STUPID SPLIT FUNCTION ERROR, MUST RETURN ONE OR TWO TERMS IN A LIST")
						raise RuntimeError
					#
					#
					new_list_of_terms += copy.deepcopy( final_terms ) 
				else:
					# This term's summations don't need to be handled.
					new_list_of_terms += [copy.deepcopy(each_term)]
		else:
			print("PLEASE AVOID GETTING SUMMATION SEPARATION WHILE IT'S STILL NOT FULLY COMPUTED!")
			raise RuntimeError
		#
		return state( new_list_of_terms )

	




	#
	# Combine Terms Helper Functions
	#
	def term_swap_index(self, term_obj, index1, index2):
		# This function only swaps summatiom objs, matrix_like term objs and vector_state objs
		# delta, op_string and coeff are not swapped due to this feature is only needed at the last step.
		new_term = copy.deepcopy(term_obj)
		# 1) Replace all index1 with 'TEMP'
		# 2) Replace all index2 with index1
		# 3) Replace all 'TEMP' with index2
		temp_index = 'TEMP'
		# step one
		sums = new_term.get_sum_objs()
		for each_sum in sums:
			each_sum.update_index(index1, temp_index)
			each_sum.update_limits(index1, temp_index)
		new_term.update_sum(sums)
		#
		mat_terms = new_term.get_mat_term()
		for each_term in mat_terms:
			each_term.update_index(index1, temp_index)
		new_term.update_mat_term(mat_terms)
		#
		vec_state = new_term.get_vec_state()
		vec_state.update_index(index1, temp_index)
		new_term.update_vec_state(vec_state)
		#
		# step two
		sums = new_term.get_sum_objs()
		for each_sum in sums:
			each_sum.update_index(index2, index1)
			each_sum.update_limits(index2, index1)
		new_term.update_sum(sums)
		#
		mat_terms = new_term.get_mat_term()
		for each_term in mat_terms:
			each_term.update_index(index2, index1)
		new_term.update_mat_term(mat_terms)
		#
		vec_state = new_term.get_vec_state()
		vec_state.update_index(index2, index1)
		new_term.update_vec_state(vec_state)
		#
		# step three
		sums = new_term.get_sum_objs()
		for each_sum in sums:
			each_sum.update_index(temp_index, index2)
			each_sum.update_limits(temp_index, index2)
		new_term.update_sum(sums)
		#
		mat_terms = new_term.get_mat_term()
		for each_term in mat_terms:
			each_term.update_index(temp_index, index2)
		new_term.update_mat_term(mat_terms)
		#
		vec_state = new_term.get_vec_state()
		vec_state.update_index(temp_index, index2)
		new_term.update_vec_state(vec_state)
		#
		return new_term


	def term_occ_to_target(self, current_term, target_occ_idx):
		term_obj = copy.deepcopy(current_term)
		term_vec = term_obj.get_vec_state()
		vec_occ  = term_vec.get_occ_idx()
		if len(target_occ_idx) == 1:
			# Easy Stuff: Single Exicitation State
			if vec_occ[0] == target_occ_idx[0]:
				# Nothing needs to be changed.
				pass
			else:
				# Swap the index in vec_occ to target_occ
				term_obj = self.term_swap_index(term_obj, vec_occ[0], target_occ_idx[0])
				# Done... b/c the rest of items are emtpy by now

		elif len(target_occ_idx) == 2:
			# Double Excitation vector_state
			if target_occ_idx[0] in vec_occ:
				target_first_in_current = True
			else:
				target_first_in_current = False
			if target_occ_idx[1] in vec_occ:
				target_second_in_current = True
			else:
				target_second_in_current = False
			# print("Target First  =", target_first_in_current, target_occ_idx, vec_occ)
			# print("Target Second =", target_second_in_current)
			if target_first_in_current or target_second_in_current:
				# This must be true, or it's an error
				if target_first_in_current and target_second_in_current:
					# check if positions are the same
					if target_occ_idx[0] == vec_occ[0]:
						# Same Vector_state, do nothing here
						pass
					else:
						# A SWAP is carried out.
						term_obj = self.term_swap_index(term_obj, vec_occ[0], vec_occ[1])

				else:
					# Only One of the two is the same
					# check where they are the same:
					#     same position => Swap the different one Only; 
					#     different Position => Swap the different one first, swap the indices within state_vector
					if target_first_in_current and not target_second_in_current:
						if target_occ_idx[0] == vec_occ[0]:
							# e.g., ['m','n'] vs ['m','i']
							# swap position two
							term_obj = self.term_swap_index(term_obj, vec_occ[1], target_occ_idx[1])
						else:
							# e.g., ['m','n'] vs ['i','m']
							# swap two times
							term_obj = self.term_swap_index(term_obj, vec_occ[0], target_occ_idx[1])
							term_obj = self.term_swap_index(term_obj, target_occ_idx[0], target_occ_idx[1])
					elif not target_first_in_current and target_second_in_current:
						# target_second must be in current
						if target_occ_idx[1] == vec_occ[1]:
							# e.g., ['m', 'n'] vs ['i', 'n']
							# swap position one
							term_obj = self.term_swap_index(term_obj, vec_occ[0], target_occ_idx[0])
						else:
							# e.g., ['m','n'] vs ['n','i']
							# swap two times
							term_obj = self.term_swap_index(term_obj, vec_occ[1], target_occ_idx[0])
							term_obj = self.term_swap_index(term_obj, target_occ_idx[0], target_occ_idx[1])
					else:
						print("THROW AN ERROR HERE, THIS MAY NOT BE AN ERROR BUT NOT YET IMPLEMENTED.")
						raise RuntimeError
							
			else:
				# print("COMPLETELY DIFFERENT INDICES SHOW UP IN THE RESULTS, TERRIBLE ERROR")
				# raise RuntimeError
				# e.g. ['m','n'] (target) vs ['i','j'] (to replace)
				term_obj = self.term_swap_index(term_obj, target_occ_idx[0], vec_occ[0])
				term_obj = self.term_swap_index(term_obj, target_occ_idx[1], vec_occ[1])
		else:
			# It is impossbile to reach this branch, so throw an error
			print("IMPOSSIBLE RUNTIME ERROR, CHECK SWAP/REPLACE VECTOR_STATE INDEX FUNCTIONS")
			raise RuntimeError 
		return term_obj  # Already A Deep Copy


	def term_vrt_to_target(self, current_term, target_vrt_idx):
		term_obj = copy.deepcopy(current_term)
		term_vec = term_obj.get_vec_state()
		vec_vrt  = term_vec.get_vrt_idx()
		if len(target_vrt_idx) == 1:
			# Easy Stuff: Single Exicitation State
			if vec_vrt[0] == target_vrt_idx[0]:
				# Nothing needs to be changed.
				pass
			else:
				# Swap the index in vec_vrt to target_vrt
				term_obj = self.term_swap_index(term_obj, vec_vrt[0], target_vrt_idx[0])
				# Done... b/c the rest of items are emtpy by now

		elif len(target_vrt_idx) == 2:
			# Double Excitation vector_state
			if target_vrt_idx[0] in vec_vrt:
				target_first_in_current = True
			else:
				target_first_in_current = False
			if target_vrt_idx[1] in vec_vrt:
				target_second_in_current = True  
			else:
				target_second_in_current = False

			if target_first_in_current or target_second_in_current:
				# This must be true, or it's an error
				if target_first_in_current and target_second_in_current:
					# check if positions are the same
					if target_vrt_idx[0] == vec_vrt[0]:
						# Same Vector_state, do nothing here
						pass
					else:
						# A SWAP is carried out.
						term_obj = self.term_swap_index(term_obj, vec_vrt[0], vec_vrt[1])

				else:
					# Only One of the two is the same
					# check where they are the same:
					#     same position => Swap the different one Only; 
					#     different Position => Swap the different one first, swap the indices within state_vector
					if target_first_in_current and not target_second_in_current:
						if target_vrt_idx[0] == vec_vrt[0]:
							# e.g., ['e','f'] vs ['e','a']
							# swap position two
							term_obj = self.term_swap_index(term_obj, vec_vrt[1], target_vrt_idx[1])
						else:
							# e.g., ['e','f'] vs ['a','e']
							# swap two times
							term_obj = self.term_swap_index(term_obj, vec_vrt[0], target_vrt_idx[1])
							term_obj = self.term_swap_index(term_obj, target_vrt_idx[0], target_vrt_idx[1])
					elif not target_first_in_current and target_second_in_current:
						# target_second must be in current
						if target_vrt_idx[1] == vec_vrt[1]:
							# e.g., ['m', 'n'] vs ['i', 'n']
							# swap position one
							term_obj = self.term_swap_index(term_obj, vec_vrt[0], target_vrt_idx[0])
						else:
							# e.g., ['m','n'] vs ['n','i']
							# swap two times
							term_obj = self.term_swap_index(term_obj, vec_vrt[1], target_vrt_idx[0])
							term_obj = self.term_swap_index(term_obj, target_vrt_idx[0], target_vrt_idx[1])
					else:						
						print("THROW AN ERROR HERE, THIS MAY NOT BE AN ERROR BUT NOT YET IMPLEMENTED.")
						raise RuntimeError

							
			else:
				# print("COMPLETELY DIFFERENT INDICES SHOW UP IN THE RESULTS, TERRIBLE ERROR")
				# raise RuntimeError
				term_obj = self.term_swap_index(term_obj, vec_vrt[0], target_vrt_idx[0])
				term_obj = self.term_swap_index(term_obj, vec_vrt[1], target_vrt_idx[1])

		else:
			# It is impossbile to reach this branch, so throw an error
			print("IMPOSSIBLE RUNTIME ERROR, CHECK SWAP/REPLACE VECTOR_STATE INDEX FUNCTIONS")
			raise RuntimeError 
		return term_obj  # Already A Deep Copy



	def combine_terms(self):
		new_terms = []
		if self.state_collapsed:
			num_terms = len(self.list_of_term_obj)
			# Use the vector_state in the first term as the target for the rest of transformations.
			if num_terms > 0:
				target_state = self.list_of_term_obj[0].get_vec_state()
				if target_state.get_occ_idx() == None:				
					# No Need to Check Further, Reference State...
					new_terms = copy.deepcopy( self.list_of_term_obj )
				else:
					# Make sure the new states are not reference state
					target_occ_idx = target_state.get_occ_idx()   # e.g., ['m','n']
					target_vrt_idx = target_state.get_vrt_idx()   # e.g., ['e','f']
					#
					# Loop Over the Rest of the Terms to Transform them into this target_state vector
					new_terms += [ copy.deepcopy(self.list_of_term_obj[0]) ]
					for i in range(1, num_terms):
						current_term = copy.deepcopy(self.list_of_term_obj[i])
						current_term = self.term_occ_to_target(current_term, target_occ_idx)
						current_term = self.term_vrt_to_target(current_term, target_vrt_idx)
						new_terms += [copy.deepcopy(current_term)]

		else:
			print("PLEASE AVOID GETTING SUMMATION SEPARATION WHILE IT'S STILL NOT FULLY COMPUTED!")
			raise RuntimeError
		#
		return state( new_terms )


	def sort_list_of_sum(self, list_of_sum):
		sorted_list = copy.deepcopy(list_of_sum)
		num_sum = len(list_of_sum)
		if num_sum < 2:
			# No need to sort: empty or with One Element.
			pass
		else:
			for i in range(num_sum-1):
				for j in range(num_sum-i-1):
					if sorted_list[j].get_index() > sorted_list[j+1].get_index():
						temp_sum = copy.deepcopy(sorted_list[j+1]) 
						sorted_list[j+1] = copy.deepcopy( sorted_list[j] )
						sorted_list[j]   = copy.deepcopy( temp_sum ) 

		return sorted_list




	def sort_sum_by_index(self):
		# ONLY CALL THIS IN THE END (When Everything is Computed)
		new_terms = []
		if self.state_collapsed:
			for each_term in self.list_of_term_obj:
				current_term = copy.deepcopy(each_term)
				list_of_sum  = current_term.get_sum_objs()
				sorted_sum   = self.sort_list_of_sum(list_of_sum)
				current_term.update_sum(sorted_sum)
				new_terms   += [copy.deepcopy(current_term)]
		else:
			print("PLEASE AVOID SORTING SUMMATIONS WHILE IT'S STILL NOT FULLY COMPUTED!")
			raise RuntimeError
		return state( new_terms )



	# def condition_sum_in_term(self, current_term):
	# 	# NO SEPARATION OF INDEX IS NEEDED. JUST TAKE THE DOUBLE VECTORS AND TRY TO SET UP THE LOWER_LIMITS,
	# 	# IF ALREADY SET, SWAP TWO INDICES, THEN SET THE LOWER_LIMIT ON THE OTHER FREE INDEX
	# 	#
	# 	this_term = copy.deepcopy(current_term)
	# 	#
	# 	# 1) Fix Limits for Vector State Indices
	# 	vec_state = this_term.get_vec_state()
	# 	occ_idx   = vec_state.get_occ_idx()
	# 	vrt_idx   = vec_state.get_vrt_idx()

	# 	if occ_idx == None:
	# 		# Result is Reference state, do nothing.
	# 		pass
	# 	elif len(occ_idx) == 1:
	# 		# Result is Single State, no summation limits at all.
	# 		pass
	# 	elif len(occ_idx) == 2:
	# 		# Doubles, check if occ and vrt limits are enforced. 
	# 		sum_obj0 = [ copy.deepcopy(sum_obj) for sum_obj in this_term.get_sum_objs() if sum_obj.get_index() == occ_idx[0] ][0] 
	# 		sum_obj1 = [ copy.deepcopy(sum_obj) for sum_obj in this_term.get_sum_objs() if sum_obj.get_index() == occ_idx[1] ][0]
	# 		#
	# 		#
	# 		if sum_obj1.get_lower_limit() == occ_idx[0]:
	# 			# best case, pass
	# 			pass
	# 		elif sum_obj1.get_lower_limit() == None:
	# 			# Add a lower limit to it.
	# 			list_of_sum = this_term.get_sum_objs()
	# 			for each_sum in list_of_sum:
	# 				if each_sum.get_index() == occ_idx[1]:
	# 					each_sum.update_lower_limit(occ_idx[0])
	# 			this_term.update_sum(list_of_sum)
	# 		elif sum_obj1.get_lower_limit() != None:
	# 			# Worst case, swap vector state.
	# 			new_occ_idx = [ occ_idx[1], occ_idx[0] ]
	# 			new_vec_state = copy.deepcopy(vec_state)
	# 			new_vec_state.update_occ_idx( new_occ_idx )
	# 			this_term.update_vec_state(new_vec_state)
	# 			this_term.update_coeff( this_term.get_coeff() * -1.0 )

	# 			list_of_sum = this_term.get_sum_objs()
	# 			for each_sum in list_of_sum:
	# 				if each_sum.get_index() == new_occ_idx[1]:
	# 					if each_sum.get_lower_limit() == None:
	# 						each_sum.update_lower_limit(new_occ_idx[0])
	# 					else:
	# 						print("Fundamental Algrithm Flaw!")
	# 						raise NotImplementedError
	# 			this_term.update_sum(list_of_sum)
	# 		# Now occupied portion fixed.

	# 		sum_obj0 = [ copy.deepcopy(sum_obj) for sum_obj in this_term.get_sum_objs() if sum_obj.get_index() == vrt_idx[0] ][0] 
	# 		sum_obj1 = [ copy.deepcopy(sum_obj) for sum_obj in this_term.get_sum_objs() if sum_obj.get_index() == vrt_idx[1] ][0]
	# 		#
	# 		#
	# 		if sum_obj1.get_lower_limit() == vrt_idx[0]:
	# 			# best case, pass
	# 			pass
	# 		elif sum_obj1.get_lower_limit() == None:
	# 			# Add a lower limit to it.
	# 			list_of_sum = this_term.get_sum_objs()
	# 			for each_sum in list_of_sum:
	# 				if each_sum.get_index() == vrt_idx[1]:
	# 					each_sum.update_lower_limit(vrt_idx[0])
	# 			this_term.update_sum(list_of_sum)
	# 		elif sum_obj1.get_lower_limit() != None:
	# 			# Worst case, swap vector state.
	# 			new_vrt_idx = [ vrt_idx[1], vrt_idx[0] ]
	# 			new_vec_state = copy.deepcopy( this_term.get_vec_state() )
	# 			new_vec_state.update_vrt_idx( new_vrt_idx )
	# 			this_term.update_vec_state(new_vec_state)
	# 			this_term.update_coeff( this_term.get_coeff() * -1.0 )

	# 			list_of_sum = this_term.get_sum_objs()
	# 			for each_sum in list_of_sum:
	# 				if each_sum.get_index() == new_vrt_idx[1]:
	# 					if each_sum.get_lower_limit() == None:
	# 						each_sum.update_lower_limit(new_vrt_idx[0])
	# 					else:
	# 						print("Fundamental Algrithm Flaw!")
	# 						raise NotImplementedError
	# 			this_term.update_sum(list_of_sum)			
	# 	else:
	# 		print("TRIPLES ARE NOT STORED IN THIS SCOPE, ERROR!")
	# 		raise TypeError
	# 	#
	# 	#
	# 	return this_term



	# def condition_sum_limits(self):
	# 	# ONLY CALL THIS IN AFTER CALLING sort_sum_by_index()
	# 	# Will Erase Upper Limits of preceding sums
	# 	# Will Swap sums if limit index (which will be undefined) appears earlier than the loop index
	# 	new_terms = []
	# 	if self.state_collapsed:
	# 		for each_term in self.list_of_term_obj:
	# 			current_term     = copy.deepcopy(each_term)
	# 			conditioned_term = self.condition_sum_in_term(current_term)
	# 			new_terms   += [copy.deepcopy(conditioned_term)]
	# 	else:
	# 		print("PLEASE AVOID CONDITIONING SUMMATIONS WHILE IT'S STILL NOT FULLY COMPUTED!")
	# 		raise RuntimeError
	# 	return state( new_terms )			


	def sort_term_by_limits(self, current_term):
		# Sort the sum indices to make sure hihger index always on the right.
		this_term = copy.deepcopy(current_term)
		list_of_sum = this_term.get_sum_objs()
		num_sums = len(list_of_sum)
		if num_sums > 1:
			i = 0
			all_sum_idx = [ item.get_index() for item in list_of_sum ]
			while i < num_sums:
				# lower_limit appear on the left first than its actual loop
				lower_limit = list_of_sum[i].get_lower_limit()
				upper_limit = list_of_sum[i].get_upper_limit()
				if lower_limit in all_sum_idx[i+1:] or upper_limit in all_sum_idx[i+1:]:
					list_of_sum  = copy.deepcopy(list_of_sum[:i]) + copy.deepcopy(list_of_sum[i+1:]) + [copy.deepcopy(list_of_sum[i])]
					all_sum_idx = [ item.get_index() for item in list_of_sum ]
				else:
					i += 1
		#
		#
		this_term.update_sum(list_of_sum)
		return this_term

	def sort_index_by_limit(self):
		new_terms = []
		if self.state_collapsed:
			for each_term in self.list_of_term_obj:
				current_term = copy.deepcopy(each_term)
				sorted_term  = self.sort_term_by_limits(current_term)
				new_terms   += [copy.deepcopy(sorted_term)]
		else:
			print("PLEASE AVOID CONDITIONING SUMMATIONS WHILE IT'S STILL NOT FULLY COMPUTED!")
			raise RuntimeError
		return state( new_terms )	



	# def terms_are_opposite(self, term1, term2):
	# 	# This function now ONLY SUPPORT SAME INDEX AS EQUAL, MAY UPGRADE TO SAME FIELD AS EQUAL LATER, WE WILL SEE.
	# 	#
	# 	#  term( list_of_sum, list_of_delta, op_string, mat_term_obj, vec_term, coeff )
	# 	# 
	# 	# print(term1.get_coeff(), '+', term2.get_coeff(), '=',term1.get_coeff()+term2.get_coeff())
	# 	if term1.get_coeff() == -1 * term2.get_coeff():
	# 		# Continue checking
	# 		# (1) Check Summations:
	# 		term1_sums = term1.get_sum_objs()
	# 		term2_sums = term2.get_sum_objs()
	# 		if len(term1_sums) != len(term2_sums):
	# 			return False
	# 		else:
	# 			for i in range(len(term1_sums)):
	# 				if term1_sums[i].get_index() != term2_sums[i].get_index():
	# 					return False
	# 				# field is bounded to index for now, no need to check.
	# 				# if term1_sums[i].get_field() != term2_sums[i].get_field():
	# 					# return False
	# 				if term1_sums[i].get_lower_limit() != term2_sums[i].get_lower_limit():
	# 					return False
	# 				if term1_sums[i].get_upper_limit() != term2_sums[i].get_upper_limit():
	# 					return False
	# 		#
	# 		# (2) Check matrix terms
	# 		term1_mats = term1.get_mat_term()
	# 		term2_mats = term2.get_mat_term()
	# 		if len(term1_mats) != len(term2_mats):
	# 			return False
	# 		else:
	# 			for i in range(len(term1_mats)):
	# 				# term1_mats[i].print_info()
	# 				# term2_mats[i].print_info()
	# 				if term1_mats[i].get_symbol() != term2_mats[i].get_symbol():
	# 					return False
	# 				if term1_mats[i].get_index()  != term2_mats[i].get_index():
	# 					return False
	# 		#
	# 		# (3) Check vector states
	# 		vec_state1 = term1.get_vec_state()
	# 		vec_state2 = term2.get_vec_state()
	# 		if vec_state1.get_occ_idx() != vec_state2.get_occ_idx():
	# 			return False
	# 		if vec_state1.get_vrt_idx() != vec_state2.get_vrt_idx():
	# 			return False
	# 		print("Matching Found")
	# 		return True
	# 	else:
	# 		return False




	# def remove_opposite_sign_terms(self):
	# 	list_of_non_zero_terms = []
	# 	if self.state_collapsed:
	# 		# make sure it's in a collapsed term.
	# 		# use a while loop here, keep changing the list until no more term to remove.
	# 		current_terms = self.list_of_term_obj
	# 		collect_flags = [ True for i in range(len(current_terms)) ] # if this term will be kept 
	# 		for i in range(len(current_terms)):
	# 			if collect_flags[i] == True:
	# 				# Make sure this term has not been marked to be removed.
	# 				match_not_found = True
	# 				max_check_index = len(current_terms)
	# 				j = i + 1
	# 				while match_not_found and j < max_check_index:
	# 					if self.terms_are_opposite( current_terms[i], current_terms[j] ):
	# 						collect_flags[i] = False
	# 						collect_flags[j] = False
	# 						match_not_found  = False
	# 					else:
	# 						j += 1
	# 		list_of_non_zero_terms = [ copy.deepcopy(current_terms[i]) for i in range(len(current_terms)) if collect_flags[i] == True ]
	# 	else:
	# 		print("PLEASE AVOID GETTING CONTRACTION RELATION WHILE IT'S STILL NOT FULLY COMPUTED!")
	# 		raise RuntimeError
	# 	return state( list_of_non_zero_terms )


	def get_contraction_relation(self):
		list_of_contraction_terms = []
		if self.state_collapsed:
			# make sure it's in a collapsed term.
			for each_term in self.list_of_term_obj:
				vec_state = each_term.get_vec_state()
				list_of_occ_idx = vec_state.get_occ_idx()
				list_of_vrt_idx = vec_state.get_vrt_idx()
				if list_of_occ_idx != None and list_of_vrt_idx != None: # Singles and Doubles
					unwanted_sum_idx = list_of_occ_idx + list_of_vrt_idx
				else:  # This is for the reference state
					unwanted_sum_idx = []
				# print("Unwanted Index =", unwanted_sum_idx)
				contraction_term = copy.deepcopy(each_term)
				contraction_sums = [ item for item in contraction_term.get_sum_objs() if item.get_index() not in unwanted_sum_idx ]
				# print("surviving index =",contraction_sums)
				contraction_term.update_sum( contraction_sums )
				# Just in case, erase everything...
				contraction_term.update_delta([])
				contraction_term.update_vec_state(None)
				contraction_term.update_op_string([])
				list_of_contraction_terms += [contraction_term]
		else:
			print("PLEASE AVOID GETTING CONTRACTION RELATION WHILE IT'S STILL NOT FULLY COMPUTED!")
			raise RuntimeError
		return list_of_contraction_terms


	def get_return_vec_state(self):
		# vector_state contains summation indices inside
		list_of_vec_state = []
		if self.state_collapsed:
			for each_term in self.list_of_term_obj:
				list_of_vec_state += [ each_term.get_vec_state() ]
		else:
			print("PLEASE AVOID GETTING CONTRACTION RELATION WHILE IT'S STILL NOT FULLY COMPUTED!")
			raise RuntimeError
		return list_of_vec_state


	def return_info(self):
		output_buf = ""
		for item in self.list_of_term_obj:
			output_buf += item.return_info() + '\n'
		return output_buf

	def print_info(self):
		print("============= State Object =====================")
		for item in self.list_of_term_obj:
			item.print_info()
		print("================================================")








def concatenate(term_obj1, term_obj2):
	# def __init__(self, list_of_sum, list_of_delta, op_string, list_of_mat_term, vector_state_term, coeff=1.0):
	new_list_of_sum      =  term_obj1.get_sum_objs()   + term_obj2.get_sum_objs()
	new_list_of_delta    =  term_obj1.get_delta()      + term_obj2.get_delta()
	new_op_string        =  term_obj1.get_op_string()  + term_obj2.get_op_string()
	new_list_of_mat_term =  term_obj1.get_mat_term()   + term_obj2.get_mat_term()
	if term_obj1.get_vec_state() == None:
		new_vector_state_term = term_obj2.get_vec_state() # it's OK if obj2 has a None state
	elif term_obj2.get_vec_state() == None:
		new_vector_state_term = term_obj1.get_vec_state() # This is TBH really weird and should not happen
	else:
		print("ONLY ONE VECTOR STATE IS ALLOWED.")  # Definitely not good for now
		raise TypeError
	new_coeff = term_obj1.get_coeff() * term_obj2.get_coeff()
	return term(new_list_of_sum, new_list_of_delta, new_op_string, new_list_of_mat_term, new_vector_state_term, new_coeff)



def compute_state(primitive_state_obj):
	state_obj = copy.deepcopy(primitive_state_obj)  # Make a copy


	# print("Collapsed Form")
	# state_obj.print_info()
	

	state_obj.populate_operators()
	# print("Populated Form")
	# state_obj.print_info()


	state_obj = state_obj.reorder()
	# print("RE-ORDER Results")
	# state_obj.print_info()


	state_obj = state_obj.remove_delta()
	# print("Remove Delta")
	# state_obj.print_info()


	state_obj.remove_redundant_sum()
	# print("Clean Redundant Sum")
	# state_obj.print_info()


	state_obj = state_obj.sort_state_op_string()
	# print("Sorting Operators")
	# state_obj.print_info()


	state_obj = state_obj.collapse_operators()
	# print("Collapse Operators")
	# state_obj.print_info()

	state_obj = state_obj.remove_higher_order_term()

	# print("+++++++++++++++++++ DEBUG RESULTS ++++++++++++++++++++++")
	# state_obj.print_info()

	# print("+++++++++++++++++++ DEBUG RESULTS ++++++++++++++++++++++>>>>>>")
	# state_obj.print_info()
	state_obj = state_obj.split_doubles_sum_index()
	# state_obj.print_info()
	# print("+++++++++++++++++++ DEBUG RESULTS ++++++++++++++++++++++")

	# state_obj = state_obj.combine_terms()
	# state_obj.print_info()

	# state_obj = state_obj.sort_sum_by_index()
	# state_obj.print_info()

	state_obj = state_obj.sort_index_by_limit()
	# state_obj.print_info()
	# print("+++++++++++++++++++ DEBUG RESULTS ++++++++++++++++++++++<<<<<<<")
	# state_obj.print_info()
	# print("+++++++++++++++++++ DEBUG RESULTS ++++++++++++++++++++++")


	return state_obj










