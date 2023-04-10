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
from base import indexed_base, copy
from scalars import condense_deltas, Coeff, delta, minus, zero
from containers import container, add, mult, comm, do_commutators, reorder_ops, reorder_pdts, nest_summations, resolve_deltas, expand, sum_case, trim_high_order
from operators import Operator
from indices import Index, Indices, case_indices, nest_lists, IndexRelationInconsistency



class Sum(container):
	def __init__(self, *index_sets):
		container.__init__(self)
		self.index_sets = list(index_sets)		#shallow copy
		self.summand = None
	def __call__(self, summand):
		self.summand = summand
		return self
	def _pass_through_helper(self, function, index_replacements):
		#print("hi")
		#for index_set in self.index_sets: print(len(index_set()))
		index_sets = [index_set.copy(index_replacements) for index_set in self.index_sets]
		old_indices = [index for index_set in self.index_sets for index in index_set()]		# flat list of old indices to replace
		new_indices = [index for index_set in      index_sets for index in index_set()]		# flat list of new indices to replace old indices
		new_index_replacements = list(zip(old_indices,new_indices))				# list of (old,new) index pairs for replacement ...
		new_index_replacements += index_replacements						# ... concatenated with replacements of indices external to this sum
		if isinstance(self.summand,container):  summand = function(self.summand, new_index_replacements)
		else:                                   summand =     copy(self.summand, new_index_replacements)
		return index_sets, summand
	def _pass_through(self, function, index_replacements):
		index_sets, summand = self._pass_through_helper(function, index_replacements)
		return Sum(*index_sets)(summand)
	def _do_commutators(self, index_replacements):
		#print("                                                 COMMUTE")
		return self._pass_through(do_commutators, index_replacements)
	def _reorder_ops(self, index_replacements):      return self._pass_through(reorder_ops, index_replacements)
	def _reorder_pdts(self, index_replacements):
		#print("XXXXXXXXX")
		return self._pass_through(reorder_pdts, index_replacements)
		#print("---------")
	def _trim_high_order(self, index_replacements):  return self._pass_through(trim_high_order, index_replacements)
	def _copy(self, index_replacements):             return self._pass_through(copy, index_replacements)			# calls copy on container and primitive terms
	def _expand(self, index_replacements):
		#print("                                                 EXPAND")
		index_sets, summand = self._pass_through_helper(expand, index_replacements)
		if summand is zero:
			return zero
		elif isinstance(summand, add):
			terms = [Sum(*index_sets)(term) for term in summand.terms]
			return add(*terms)
		else:
			return Sum(*index_sets)(summand)
	def _sum_case(self, index_replacements):
		#print("                                                 CASE")
		index_sets, summand = self._pass_through_helper(copy, index_replacements)		# unusual case where we insist on working outside to inside, so just pass "copy" through here
		list_of_cased_sets = []
		for index_set in index_sets:  list_of_cased_sets += [case_indices(index_set)]		# each set now represented by a list of cases
		all_cases = nest_lists(list_of_cased_sets)						# tensor pdt of cases from each set
		terms = []
		#print("number of cases=",len(all_cases))
		#for case in all_cases:
		#	print("case")
		#	the_indices = [index for index_set in case       for index in index_set()]
		#	for index in the_indices:
		#		print("         "+str(index.relationships))
		for case in all_cases:
			old_indices = [index for index_set in index_sets for index in index_set()]		# flat list of old indices to replace
			new_indices = [index for index_set in case       for index in index_set()]		# flat list of new indices to replace old indices
			new_index_replacements = list(zip(old_indices,new_indices))				# list of (old,new) index pairs for replacement ...
			#info = "replacements:\n        "
			#for old,new in new_index_replacements:
			#	info += "(" + str(id(old)) + "," + str(id(new)) + ")   "
			#print(info)
			#print("                                                 begin subcase")
			if isinstance(summand,container):  new_summand = sum_case(summand, new_index_replacements)
			else:                              new_summand =     copy(summand, new_index_replacements)
			#print("                                                 end subcase")
			terms += [Sum(*case)(new_summand)]
		return add(*terms)
	def _nest_summations(self,index_replacements):
		index_sets, summand = self._pass_through_helper(nest_summations, index_replacements)
		if isinstance(summand, mult):
			if len(summand.factors)==1:
				coeff_pdt = summand.factors[0]
				if isinstance(coeff_pdt, mult):		# a mult in a mult
					correct_form = True
					for factor in coeff_pdt.factors:
						if not isinstance(factor, Coeff):     correct_form = False
					if correct_form:
						return Sum(*index_sets)(coeff_pdt)	# basically just unwraps double-enclosure in two mults
			if len(summand.factors)==2:
				coeff_pdt, op_pdt = summand.factors
				if isinstance(coeff_pdt, mult) and isinstance(op_pdt, mult):
					correct_form = True
					for factor in coeff_pdt.factors:
						if not isinstance(factor, Coeff):     correct_form = False
					for factor in op_pdt.factors:
						if not isinstance(factor, Operator):  correct_form = False
					if correct_form:
						outer_indices = []
						for factor in op_pdt.factors:  outer_indices += list(factor.transition())
						outer_index_sets = []
						inner_index_sets = []
						for index_set in index_sets:
							outer_index_list = []
							inner_index_list = []
							for index in index_set.indices:
								if index in outer_indices:  outer_index_list += [index]
								else:                       inner_index_list += [index]
							n = len(outer_index_list)
							if outer_index_list:
								outer_index_set = Indices(outer_index_list, bound=index_set.bound, letters=index_set.letters[:n])
								outer_index_sets += [outer_index_set]
							if inner_index_list:
								inner_index_set = Indices(inner_index_list, bound=index_set.bound, letters=index_set.letters[n:])
								inner_index_sets += [inner_index_set]
						if inner_index_sets:  inner_sum = Sum(*inner_index_sets)(coeff_pdt)
						else:                 inner_sum = coeff_pdt
						return Sum(*outer_index_sets)(mult(inner_sum,op_pdt))
		return Sum(*index_sets)(summand)		# either this function has already returned or it was not a primitive sum over the right kind of product
	def _resolve_deltas(self,index_replacements):
		index_sets, summand = self._pass_through_helper(resolve_deltas, index_replacements)
		if isinstance(summand, mult):
			if len(summand.factors)==2:									# this block same as one below with any reference to op_pdt removed ... too lazy to condense now
				delta_pdt, coeff_pdt = summand.factors
				if isinstance(delta_pdt, mult) and isinstance(coeff_pdt, mult):
					correct_form = True
					for factor in delta_pdt.factors:
						if not isinstance(factor, delta):     correct_form = False
					for factor in coeff_pdt.factors:
						if not isinstance(factor, Coeff):     correct_form = False
					if correct_form:
						delta_pairs = [factor.bottom for factor in delta_pdt.factors]
						for index_set in index_sets:
							for index in index_set.indices:
								for relate,other in index.relationships:
									if relate=="=":
										delta_pairs += [(index,other)]
						equality_groups = condense_deltas(delta_pairs)
						for parents in equality_groups:
							try:
								child,(coeff_pdt,) = Index.inherit_relationships(parents, (coeff_pdt,))		# naming made more sense when child was returned and replaced from outside
							except IndexRelationInconsistency:
								return zero
						return Sum(*index_sets)(coeff_pdt)
			elif len(summand.factors)==3:
				delta_pdt, coeff_pdt, op_pdt = summand.factors
				if isinstance(delta_pdt, mult) and isinstance(coeff_pdt, mult) and isinstance(op_pdt, mult):
					correct_form = True
					for factor in delta_pdt.factors:
						if not isinstance(factor, delta):     correct_form = False
					for factor in coeff_pdt.factors:
						if not isinstance(factor, Coeff):     correct_form = False
					for factor in op_pdt.factors:
						if not isinstance(factor, Operator):  correct_form = False
					if correct_form:
						delta_pairs = [factor.bottom for factor in delta_pdt.factors]
						for index_set in index_sets:
							for index in index_set.indices:
								for relate,other in index.relationships:
									if relate=="=":
										delta_pairs += [(index,other)]
						equality_groups = condense_deltas(delta_pairs)
						for parents in equality_groups:
							try:
								child,(coeff_pdt,op_pdt) = Index.inherit_relationships(parents, (coeff_pdt,op_pdt))
							except IndexRelationInconsistency:
								return zero
						return Sum(*index_sets)(mult(coeff_pdt,op_pdt))
		# if that did not work, go after deltas buried only in sum conditions
		delta_pairs = []
		for index_set in index_sets:
			for index in index_set.indices:
				for relate,other in index.relationships:
					if relate=="=":
						delta_pairs += [(index,other)]
		equality_groups = condense_deltas(delta_pairs)
		for parents in equality_groups:
			try:
				child,(summand,) = Index.inherit_relationships(parents, (summand,))
			except IndexRelationInconsistency:
				return zero
		return Sum(*index_sets)(summand)
	@staticmethod
	def _index_ordering(redundant_ordered_indices, index_sets, should_match=False):
		#print("==")
		#print(len(redundant_ordered_indices))
		sum_indices = []				# all the indices that occur under the summation sign
		for index_set in index_sets:
			#print(index_set.bound,len(index_set.indices))
			sum_indices += index_set.indices
		#print("sum indices", len(sum_indices))
		#
		nonredundant = []				# unique summand indices in order of appearance withouth regard to which set they are in
		for index in redundant_ordered_indices:
			if index not in nonredundant:  nonredundant += [index]		# keeping first occurance is optimal for loop ordering (otherwise get big memory-space jumps when inner index changes)
		#print("nonredundant", len(nonredundant))
		checked = []					# unique indices that occur both in the summand and under the summation sign
		for index in nonredundant:
			if index in sum_indices:  checked += [index]
			elif should_match:        raise Exception("assertion that summation and summand indices match is not correct")		# summands of inner sums need not obey this
		N = len(checked)
		#print("checked", len(checked))
		if len(sum_indices)!=N:  raise Exception("mismatch of summation and summand indices")
		#
		ordered_char_bank = [None]*N			# the "target" to populate in order to print (might get discarded if next condition not true)
		for index_set in index_sets:
			ordered_set_list = []			# for each set make a list of they order those will appear in ...
			positions = []				# ... and keep track of where they will appear in the global ordering
			for n,index in enumerate(checked):	# Now populate the above two lists
				if index in index_set.indices:
					ordered_set_list += [index]
					positions += [n]
			new_index_set = Indices(ordered_set_list, bound=index_set.bound, letters=index_set.letters)	# make a shadow index_set with ordered indices ...
			char_bank = new_index_set.assign_letters()							# ... and resolve the characters ...
			for n,chars in zip(positions,char_bank):  ordered_char_bank[n] = chars				# ... and put those into the globally ordered list of characters (the target)
		return ordered_char_bank
	@staticmethod
	def _string_helper(char_bank, summand):
		string = "\\left\\{\\sum_{"
		upper = ""
		sep1 = ""
		for chars in char_bank:
			string += sep1 + chars[0]
			upper  += sep1 + chars[1][0]+"_\\text{"+chars[1][1]+"}" 	# pretty strict assumption on how this is formatted
			sep1 = ", "
			sep2 = ""
			for relate,other in chars[2:]:
				string += sep2 + relate + other
				sep2 = "\\text{\\&}"
		string += "}^{" + upper + "} " + str(summand) + "\\right\\}"
		return string
	def __str__(self):
		#print("                                                 STRING")
		index_sets, summand = self._pass_through_helper(copy, [])
		#print("                                                 helper done")
		if isinstance(summand, mult):
			if len(summand.factors)==2:
				_, op_pdt = summand.factors
				if isinstance(op_pdt, mult):
					correct_form = True
					for factor in op_pdt.factors:
						if not isinstance(factor, Operator):  correct_form = False
					if correct_form:
						op_indices = []
						for factor in op_pdt.factors:  op_indices += [factor.bottom, factor.top]
						ordered_char_bank = Sum._index_ordering(op_indices, index_sets, should_match=True)
						return Sum._string_helper(ordered_char_bank, summand)
		# if we did not hit the return above, then we try this option
		if isinstance(summand, mult):
			if (len(summand.factors)==1) or (len(summand.factors)==2 and summand.factors[0] is minus):
				h_or_V = summand.factors[-1]
				possible_minus = summand.factors[0] 	# possibly the same as h_or_V (ok, b/c just testing that it is Coeff)
				correct_form = True
				for factor in (possible_minus, h_or_V):
					if not isinstance(factor, Coeff):  correct_form = False
				if correct_form:
					coeff_indices = []				# all the index pairs of the Coeffs, with lowest order indices at the end
					for pair in range(len(h_or_V.top)):		# bottom and top of each have same length
						p = -1-pair
						coeff_indices = [h_or_V.top[p]]    + coeff_indices
						coeff_indices = [h_or_V.bottom[p]] + coeff_indices
					ordered_char_bank = Sum._index_ordering(coeff_indices, index_sets)
					return Sum._string_helper(ordered_char_bank, summand)
			elif len(summand.factors) in (2,3):
				X,T = summand.factors[-2], summand.factors[-1]
				possible_minus = summand.factors[0] 	# possibly the same as X (ok, b/c just testing that it is Coeff)
				correct_form = True
				for factor in (possible_minus, X, T):
					if not isinstance(factor, Coeff):  correct_form = False
				if correct_form:
					coeff_indices = []					# all the index pairs of the Coeffs, with lowest order indices at the end
					for pair in range(max(len(X.top),len(T.top))):		# bottom and top of each have same length
						p = -1-pair
						if pair<len(T.top):     coeff_indices = [T.top[p]]    + coeff_indices
						if pair<len(X.top):     coeff_indices = [X.top[p]]    + coeff_indices
						if pair<len(T.bottom):  coeff_indices = [T.bottom[p]] + coeff_indices
						if pair<len(X.bottom):  coeff_indices = [X.bottom[p]] + coeff_indices
					ordered_char_bank = Sum._index_ordering(coeff_indices, index_sets)
					return Sum._string_helper(ordered_char_bank, summand)
		# if we did not hit return on either of the above, do a generic unordered print of the summation indices
		overall_char_bank = []
		for index_set in index_sets:
			char_bank = index_set.assign_letters()
			for chars in char_bank:  overall_char_bank += [chars]
		return Sum._string_helper(overall_char_bank, summand)
	def multiline(self,_):
		return self
	def code(self, func_num, prefix):
		# print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		from c_generator import full_case_gen, two_mult_gen, three_coeff_gen, one_term_gen
		outmost_index_sets, outmost_summand = self._pass_through_helper(copy, [])
		#
		# May have to send them to the detailed handlers or do it here!! Assign letters of Molecules is 
		# DIFFERENT THAN the assignment of indices
		#
		outmost_char_bank = []
		for index_set in outmost_index_sets:
			char_bank = index_set.assign_letters()
			for chars in char_bank:  
				outmost_char_bank += [chars]
		#
		#
			index_sets, summand = outmost_summand._pass_through_helper(copy, [])
			if isinstance(summand, mult):
				num_factors = len(summand.factors)
				if num_factors == 2:
					# print("TWO-TERM CASE")
					if isinstance(summand.factors[0], Sum) and isinstance(summand.factors[1], mult):
						# print("The Full Pack Case")
						full_case_gen(outmost_char_bank, index_sets, summand, func_num, prefix)
					elif isinstance(summand.factors[0], mult) and isinstance(summand.factors[1], mult):
						# print("The Two Mult Case")
						two_mult_gen(outmost_char_bank, index_sets, summand, func_num, prefix)
					elif isinstance(summand.factors[0], Coeff) and isinstance(summand.factors[1], Coeff):
						# print("XT Case")
						three_coeff_gen(outmost_char_bank, index_sets, summand, func_num, prefix)
					else:
						raise AssertionError
				elif num_factors == 3:
					# The (-1)XT case
					# print("(-1)XT CASE")
					if isinstance(summand.factors[0], Coeff) and isinstance(summand.factors[1], Coeff) and \
					       isinstance(summand.factors[2], Coeff):
						three_coeff_gen(outmost_char_bank, index_sets, summand, func_num, prefix)
					else:
						raise AssertionError
				elif num_factors == 1:
					one_term_gen(outmost_char_bank, index_sets, summand, func_num, prefix)
				else:
					raise AssertionError
