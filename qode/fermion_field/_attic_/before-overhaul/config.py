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



class configuration(object):
	"""\
	The main job of this class is to store an occupation string and be able to resolve the action of a
	creation or annihilation operator on that string.  Since this may involve changing the phase of the
	configuration, we generalize and allow a coefficient to be passed along, since that will also be a
	convenience for functions that use this.

	In the interest of being explicit, this class expects you to test whether or not it represents the Null
	configuration from the outside, because it will raise an Exception if you try to look inside Null.
	"""
	def __init__(self, occ_string, coefficient=1):
		# occ_string = [T,F,T,...F],  ... a list of True/False to represent the occupation information.
		self._occ_string  = occ_string
		self._coefficient = coefficient
	def __str__(self):
		string = "{:+6.3f}   x   | ".format(self._coefficient)		# do for Null too (coeff=0), as cheap way of keeping justification
		if self.is_not_null():
			for i,occupancy in enumerate(self._occ_string):
				if occupancy:  string += "{} ".format(i)
				else:          string += ". "
			string += ">"
		else:
			string += "NULL VECTOR >"
		return string
	def value(self):
		if self.is_not_null():  return (copy.copy(self._occ_string),self._coefficient)
		else:                   raise  Exception("attempted access of occupation information for Null configuration")
	def create(self,index):
		""" operates on configuration in place! """
		if self.is_not_null():			# Action onto Null leaves Null.
			if self._occ_string[index]:	# If occupied, then it becomes Null.
				self._make_null()
			else:				# Not occupied.   Good, we can create on this slot!
				self._occ_string[index]  = True
				self._coefficient       *= pow(-1, self._num_occ_before(index))
	def annihilate(self,index):
		""" operates on configuration in place! """
		if self.is_not_null():			# Action onto Null leaves Null.
			if self._occ_string[index]:	# Occupied.   Good, we can annihilate on this slot!
				self._occ_string[index]  = False
				self._coefficient       *= pow(-1, self._num_occ_before(index))
			else:				# If not occupied, then it becomes Null.
				self._make_null()
	def is_not_null(self):
		return (self._occ_string is not None)
	def _make_null(self):
		self._occ_string = None
		self._coefficient   = 0
	def _num_occ_before(self, index):
		""" How many electrons are in orbitals with indices lower than the one of interest """
		n = 0
		for occupancy in self._occ_string[:index]:
			if occupancy:  n += 1
		return n
	def _op_string_action(self,ops):
		new_config = copy.deepcopy(self)
		for op in reversed(ops.operators):
			if op.is_create:  new_config.create(op.index)
			else:             new_config.annihilate(op.index)
		return new_config
	def __or__(self,other):
		if self._occ_string==other._occ_string:  return self._coefficient * other._coefficient
		else:  return 0



class create(object):
	def __init__(self,index):
		self.index = index
		self.is_create     = True
		self.is_annihilate = False
	def __call__(self,the_config):
		new_config = copy.deepcopy(the_config)
		new_config.create(self.index)
		return new_config

class annihilate(object):
	def __init__(self,index):
		self.index = index
		self.is_create     = False
		self.is_annihilate = True
	def __call__(self,the_config):
		new_config = copy.deepcopy(the_config)
		new_config.annihilate(self.index)
		return new_config

class op_string(object):
	""" historically this was defined in state.py as a means by which to safely act number-conserving strings on a state, but then I wanted it here to ... to avoid bugs, the old interface is entirely deprecated (used to used __call__) """
	def __init__(self,*operators):
		self.operators = operators
	def __or__(self,config_or_state):
		return config_or_state._op_string_action(self)
