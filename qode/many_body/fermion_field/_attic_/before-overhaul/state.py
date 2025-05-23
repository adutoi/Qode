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
from copy import deepcopy
import numpy as np
from . import config
from .config import configuration, create, annihilate, op_string	# so they may be imported directly from the state module without importing config



class state(object):
	def __init__(self, occ_strings, basis_vec=None):
		coefficients = np.zeros(( len(occ_strings) ))	# Numpy 1D Array
		if basis_vec is not None:  coefficients[basis_vec] = 1.
		self._data = ( deepcopy(occ_strings), coefficients )
	def __str__(self):
		occ_strings, coefficients = self._data
		string = ""
		for i in range(len(occ_strings)):  string += str(configuration(occ_strings[i],coefficients[i])) + "\n"
		return string
	def value(self):
		occ_strings, coefficients = self._data
		return (deepcopy(occ_strings), deepcopy(coefficients))
	def add_config(self, config):
		""" Add a single configuration to this state.  If config not represented in occ_strings, it will be ignored (projected off). """
		if config.is_not_null():
			occ_string,  coefficient  = config.value()
			occ_strings, coefficients = self._data
			# .index is better than searching the "list" manually because occ_strings may not be a list (though it should be indexable and iterable and have a length)
			try:                i = occ_strings.index(occ_string)
			except ValueError:  pass		# means that occ_string was not found in occ_strings
			else:               coefficients[i] += coefficient
	def increment(self, other, scaling_factor=1):
		""" Scale and add another state to this state """
		_, these_coefficients =  self._data
		_, those_coefficients = other._data
		these_coefficients += scaling_factor * those_coefficients
	def scale(self, scaling_factor):
		""" scale this state in place """
		_, coefficients =  self._data
		coefficients *= scaling_factor
	def dot(self,other):
		_, these_coefficients =  self._data
		_, those_coefficients = other._data
		return np.dot(these_coefficients,those_coefficients)
	def _op_string_action(self,ops):
		occ_strings, coefficients = self._data
		new_state = state(occ_strings)
		for pair in zip(occ_strings,coefficients):  new_state.add_config( ops | configuration(*pair) )
		return new_state
	############## DIRTY DIRTY DIRTY BIRDY ##############
	def update_coeffs(self, new_coeffs):
		occ_strings, coefficients = self._data
		coefficients = deepcopy(new_coeffs)
		self._data = ( occ_strings, np.array(coefficients) )

def dot(state1, state2):  return state1.dot(state2)



def resolvent(operator):
	""" Generates a function that acts the inverse of the diagonal of some operator on the state """
	def action(the_state):
		occ_strings, coefficients = the_state.value()
		new_state = state(occ_strings)
		for occ_string,coefficient in zip(occ_strings,coefficients):
			temp_state = state(occ_strings)
			temp_state.add_config(configuration(occ_string))
			diagonal = dot(temp_state,operator(temp_state))
			new_state.add_config(configuration(occ_string,coefficient/diagonal))
		return new_state
	return action
