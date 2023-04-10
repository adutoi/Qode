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
from . import config
import numpy as np
#
#
#
#
#

class state(object):
	def __init__(self, configs, basis=None):
		self.configs    = configs
		self.num_config = len(self.configs)
		self.coeffs     = np.zeros(( self.num_config ))  # Numpy 1D Array
		if basis is not None:  self.coeffs[basis] = 1.
	def add_config(self, config_obj):
		#
		# Add a single configuration to this state object
		#
		if not config_obj.is_null():
			config = config_obj.get_config()
			#
			# Loop over all config inside the configs to look for a match, 
			# then add the coeff to that slot
			#
			# if config not found in configs, it will be ignored.
			#
			keep_checking = True
			i = 0
			while i < self.num_config and keep_checking:
				if config == self.configs[i]:
					self.coeffs[i] += config_obj.coeff
					keep_checking   = False
				else:
					i += 1
		else:
			# Null config can be ignored.
			pass	
		#
		#
		#
	def scale(self, value):
		self.coeffs = self.coeffs * value

	def increment(self, state_obj, scaling_factor=1):
		# Scale and Add state_obj to self object. 
		#
		state_obj.scale(scaling_factor)
		self.coeffs = self.get_coeffs() + state_obj.get_coeffs()

	def update_coeffs(self, new_coeffs):
		# new_coeffs is Numpy 1D Array
		#
		if new_coeffs.shape[0] == self.num_config:
			self.coeffs = copy.deepcopy(new_coeffs)
		else:
			raise ValueError

	def get_coeffs(self):
		return copy.deepcopy(self.coeffs)

	def get_configs(self):
		return copy.deepcopy(self.configs)

	def print_info(self):
		for i in range(self.num_config):
			print("Coeff =", self.coeffs[i], "Config =", self.configs[i])





def operate(op_list, state_obj):
	#
	# op_list creation/annihilatio operators/index in normal written order
	# actual loop in a reversed order.
	#
	num_operator = len(op_list)
	this_state   = state(state_obj.get_configs())
	#
	# 
	for this_coeff, this_config in zip(state_obj.get_coeffs(), state_obj.get_configs()):
		temp_config = config.configuration(this_config, this_coeff)
		# Loop over all operators to make sure not to lose any intermediate configs
		for op, idx in reversed(op_list):
			temp_config  = op(idx, temp_config)
		this_state.add_config(temp_config)
	#
	#
	return this_state





def _add_state(state_obj1, state_obj2):
	new_state = copy.deepcopy(state_obj1)
	# NO SERIOUS CHECKING IS DONE HERE, USE WITH CAUTIONS
	#
	if state_obj1.get_coeffs().shape == state_obj2.get_coeffs().shape:
		new_state.update_coeffs( copy.deepcopy( state_obj1.get_coeffs() + state_obj2.get_coeffs() ) )
	else:
		raise ValueError
	return new_state




def dot(state_obj1, state_obj2):
	dot_prod = 0.0
	size = len(state_obj1.get_configs())
	# print("vector size =", size)
	for i in range(size):
		dot_prod += state_obj1.coeffs[i] * state_obj2.coeffs[i]
	return dot_prod





def scaled(state_obj, value):
	scaled_state = copy.deepcopy(state_obj)
	scaled_state.scale(value)
	return scaled_state













