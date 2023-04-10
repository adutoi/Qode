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
#
#
#
#

class configuration(object):
	def __init__(self, config, coeff=1.0):
		# config = [T,F,T,F...]
		# It's a list of True/False to represent the occupation information.
		self.config = config
		self.coeff  = coeff


	def is_null(self):
		# The VECTOR is considered as NULL if 
		# the coeff == 0 and config == None 
		#
		if self.coeff is 0 and self.config == None:
			return True
		else:
			return False


	def is_created(self, orb_index):
		# Check if this orbital of interest has been created.
		#
		if self.config[orb_index] == True:
			return True
		else:
			return False

	def make_null(self):
		#
		# This make_null set coeff to be INTEGER ZERO.
		self.config = None
		self.coeff  = 0


	def num_created_before(self, orb_index):
		num_orb_created_before = 0
		for orb in self.config[:orb_index]:
			if orb == True:
				num_orb_created_before += 1
		return num_orb_created_before



	def get_config(self):
		return copy.deepcopy(self.config)


	def get_coeff(self):
		return self.coeff


	def print_info(self):
		print("COEFF =",self.coeff,"CONFIG =", self.config)







def create(orb_index, config_obj):
	#
	#
	this_config = copy.deepcopy(config_obj) 
	if this_config.is_null():
		# Do nothing if it is Null
		pass
	#
	#
	elif not this_config.is_created(orb_index):
		# Not Created, Good, we can create on this slot!
		this_config.config[orb_index]  = True
		this_config.coeff             *= pow(-1.0, this_config.num_created_before(orb_index))
	else:
		# If Created, then it becomes Null.
		this_config.make_null()
	#
	#
	return this_config



def annihilate(orb_index, config_obj):
	#
	#
	this_config = copy.deepcopy(config_obj) 
	if this_config.is_null():
		# Do nothing if it is Null
		pass
	#
	#
	elif this_config.is_created(orb_index):
		# Not Created, Good, we can create on this slot!
		this_config.config[orb_index]  = False
		this_config.coeff             *= pow(-1.0, this_config.num_created_before(orb_index))
	else:
		# If Not Created, then it becomes Null.
		this_config.make_null()
	#
	#
	return this_config




















if __name__ == "__main__":
	pass







