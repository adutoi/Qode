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
import pickle
import numpy as np
from math import sqrt

# 1) For printing errors
#
def get_norm(cc_obj):
	dict_obj = cc_obj.__dict__
	norm2    = 0.0
	for array in  dict_obj.values():
		if array != None:
			for item in array:
				norm2 += item**2
	return sqrt(norm2)

def err(vec1, vec2):
	diff_vec  = add_cluster_ops(deepcopy(vec1), scale_cluster_op(deepcopy(vec2), -1.0))
	diff = get_norm(diff_vec) 
	norm_vec2 = get_norm(vec2)
	return diff/norm_vec2

def print_err(vec1, vec2):
	print("err =", err(vec1,vec2), "norm =", get_norm(vec2))



# 2) Indexing Utilities
#
def abs_to_relative_pos(this_index, occ_orbs, vrt_orbs, No, Nv):
	# print(occ_orbs,vrt_orbs,this_index)
	if this_index in vrt_orbs:
		return vrt_orbs.index(this_index), Nv
	elif this_index in occ_orbs:
		return occ_orbs.index(this_index), No
	else:
		raise AssertionError


def index_to_position(indices, multipliers):
	if len(indices) % 2 == 0 and len(indices) == len(multipliers):
		position = 0
		for i in range(len(indices)):
			temp = indices[i]
			for j in range(i+1, len(indices)):
				temp *= multipliers[j]
			position += temp
		return position
	else:
		raise AssertionError




# 3) For loading pickled data

def loadpickled(filename):
	return pickle.load(open(filename,'rb'))

def add_value_to_block(key, indices, multipliers, value ,cc_op_obj, No, Nv):
	dict_dim = { 'K': 1, 'Ex': No*Nv, 'Fo':No*No, 'Fv': Nv*Nv, 'Dx': Nv*No, 'ExEx': (No*Nv)**2, 'ExFo': Nv*No**3, \
                 'ExFv': No*Nv**3, 'ExDx': (No*Nv)**2, 'FoFo': No**4, 'FvFv': Nv**4, 'FoDx': Nv*No**3, 'FvDx': No*Nv**3, \
                 'DxDx': (No*Nv)**2, 'ExFoFo': No**5*Nv, 'ExFvFv': Nv**5*No, 'ExFoDx': (No*Nv)**2*No**2, 'ExFvDx': (No*Nv)**2*Nv**2, \
                 'ExExDx': (No*Nv)**3 }
	current_dict = cc_op_obj.__dict__
	current_keys = current_dict.keys()
	if key in current_keys: # sanity check
		if current_dict[key] == None: # NoneType Not initialized
			if indices == None:
				# Single Value K
				current_dict[key] = np.array([value])
			else:
				current_dict[key] = np.zeros((dict_dim[key]))
				current_dict[key][index_to_position(indices, multipliers)] = value
		else:   # initialized...
			if indices == None: 
				# Single Value K
				current_dict[key][0] = value
			else:
				current_dict[key][index_to_position(indices, multipliers)] = value 



def load_pickle_into_blocks(unpickled_data, cc_op_obj, occ_orbs, vrt_orbs):
	# unpicked_data is in Tony's format!
	#
	occ_orbs = sorted(occ_orbs)
	vrt_orbs = sorted(vrt_orbs)
	No = len(occ_orbs)
	Nv = len(vrt_orbs)
	# 
	for item in unpickled_data:
		identifier, value = item
		if len(identifier) == 0:
			key     = 'K'
			indices = None
			multipliers = None
		else:
			key = ''
			abs_indices = []
			for op, loc in identifier:
				key += op
				if op == 'Ex' or op == 'Fv':
					abs_indices += reversed(list(loc))
				else:
					abs_indices += list(loc)
			indices = []
			multipliers = []
			for item in abs_indices:
				rel_index, this_multiplier = abs_to_relative_pos(item, occ_orbs, vrt_orbs, No, Nv)
				indices     += [rel_index]
				multipliers += [this_multiplier]
		add_value_to_block(key, indices, multipliers, value, cc_op_obj, No, Nv)
		# if key == 'Fv':
		# 	print(key ,'new', indices, multipliers, value)
		# if key == 'ExDx':
		# 	print(key ,'new', indices, multipliers, value)
			
	return cc_op_obj


# 4) Printing Utilities
#
def print_array(cc_op_obj):
	print("++++++++++++++++++ PRINTING CLUSTER OPERATOR BLOCKS ++++++++++++++++++++++++++++++++")
	this_dict = cc_op_obj.__dict__
	keys = this_dict.keys()
	for key in sorted(keys):
		print(key, this_dict[key])
	print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

