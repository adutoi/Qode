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
import numpy
from .local_operator import cc_operator, init_K, init_blocks



def dot_cluster_ops(obj1, obj2):
	dict_obj1 = obj1.__dict__
	dict_obj2 = obj2.__dict__
	result = 0.0
	for key in dict_obj1.keys():
		if '_' not in key:
			# print("dot func key =",key)
			if dict_obj1[key] is not None and dict_obj2[key] is not None:
				for i in range( dict_obj1[key].shape[0] ):
					result += dict_obj1[key][i] * dict_obj2[key][i]
	# print("DOT PRODUCT RESULT =", result)
	return result

def add_cluster_ops(destination_op_obj, cc_op_obj, scalar):
	dest_dict  = destination_op_obj.__dict__
	cc_op_dict = cc_op_obj.__dict__
	#
	for key, array in cc_op_dict.items():
		if '_' not in key:
			if array is not None:
				if dest_dict[key] is None:
					dest_dict[key]  = deepcopy(array)
					dest_dict[key] *= scalar
				else:
					dest_dict[key] += scalar * array

def scale_cluster_op(cc_obj, factor):
	if factor == 1.0:
		pass
	else:
		obj_dict = cc_obj.__dict__
		for key in obj_dict.keys():
			if '_' not in key:
				if obj_dict[key] is not None: 
					obj_dict[key] *= factor

# def two_idx(i,j):  # Molecule i and j, find the absolute index in the storage.
# 	abs_idx = 0
# 	for x in range(i):
# 		abs_idx += x + 1
# 	abs_idx += j
# 	return abs_idx


def contiguous_deepcopy(cc_obj):
	input_dict   = cc_obj.__dict__
	new_obj      = cc_operator()
	new_obj_dict = new_obj.__dict__
	for key in input_dict.keys():
		if input_dict[key] is not None:
			if '_' not in key: # Storage Blocks
				new_obj_dict[key] = numpy.ascontiguousarray( input_dict[key].copy(), dtype=numpy.float64 ) 
			elif'starter' in key:  # Starters-position array
				new_obj_dict[key] = numpy.ascontiguousarray( input_dict[key].copy(), dtype=numpy.int64 )
			else:
				new_obj_dict[key] = input_dict[key]        # They are int64 types
		else:
			new_obj_dict[key] = None
	return new_obj

class space_traits_class(object):
	def __init__(self):
		self.field = numpy.float64
	@staticmethod
	def dot(v,w):
		return dot_cluster_ops(v,w)
	@staticmethod
	def add_to(v,w,n=1.0):
		add_cluster_ops(v,w,n)
	@staticmethod
	def scale(n,v):
		scale_cluster_op(v,n)
	@staticmethod
	def copy(v):
		return contiguous_deepcopy(v)
	@staticmethod
	def act_on_vec(op,v):
		return op(v)	# This is a single tasker, only here for acting orb energy denominators
	#
	def check_member(self,v):
		pass
	def check_lin_op(self,op):
		return False
	@staticmethod
	def function_on_diags(func,op):
		raise NotImplementedError
	@staticmethod
	def back_act_on_vec(v,op):
		raise NotImplementedError
	@staticmethod
	def diagonal(op):
		raise NotImplementedError

	# occ_orb_idx = [ i for i in range(num_alpha_elec)] + [ i+num_spatial_orb for i in range(num_beta_elec)]
	# vrt_orb_idx = [ i+num_alpha_elec for i in range(num_spatial_orb - num_alpha_elec) ] + [ i+num_spatial_orb+num_beta_elec for i in range(num_spatial_orb - num_beta_elec) ]

class orbital_energy(object):
	"""Class for modifying Omega operators with orbital energy differences"""
	def __init__(self, mol_eigen_energy, rec_num_states):
		self.eigen_energy = mol_eigen_energy
		self.rec = rec_num_states
	def __call__(self, omega_obj):
		inverted_omega = contiguous_deepcopy(omega_obj)
		inverted_omega.K[0] = 0.0
		# print("PRINTING FROM INVF OBJ--------------------------------------------------------------")
		# print(inverted_omega.Ex)
		# print(inverted_omega.ExEx)
		# print('------------------------------------------------------------------------------------')
		num_mol = len(self.rec)
		No = 1
		Nv = [ num - 1 for num in self.rec ]
		#
		# Ex divisions:
		for n in range(num_mol):
			for a in range(Nv[n]):
				# print("InvF Ex", self.eigen_energy[n][ a+1 ] -  self.eigen_energy[n][0])
				inverted_omega.Ex[ inverted_omega.Ex_starters[n] + a ] /=  (self.eigen_energy[n][ a+1 ] -  self.eigen_energy[n][0])
				# print(inverted_omega.Ex[ inverted_omega.Ex_starters[n] + a ] )
		#
		# ExEx Divisions:  ExEx[i,a,j,b] -->  <ij|V|ab> / ( Em[a+1] + En[b+1] - Em[0] - En[0] )
		for m in range(num_mol):
			for n in range(num_mol):  # Lower Triangle of the matrix
				if m != n:
					for a in range(Nv[m]):
						for b in range(Nv[n]):
							# print("InvF ExEx", self.eigen_energy[m][a+1] + self.eigen_energy[n][b+1] -  self.eigen_energy[m][0] - self.eigen_energy[n][0])
							inverted_omega.ExEx[ inverted_omega.ExEx_starters[ m*num_mol+n ] + a*Nv[n] + b]  /= \
									       (self.eigen_energy[m][a+1] + self.eigen_energy[n][b+1] -  self.eigen_energy[m][0] - self.eigen_energy[n][0] )
							# print(inverted_omega.ExEx[ inverted_omega.ExEx_starters[ two_idx(m,n) ] + a*Nv[n] + b])
		return inverted_omega



# for now, this acts like a singleton, but, eventually, I will want to put in some kind of dimension checking, in which case this will be deleted, the above will have "_class" removed, and that will be used

space_traits = space_traits_class()
