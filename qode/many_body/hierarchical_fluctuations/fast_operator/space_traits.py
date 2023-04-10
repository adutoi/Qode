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
from copy import deepcopy

def dot_cluster_ops(obj1, obj2):
	dict_obj1 = obj1.__dict__
	dict_obj2 = obj2.__dict__
	result = 0.0
	for key in dict_obj1.keys():
		if dict_obj1[key] != None and dict_obj2[key] != None:
			array1 = dict_obj1[key]
			array2 = dict_obj2[key]
			for i in range( array1.shape[0] ):
				result += array1[i] * array2[i]
	# print("RESULT =", result)
	return result

def add_cluster_ops(destination_op_obj, cc_op_obj):
	dest_dict  = destination_op_obj.__dict__
	cc_op_dict = cc_op_obj.__dict__
	#
	for key, array in cc_op_dict.items():
		if array != None:
			if dest_dict[key] == None:
				dest_dict[key] = deepcopy(array)
			else:
				dest_dict[key] = dest_dict[key] + array

def scale_cluster_op(cc_op_obj, factor):
	if factor == 1.0:
		return cc_op_obj
	else:
		obj_dict = cc_op_obj.__dict__
		for key in obj_dict.keys():
			if obj_dict[key] != None:
				obj_dict[key] = obj_dict[key] * factor
		return cc_op_obj


class space_traits_class(object):
	def __init__(self):
		self.field = numpy.float64
	@staticmethod
	def dot(v,w):
		return dot_cluster_ops(v,w)
	@staticmethod
	def add_to(v,w,n=1.0):
		add_cluster_ops( v, scale_cluster_op(w,n) )
	@staticmethod
	def scale(n,v):
		scale_cluster_op(v,n)
	@staticmethod
	def copy(v):
		return deepcopy(v)
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
	def __init__(self, alpha_beta_orb_energy, occ_orbs, vrt_orbs):
		self.No = len(occ_orbs)
		self.Nv = len(vrt_orbs)
		self.Nspin = self.No + self.Nv
		#if len(alpha_beta_orb_energy) % 2 != 0:
		#	raise AssertionError
		self.orbital_energy = alpha_beta_orb_energy
		self.occ_idx = sorted(occ_orbs)
		self.vrt_idx = sorted(vrt_orbs)
	def __call__(self, omega_obj):
		inverted_omega = deepcopy(omega_obj)
		inverted_omega.K[0] = 0.0
		for i in range(self.No):
			for a in range(self.Nv):
				inverted_omega.Ex[i*self.Nv + a] /= (self.orbital_energy[self.vrt_idx[a]] -  self.orbital_energy[self.occ_idx[i]]) #* -1.0
		for i in range(self.No):
			for a in range(self.Nv):
				for j in range(i):
					for b in range(a):
						inverted_omega.ExEx[i*self.Nv*self.No*self.Nv + a*self.No*self.Nv + j* self.Nv + b] /= \
						        (self.orbital_energy[self.vrt_idx[a]] + self.orbital_energy[self.vrt_idx[b]] -  self.orbital_energy[self.occ_idx[i]] - self.orbital_energy[self.occ_idx[j]]) #*-1.0
		return inverted_omega



# for now, this acts like a singleton, but, eventually, I will want to put in some kind of dimension checking, in which case this will be deleted, the above will have "_class" removed, and that will be used

space_traits = space_traits_class()
