#    (C) Copyright 2018 Anthony D. Dutoi and Yuhong Liu
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
import sys
import numpy
from copy import deepcopy
# from qode.many_body.fast_operator.fast_operator import cc_operator, init_K\
#               init_Ex, init_Dx, init_Fo, init_Fv, init_ExEx, init_ExFo, init_ExFv, init_ExDx, init_FoFo, init_FvFv, \
#               init_FoDx, init_FvDx, init_ExExDx, init_ExFoFo, init_ExFvFv, init_ExFoDx, init_ExFvDx, init_ExFoFv, init_ExFvFo
from qode.many_body.fast_operator.fast_operator import cc_operator
from qode.many_body.fast_operator.baker_campbell_hausdorff.autoQode.c_caller import ExFvDx_ExEx, ExDx_Ex, FvDx_Ex, DxDx_Ex, ExFv_Ex, ExFoFo_Ex, Dx_ExEx, ExFoDx_Ex, \
                   Fv_ExEx, FoDx_Ex, Dx_Ex, Fo_ExEx, FoFo_ExEx, ExFvDx_Ex, Fv_Ex, FvDx_ExEx, Fo_Ex, ExExDx_Ex, FvFv_Ex, ExFoDx_ExEx, FoFo_Ex, ExDx_ExEx, FoDx_ExEx, \
                   DxDx_ExEx, FvFv_ExEx, ExFvFv_Ex, ExFo_Ex

numpy.set_printoptions(threshold=sys.maxsize)

def init_X1(cc_op_obj, No, Nv):
	cc_op_obj.K      = numpy.zeros(1)
	cc_op_obj.Ex     = numpy.zeros(No*Nv)
	cc_op_obj.Fo     = numpy.zeros(No*No)
	cc_op_obj.Fv     = numpy.zeros(Nv*Nv)
	cc_op_obj.Dx     = numpy.zeros(Nv*No)
	cc_op_obj.ExEx   = numpy.zeros((No*Nv)**2)
	cc_op_obj.ExFo   = numpy.zeros(Nv*No**3)
	cc_op_obj.ExFv   = numpy.zeros(No*Nv**3)
	cc_op_obj.ExDx   = numpy.zeros((No*Nv)**2)
	cc_op_obj.FoFo   = numpy.zeros(No**4)
	cc_op_obj.FvFv   = numpy.zeros(Nv**4)
	cc_op_obj.FoDx   = numpy.zeros(Nv*No**3)
	cc_op_obj.FvDx   = numpy.zeros(No*Nv**3)
	cc_op_obj.FoFv   = numpy.zeros((No*Nv)**2)
	# cc_op_obj.DxDx   = numpy.zeros(((No*Nv)**2))
	cc_op_obj.ExFoFo = numpy.zeros(No**5*Nv)
	cc_op_obj.ExFvFv = numpy.zeros(Nv**5*No)
	cc_op_obj.ExFoDx = numpy.zeros((No*Nv)**2*No**2)
	cc_op_obj.ExFvDx = numpy.zeros((No*Nv)**2*Nv**2)
	cc_op_obj.ExExDx = numpy.zeros((No*Nv)**3)
	return cc_op_obj

def init_X2(cc_op_obj, No, Nv):
	cc_op_obj.K      = numpy.zeros(1)
	cc_op_obj.Ex     = numpy.zeros(No*Nv)
	cc_op_obj.Fo     = numpy.zeros(No*No)
	cc_op_obj.Fv     = numpy.zeros(Nv*Nv)
	# cc_op_obj.Dx     = numpy.zeros((Nv*No))
	cc_op_obj.ExEx   = numpy.zeros((No*Nv)**2)
	cc_op_obj.ExFo   = numpy.zeros(Nv*No**3)
	cc_op_obj.ExFv   = numpy.zeros(No*Nv**3)
	cc_op_obj.ExDx   = numpy.zeros((No*Nv)**2)
	cc_op_obj.FoFo   = numpy.zeros(No**4)
	cc_op_obj.FvFv   = numpy.zeros(Nv**4)
	cc_op_obj.FoFv   = numpy.zeros(No*Nv)**2
	# cc_op_obj.FoDx   = numpy.zeros(Nv*No**3)
	# cc_op_obj.FvDx   = numpy.zeros(No*Nv**3)
	# cc_op_obj.DxDx   = numpy.zeros((No*Nv)**2)
	cc_op_obj.ExFoFo = numpy.zeros(No**5*Nv)
	cc_op_obj.ExFvFv = numpy.zeros(Nv**5*No)
	# cc_op_obj.ExFoDx = numpy.zeros((No*Nv)**2*No**2)
	# cc_op_obj.ExFvDx = numpy.zeros(((No*Nv)**2*Nv**2))
	cc_op_obj.ExExDx = numpy.zeros(((No*Nv)**3))
	return cc_op_obj

def init_X3(cc_op_obj, No, Nv):
	# cc_op_obj.K      = numpy.zeros(1)
	cc_op_obj.Ex     = numpy.zeros((No*Nv))
	# cc_op_obj.Fo     = numpy.zeros((No*No))
	# cc_op_obj.Fv     = numpy.zeros((Nv*Nv))
	# cc_op_obj.Dx     = numpy.zeros((Nv*No))
	cc_op_obj.ExEx   = numpy.zeros(((No*Nv)**2))
	cc_op_obj.ExFo   = numpy.zeros((Nv*No**3))
	cc_op_obj.ExFv   = numpy.zeros((No*Nv**3))
	# cc_op_obj.ExDx   = numpy.zeros(((No*Nv)**2))
	# cc_op_obj.FoFo   = numpy.zeros((No**4))
	# cc_op_obj.FvFv   = numpy.zeros((Nv**4))
	# cc_op_obj.FoDx   = numpy.zeros((Nv*No**3))
	# cc_op_obj.FvDx   = numpy.zeros((No*Nv**3))
	# cc_op_obj.DxDx   = numpy.zeros(((No*Nv)**2))
	# cc_op_obj.ExFoFo = numpy.zeros((No**5*Nv))
	# cc_op_obj.ExFvFv = numpy.zeros((Nv**5*No))
	# cc_op_obj.ExFoDx = numpy.zeros(((No*Nv)**2*No**2))
	# cc_op_obj.ExFvDx = numpy.zeros(((No*Nv)**2*Nv**2))
	cc_op_obj.ExExDx = numpy.zeros(((No*Nv)**3))
	return cc_op_obj

def init_X4(cc_op_obj, No, Nv):
	#
	# Only ExEx
	cc_op_obj.ExEx   = numpy.zeros(((No*Nv)**2))
	return cc_op_obj


def init_Omega_block_zero(cc_op_obj, No, Nv):
	cc_op_obj.K      = numpy.zeros(1)
	cc_op_obj.Ex     = numpy.zeros((No*Nv))	
	cc_op_obj.ExEx   = numpy.zeros(((No*Nv)**2))
	return cc_op_obj

def scale_cluster_op(cc_op_obj, factor):
	obj_dict = cc_op_obj.__dict__
	for key in obj_dict.keys():
		if obj_dict[key] != None:
			obj_dict[key] = obj_dict[key] * factor
	return cc_op_obj

def add_Ex_block_ops(destination_op_obj, cc_op_obj):
	# __dict__ are shallow copies
	#
	dest_dict  = destination_op_obj.__dict__
	cc_op_dict = cc_op_obj.__dict__
	#
	for key in ['K', 'Ex', 'ExEx']:
		if cc_op_dict[key] != None:
			if dest_dict[key] == None:
				dest_dict[key] = deepcopy(cc_op_dict[key])
			else:
				dest_dict[key] = dest_dict[key] + cc_op_dict[key]
	return destination_op_obj



dict_Xn = {'1': init_X1, '2': init_X2, '3': init_X3, '4': init_X4}
dict_func = { 'ExFvDx_ExEx':ExFvDx_ExEx, 'ExDx_Ex':ExDx_Ex, 'FvDx_Ex':FvDx_Ex, 'DxDx_Ex':DxDx_Ex, 'ExFv_Ex':ExFv_Ex,\
              'ExFoFo_Ex':ExFoFo_Ex, 'Dx_ExEx':Dx_ExEx, 'ExFoDx_Ex':ExFoDx_Ex, 'Fv_ExEx':Fv_ExEx, 'FoDx_Ex':FoDx_Ex,\
	          'Dx_Ex':Dx_Ex, 'Fo_ExEx':Fo_ExEx, 'FoFo_ExEx':FoFo_ExEx, 'ExFvDx_Ex':ExFvDx_Ex, 'Fv_Ex':Fv_Ex, 'FvDx_ExEx':FvDx_ExEx,\
    	      'Fo_Ex':Fo_Ex, 'ExExDx_Ex':ExExDx_Ex, 'FvFv_Ex':FvFv_Ex, 'ExFoDx_ExEx':ExFoDx_ExEx, 'FoFo_Ex':FoFo_Ex, 'ExDx_ExEx':ExDx_ExEx,\
        	  'FoDx_ExEx':FoDx_ExEx, 'DxDx_ExEx':DxDx_ExEx, 'FvFv_ExEx':FvFv_ExEx, 'ExFvFv_Ex':ExFvFv_Ex, 'ExFo_Ex':ExFo_Ex }

def arraySizeOf(cc_op_obj):
	mem_in_bytes = 0
	for value in cc_op_obj.__dict__.values():
		if isinstance(value, numpy.ndarray):
			mem_in_bytes += value.nbytes
	print("\tTOTAL MEMORY USAGE = %d Bytes = %f MB = %f GB"  %(mem_in_bytes, mem_in_bytes/1024**2, mem_in_bytes/1024**3))
	if mem_in_bytes > 15 * 1024**3:
		raise AssertionError

def print_amps(X_obj):
	print('========================================')
	X_dict = X_obj.__dict__
	for key in sorted(X_dict.keys()):
		if '_' not in key and len(key) < 5:
			print(key)
			print(X_dict[key])
	print('========================================')



def commute(X, T, nth, No, Nv, resources):
	X_new = dict_Xn[ str(nth) ]( cc_operator(), No, Nv )
	X_dict_old  = X.__dict__
	T_dict      = T.__dict__
	for key, value in dict_func.items():
		key_X, key_T = key.split('_')
		if X_dict_old[key_X] != None and T_dict[key_T] != None:
			block_X = X_dict_old[key_X]
			block_T = T_dict[key_T]
			X_new = dict_func[key]( X_new, block_X, block_T, No, Nv, resources.n_cores )
	return X_new


def make_omega(H,T, No, Nv, resources):
	omega = init_Omega_block_zero(cc_operator(), No, Nv)
	X = H   # X_0 = H
	omega = add_Ex_block_ops(omega, X)
	for n in [1,2,3,4]:
		X     = commute(X, T, n, No, Nv, resources)
		X     = scale_cluster_op(X, 1.0/n)
		omega = add_Ex_block_ops(omega, X)
	return omega


class BakerCampbellHausdorff(object):
	""" This class is a wrapper for the nested-operator-specific dirty-work.  It has an interface needed by the coupled_cluster module, but hides the internal nature of the operator.  """
	def __init__(self, No, Nv, resources):
		self.No = No
		self.Nv = Nv
		self.resources = resources
	
	def computeOmega(self, H, T, textlog):	# call must have this signature!   Computes excitation part (0th, 1st and 2nd order) of HBar
		return make_omega(H, T, self.No, self.Nv, self.resources)

	@staticmethod
	def Energy(omega, textlog):		# call must have this signature!   Retrieves energy from excitation part (0th order component) of HBar
		return omega.K[0]

