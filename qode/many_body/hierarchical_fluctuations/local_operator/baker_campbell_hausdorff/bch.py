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

from ..local_operator import cc_operator, init_K, init_blocks
try:
	from .autoQode import c_caller
except ImportError:
	# This works for now.  See notes in NotImplemented/c_caller.  The point is to keep a broken import from stopping completely unrelated code from running.
	from .NotImplemented import c_caller

numpy.set_printoptions(threshold=sys.maxsize,linewidth=270)


block_dict ={	'Ex': ['OCC','VRT'],
				'Dx': ['VRT','OCC'],
				'Fo': ['OCC','OCC'],
				'Fv': ['VRT','VRT'],
				'ExEx': ['OCC','VRT','OCC','VRT'],
				'ExFo': ['OCC','VRT','OCC','OCC'],
				'ExFv': ['OCC','VRT','VRT','VRT'],
				'ExDx': ['OCC','VRT','VRT','OCC'],
				'FoFo': ['OCC','OCC','OCC','OCC'],
				'FvFv': ['VRT','VRT','VRT','VRT'],
				'FoDx': ['OCC','OCC','VRT','OCC'],
				'FvDx': ['VRT','VRT','VRT','OCC'],
				'FoFv': ['OCC','OCC','VRT','VRT'],
				'DxDx': ['VRT','OCC','VRT','OCC'],
				'ExExDx': ['OCC','VRT','OCC','VRT','VRT','OCC'],
				'ExFoFo': ['OCC','VRT','OCC','OCC','OCC','OCC'],
				'ExFvFv': ['OCC','VRT','VRT','VRT','VRT','VRT'],
				'ExFoDx': ['OCC','VRT','OCC','OCC','VRT','OCC'],
				'ExFvDx': ['OCC','VRT','VRT','VRT','VRT','OCC'],
				'ExFoFv': ['OCC','VRT','OCC','OCC','VRT','VRT'] }


# def arraySizeOf(cc_op_obj):
# 	mem_in_bytes = 0
# 	for value in cc_op_obj.__dict__.values():
# 		if isinstance(value, numpy.ndarray):
# 			mem_in_bytes += value.nbytes
# 	print("\tTOTAL MEMORY USAGE = %d Bytes = %f MB = %f GB"  %(mem_in_bytes, mem_in_bytes/1024**2, mem_in_bytes/1024**3))
# 	if mem_in_bytes > 15 * 1024**3:
# 		raise AssertionError

def init_X1(cc_obj, rec):
	# REAL USEFUL BLOCKS = ['Dx', 'Ex', 'ExDx', 'ExEx', 'ExExDx', 'ExFo', 'ExFoDx', 'ExFoFo', 'ExFoFv', 'ExFv', 'ExFvDx', 'ExFvFv', 'Fo', 'FoDx', 'FoFo', 'FoFv', 'Fv', 'FvDx', 'FvFv', 'K']
	cc_obj.K,  cc_obj.K_starters, cc_obj.K_int_dim, cc_obj.K_dim   = init_K()
	cc_obj.Ex, cc_obj.Ex_starters, cc_obj.Ex_int_dim, cc_obj.Ex_dim   = init_blocks(rec, block_dict['Ex'])
	cc_obj.Fo, cc_obj.Fo_starters, cc_obj.Fo_int_dim, cc_obj.Fo_dim   = init_blocks(rec, block_dict['Fo'])
	cc_obj.Fv, cc_obj.Fv_starters, cc_obj.Fv_int_dim, cc_obj.Fv_dim   = init_blocks(rec, block_dict['Fv'])
	cc_obj.Dx, cc_obj.Dx_starters, cc_obj.Dx_int_dim, cc_obj.Dx_dim   = init_blocks(rec, block_dict['Dx'])
	
	cc_obj.ExEx, cc_obj.ExEx_starters, cc_obj.ExEx_int_dim, cc_obj.ExEx_dim = init_blocks(rec, block_dict['ExEx'])
	cc_obj.ExFo, cc_obj.ExFo_starters, cc_obj.ExFo_int_dim, cc_obj.ExFo_dim = init_blocks(rec, block_dict['ExFo'])
	cc_obj.ExFv, cc_obj.ExFv_starters, cc_obj.ExFv_int_dim, cc_obj.ExFv_dim = init_blocks(rec, block_dict['ExFv'])
	cc_obj.ExDx, cc_obj.ExDx_starters, cc_obj.ExDx_int_dim, cc_obj.ExDx_dim = init_blocks(rec, block_dict['ExDx'])
	cc_obj.FoFo, cc_obj.FoFo_starters, cc_obj.FoFo_int_dim, cc_obj.FoFo_dim = init_blocks(rec, block_dict['FoFo'])
	cc_obj.FvFv, cc_obj.FvFv_starters, cc_obj.FvFv_int_dim, cc_obj.FvFv_dim = init_blocks(rec, block_dict['FvFv'])
	cc_obj.FoDx, cc_obj.FoDx_starters, cc_obj.FoDx_int_dim, cc_obj.FoDx_dim = init_blocks(rec, block_dict['FoDx'])
	cc_obj.FvDx, cc_obj.FvDx_starters, cc_obj.FvDx_int_dim, cc_obj.FvDx_dim = init_blocks(rec, block_dict['FvDx'])
	cc_obj.FoFv, cc_obj.FoFv_starters, cc_obj.FoFv_int_dim, cc_obj.FoFv_dim = init_blocks(rec, block_dict['FoFv'])
	# cc_obj.DxDx, cc_obj.DxDx_starters, cc_obj.DxDx_int_dim, cc_obj.DxDx_dim = init_blocks(rec, block_dict['DxDx'])

	cc_obj.ExExDx, cc_obj.ExExDx_starters, cc_obj.ExExDx_int_dim, cc_obj.ExExDx_dim = init_blocks(rec, block_dict['ExExDx'])
	cc_obj.ExFoFo, cc_obj.ExFoFo_starters, cc_obj.ExFoFo_int_dim, cc_obj.ExFoFo_dim = init_blocks(rec, block_dict['ExFoFo'])
	cc_obj.ExFoFv, cc_obj.ExFoFv_starters, cc_obj.ExFoFv_int_dim, cc_obj.ExFoFv_dim = init_blocks(rec, block_dict['ExFoFv'])
	cc_obj.ExFvFv, cc_obj.ExFvFv_starters, cc_obj.ExFvFv_int_dim, cc_obj.ExFvFv_dim = init_blocks(rec, block_dict['ExFvFv'])
	# print("ExFvFv size = %d doubles = %f GB" %(cc_obj.ExFvFv_dim, cc_obj.ExFvFv_dim / 1024**3) )
	cc_obj.ExFoDx, cc_obj.ExFoDx_starters, cc_obj.ExFoDx_int_dim, cc_obj.ExFoDx_dim = init_blocks(rec, block_dict['ExFoDx'])
	cc_obj.ExFvDx, cc_obj.ExFvDx_starters, cc_obj.ExFvDx_int_dim, cc_obj.ExFvDx_dim = init_blocks(rec, block_dict['ExFvDx'])
	return cc_obj

def init_X2(cc_obj, rec):
	# REAL USEFUL BLOCKS = ['Ex', 'ExDx', 'ExEx', 'ExExDx', 'ExFo', 'ExFoFo', 'ExFoFv', 'ExFv', 'ExFvFv', 'Fo', 'FoFo', 'FoFv', 'Fv', 'FvFv', 'K']
	cc_obj.K,  cc_obj.K_starters, cc_obj.K_int_dim, cc_obj.K_dim   = init_K()
	cc_obj.Ex, cc_obj.Ex_starters, cc_obj.Ex_int_dim, cc_obj.Ex_dim   = init_blocks(rec, block_dict['Ex'])
	cc_obj.Fo, cc_obj.Fo_starters, cc_obj.Fo_int_dim, cc_obj.Fo_dim   = init_blocks(rec, block_dict['Fo'])
	cc_obj.Fv, cc_obj.Fv_starters, cc_obj.Fv_int_dim, cc_obj.Fv_dim   = init_blocks(rec, block_dict['Fv'])
	# cc_obj.Dx, cc_obj.Dx_starters, cc_obj.Dx_int_dim, cc_obj.Dx_dim   = init_blocks(rec, block_dict['Dx'])
	
	cc_obj.ExEx, cc_obj.ExEx_starters, cc_obj.ExEx_int_dim, cc_obj.ExEx_dim = init_blocks(rec, block_dict['ExEx'])
	cc_obj.ExFo, cc_obj.ExFo_starters, cc_obj.ExFo_int_dim, cc_obj.ExFo_dim = init_blocks(rec, block_dict['ExFo'])
	cc_obj.ExFv, cc_obj.ExFv_starters, cc_obj.ExFv_int_dim, cc_obj.ExFv_dim = init_blocks(rec, block_dict['ExFv'])
	cc_obj.ExDx, cc_obj.ExDx_starters, cc_obj.ExDx_int_dim, cc_obj.ExDx_dim = init_blocks(rec, block_dict['ExDx'])
	cc_obj.FoFo, cc_obj.FoFo_starters, cc_obj.FoFo_int_dim, cc_obj.FoFo_dim = init_blocks(rec, block_dict['FoFo'])
	cc_obj.FvFv, cc_obj.FvFv_starters, cc_obj.FvFv_int_dim, cc_obj.FvFv_dim = init_blocks(rec, block_dict['FvFv'])
	# cc_obj.FoDx
	# cc_obj.FoDx, cc_obj.FoDx_starters, cc_obj.FoDx_int_dim, cc_obj.FoDx_dim = init_blocks(rec, block_dict['FoDx'])
	# cc_obj.FvDx, cc_obj.FvDx_starters, cc_obj.FvDx_int_dim, cc_obj.FvDx_dim = init_blocks(rec, block_dict['FvDx'])
	cc_obj.FoFv, cc_obj.FoFv_starters, cc_obj.FoFv_int_dim, cc_obj.FoFv_dim = init_blocks(rec, block_dict['FoFv'])
	# cc_obj.DxDx, cc_obj.DxDx_starters, cc_obj.DxDx_int_dim, cc_obj.DxDx_dim = init_DxDx(rec, block_dict['DxDx'])
	
	cc_obj.ExExDx, cc_obj.ExExDx_starters, cc_obj.ExExDx_int_dim, cc_obj.ExExDx_dim = init_blocks(rec, block_dict['ExExDx'])
	cc_obj.ExFoFo, cc_obj.ExFoFo_starters, cc_obj.ExFoFo_int_dim, cc_obj.ExFoFo_dim = init_blocks(rec, block_dict['ExFoFo'])
	cc_obj.ExFoFv, cc_obj.ExFoFv_starters, cc_obj.ExFoFv_int_dim, cc_obj.ExFoFv_dim = init_blocks(rec, block_dict['ExFoFv'])
	cc_obj.ExFvFv, cc_obj.ExFvFv_starters, cc_obj.ExFvFv_int_dim, cc_obj.ExFvFv_dim = init_blocks(rec, block_dict['ExFvFv'])
	# cc_obj.ExFoDx
	# cc_obj.ExFvDx
	# cc_obj.ExFoDx, cc_obj.ExFoDx_starters, cc_obj.ExFoDx_int_dim, cc_obj.ExFoDx_dim = init_blocks(rec, block_dict['ExFoDx'])
	# cc_obj.ExFvDx, cc_obj.ExFvDx_starters, cc_obj.ExFvDx_int_dim, cc_obj.ExFvDx_dim = init_blocks(rec, block_dict['ExFvDx'])

	return cc_obj

def init_X3(cc_obj, rec):
	# REAL USEFUL BLOCKS = ['Ex', 'ExEx', 'ExFo', 'ExFv']
	# cc_obj.K,  cc_obj.K_starters, cc_obj.K_int_dim, cc_obj.K_dim   = init_K()
	cc_obj.Ex, cc_obj.Ex_starters, cc_obj.Ex_int_dim, cc_obj.Ex_dim   = init_blocks(rec, block_dict['Ex'])
	# cc_obj.Fo, cc_obj.Fo_starters, cc_obj.Fo_int_dim, cc_obj.Fo_dim   = init_blocks(rec, block_dict['Fo'])
	# cc_obj.Fv, cc_obj.Fv_starters, cc_obj.Fv_int_dim, cc_obj.Fv_dim   = init_blocks(rec, block_dict['Fv'])
	cc_obj.ExEx, cc_obj.ExEx_starters, cc_obj.ExEx_int_dim, cc_obj.ExEx_dim = init_blocks(rec, block_dict['ExEx'])
	cc_obj.ExFo, cc_obj.ExFo_starters, cc_obj.ExFo_int_dim, cc_obj.ExFo_dim = init_blocks(rec, block_dict['ExFo'])
	cc_obj.ExFv, cc_obj.ExFv_starters, cc_obj.ExFv_int_dim, cc_obj.ExFv_dim = init_blocks(rec, block_dict['ExFv'])
	# cc_obj.ExDx, cc_obj.ExDx_starters, cc_obj.ExDx_int_dim, cc_obj.ExDx_dim = init_blocks(rec, block_dict['ExDx'])
	# cc_obj.FvFv, cc_obj.FvFv_starters, cc_obj.FvFv_int_dim, cc_obj.FvFv_dim = init_blocks(rec, block_dict['FvFv'])
	# cc_obj.FoFv, cc_obj.FoFv_starters, cc_obj.FoFv_int_dim, cc_obj.FoFv_dim = init_blocks(rec, block_dict['FoFv'])
	# cc_obj.FoFo, cc_obj.FoFo_starters, cc_obj.FoFo_int_dim, cc_obj.FoFo_dim = init_blocks(rec, block_dict['FoFo'])
	
	# cc_obj.ExExDx, cc_obj.ExExDx_starters, cc_obj.ExExDx_int_dim, cc_obj.ExExDx_dim = init_blocks(rec, block_dict['ExExDx'])
	# cc_obj.ExFoFv, cc_obj.ExFoFv_starters, cc_obj.ExFoFv_int_dim, cc_obj.ExFoFv_dim = init_blocks(rec, block_dict['ExFoFv'])
	# cc_obj.ExFoFo, cc_obj.ExFoFo_starters, cc_obj.ExFoFo_int_dim, cc_obj.ExFoFo_dim = init_blocks(rec, block_dict['ExFoFo'])
	# cc_obj.ExFvFv, cc_obj.ExFvFv_starters, cc_obj.ExFvFv_int_dim, cc_obj.ExFvFv_dim = init_blocks(rec, block_dict['ExFvFv'])
	
	return cc_obj

def init_X4(cc_obj, rec):
	#
	# Only ExEx
	# cc_obj.Ex, cc_obj.Ex_starters, cc_obj.Ex_int_dim, cc_obj.Ex_dim   = init_blocks(rec, block_dict['Ex'])
	cc_obj.ExEx, cc_obj.ExEx_starters, cc_obj.ExEx_int_dim, cc_obj.ExEx_dim = init_blocks(rec, block_dict['ExEx'])
	# cc_obj.ExFo, cc_obj.ExFo_starters, cc_obj.ExFo_int_dim, cc_obj.ExFo_dim = init_blocks(rec, block_dict['ExFo'])
	# cc_obj.ExFv, cc_obj.ExFv_starters, cc_obj.ExFv_int_dim, cc_obj.ExFv_dim = init_blocks(rec, block_dict['ExFv'])

	return cc_obj

def init_Omega_block_zero(cc_obj, rec):
	cc_obj.K,  cc_obj.K_starters, cc_obj.K_int_dim, cc_obj.K_dim    = init_K()
	cc_obj.Ex, cc_obj.Ex_starters, cc_obj.Ex_int_dim, cc_obj.Ex_dim = init_blocks(rec, block_dict['Ex'])
	cc_obj.ExEx, cc_obj.ExEx_starters, cc_obj.ExEx_int_dim, cc_obj.ExEx_dim = init_blocks(rec, block_dict['ExEx'])
	return cc_obj


def scale_cluster_op(cc_obj, factor):
	obj_dict = cc_obj.__dict__
	for key in obj_dict.keys():
		if '_' not in key and obj_dict[key] is not None:  # NOT DIRTY BUT DEPENDS ON NAMINGS...
			obj_dict[key] *= factor
	return cc_obj


def add_Ex_block_ops(destination_op_obj, cc_op_obj):
	# __dict__ are shallow copies
	#
	dest_dict  = destination_op_obj.__dict__
	cc_op_dict = cc_op_obj.__dict__
	#
	for key in ['K', 'Ex', 'ExEx']:
		if cc_op_dict[key] is not None:
			if dest_dict[key] is None:
				dest_dict[key] = deepcopy(cc_op_dict[key])
			else:
				dest_dict[key] += cc_op_dict[key]
	return destination_op_obj


dict_Xn = {'1': init_X1, '2': init_X2, '3': init_X3, '4': init_X4}

def commute(prev_X_obj, prev_T_obj, nth, rec_num_states, resources):
	new_X_obj = dict_Xn[ str(nth) ]( cc_operator(), rec_num_states )
	c_caller.compute_next_X( prev_X_obj, prev_T_obj, new_X_obj, rec_num_states, resources.n_cores )
	return new_X_obj


def print_amps(X_obj):
	print('========================================')
	X_dict = X_obj.__dict__
	for key in sorted(X_dict.keys()):
		if '_' not in key and len(key) < 5:
			print(key)
			print(X_dict[key])
	print('========================================')



def make_omega(H, T, rec_num_states, resources):
	omega = init_Omega_block_zero(cc_operator(), rec_num_states)
	X = H   # X_0 = H
	# print("Copying X_0 to Omega " + '-' * 50)
	omega = add_Ex_block_ops(omega, X)

	for n in [1,2,3,4]:
		# print("Computing [X_%d, T] "  %(n) + '-' * 50 )
		X     = commute(X, T, n, rec_num_states, resources)
		X     = scale_cluster_op(X, 1.0/n)
		omega = add_Ex_block_ops(omega, X)
		
	return omega


class BakerCampbellHausdorff(object):
	""" This class is a wrapper for the nested-operator-specific dirty-work.  It has an interface needed by the coupled_cluster module, but hides the internal nature of the operator.  """
	def __init__(self, rec_num_states, resources):
		self.rec_num_states = rec_num_states
		self.resources      = resources
	
	def computeOmega(self, H, T, textlog):	# call must have this signature!   Computes excitation part (0th, 1st and 2nd order) of HBar
		return make_omega(H, T, self.rec_num_states, self.resources)

	@staticmethod
	def Energy(omega, textlog):		# call must have this signature!   Retrieves energy from excitation part (0th order component) of HBar
		return omega.K[0]

