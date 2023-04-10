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
# import ctypes as ct
import numpy  as np
from io import StringIO

class cc_operator(object):
	"""operator class holding H, X, T, or Omega components"""
	def __init__(self):
		self.K      = None
		self.Ex     = None
		self.Fo     = None
		self.Fv     = None
		self.Dx     = None
		self.ExEx   = None   # diag = zero
		self.ExFo   = None   # diag = zero
		self.ExFv   = None   # diag = zero
		self.ExDx   = None
		self.FoFo   = None
		self.FvFv   = None
		self.FoFv   = None   # diag = zero
		# self.FvFo   = None   # diag = zero --> removed in the latest version (Aug 29, 2016)
		self.FoDx   = None   # diag = zero
		self.FvDx   = None   # diag = zero
		self.DxDx   = None
		self.ExFoFo = None
		self.ExFvFv = None
		self.ExFoFv = None
		# self.ExFvFo = None # --> removed in the latest version (Aug 29, 2016)
		self.ExFoDx = None
		self.ExFvDx = None
		self.ExExDx = None
	def __str__(self):
		out = StringIO()
		names = (    "K",    "Ex",    "Fo",    "Fv",    "Dx",    "ExEx",    "ExFo",    "ExFv",    "ExDx",    "FoFo",    "FvFv",    "FoFv",    "FoDx",    "FvDx",    "DxDx",    "ExFoFo",    "ExFvFv",    "ExFoFv",    "ExFoDx",    "ExFvDx",    "ExExDx")
		vals  = (self.K, self.Ex, self.Fo, self.Fv, self.Dx, self.ExEx, self.ExFo, self.ExFv, self.ExDx, self.FoFo, self.FvFv, self.FoFv, self.FoDx, self.FvDx, self.DxDx, self.ExFoFo, self.ExFvFv, self.ExFoFv, self.ExFoDx, self.ExFvDx, self.ExExDx)
		for n,v in zip(names,vals):  print(n+"\n",v,file=out)
		#print_contents(out, "K", self.K)
		return out.getvalue()


# def print_contents(list_double_star): 
# 	try:
# 		print('---------- List of double* Pointers ------------------')
# 		for ptr_obj in list_double_star:
# 			print('double* pointer at',  hex(ct.addressof(ptr_obj.contents)) , '   w/ value =', ptr_obj.contents )
# 		print('------------------------------------------------------')
# 	except:
# 		print('------------------------------------------------------')


def init_K(): # ONE MOLECULAR INDEX
	# print("K     block has 1 double precision numbers.")
	return np.zeros(1), np.zeros(1, dtype=np.int64), 1, 1

def init_blocks(rec_num_states, trasitions):
	# transitions = ['OCC','VRT'] / ['OCC','VRT','OCC','VRT'], ['OCC','VRT', 'OCC','VRT','OCC','VRT']
	#
	Nmol = len(rec_num_states)
	double_loc = [0]
	Nv = [ num - 1 for num in rec_num_states ]
	if len(trasitions) == 2:
		# One Operator Block
		for n in range(Nmol):
			if trasitions[0] == 'OCC':
				size = 1
			elif trasitions[0] == 'VRT':
				size = Nv[n]
			else:
				raise AssertionError
			if trasitions[1] == 'OCC':
				pass
			elif trasitions[1] == 'VRT':
				size *= Nv[n]
			else:
				raise AssertionError
			# print("size =", size)
			double_loc += [ size + double_loc[-1]]
		# print("%d double numbers" %(double_loc[-1]))
		starter_dim = len(double_loc) - 1
		starters = np.ascontiguousarray( np.zeros(starter_dim, dtype=np.int64), dtype=np.int64 )
		for i in range(starter_dim):
			starters[i] = double_loc[i]
		return  np.ascontiguousarray( np.zeros(double_loc[-1]), dtype=np.float64 ),  starters, starter_dim, double_loc[-1]

	elif len(trasitions) == 4:
		# Two Operator Block
		for m in range(Nmol):
			if trasitions[0] == 'OCC':
				size1 = 1
			elif trasitions[0] == 'VRT':
				size1 = Nv[m]
			else:
				raise AssertionError
			if trasitions[1] == 'OCC':
				pass
			elif trasitions[1] == 'VRT':
				size1 *= Nv[m]
			else:
				raise AssertionError
			#
			for n in range(Nmol):
				if trasitions[2] == 'OCC':
					size2 = 1
				elif trasitions[2] == 'VRT':
					size2 = Nv[n]
				else:
					raise AssertionError
				if trasitions[3] == 'OCC':
					pass
				elif trasitions[3] == 'VRT':
					size2 *= Nv[n]
				else:
					raise AssertionError
				size = size1 * size2
				# print("size =", size)
				double_loc += [ size + double_loc[-1]]	
		# print("%d double numbers" %(double_loc[-1]))
		starter_dim = len(double_loc) - 1
		starters =  np.ascontiguousarray( np.zeros(starter_dim, dtype=np.int64), dtype=np.int64 )
		for i in range(starter_dim):
			starters[i] = double_loc[i]
		return  np.ascontiguousarray( np.zeros(double_loc[-1]), dtype=np.float64 ), starters, starter_dim, double_loc[-1]
		
	elif len(trasitions) == 6:
		# Three Operator Block
		for m in range(Nmol):
			if trasitions[0] == 'OCC':
				size1 = 1
			elif trasitions[0] == 'VRT':
				size1 = Nv[m]
			else:
				raise AssertionError
			if trasitions[1] == 'OCC':
				pass
			elif trasitions[1] == 'VRT':
				size1 *= Nv[m]
			else:
				raise AssertionError
			#
			for n in range(Nmol):
				if trasitions[2] == 'OCC':
					size2 = 1
				elif trasitions[2] == 'VRT':
					size2 = Nv[n]
				else:
					raise AssertionError
				if trasitions[3] == 'OCC':
					pass
				elif trasitions[3] == 'VRT':
					size2 *= Nv[n]
				else:
					raise AssertionError
				#
				for k in range(Nmol):
					if trasitions[4] == 'OCC':
						size3 = 1
					elif trasitions[4] == 'VRT':
						size3 = Nv[n]
					else:
						raise AssertionError
					if trasitions[5] == 'OCC':
						pass
					elif trasitions[5] == 'VRT':
						size3 *= Nv[n]
					else:
						raise AssertionError		
					size = size1 * size2 * size3
					# print("size =", size)
					double_loc += [ size + double_loc[-1]]	
		# print("%d double numbers" %(double_loc[-1]))
		starter_dim = len(double_loc) - 1
		starters =  np.ascontiguousarray( np.zeros(starter_dim, dtype=np.int64), dtype=np.int64 )
		for i in range(starter_dim):
			starters[i] = double_loc[i]
		return  np.ascontiguousarray( np.zeros(double_loc[-1]), dtype=np.float64 ), starters, starter_dim, double_loc[-1]

	else:
		raise AssertionError



def T_operator(rec):
	cc_obj = cc_operator()
	cc_obj.K,  cc_obj.K_starters, cc_obj.K_int_dim, cc_obj.K_dim    = init_K()
	cc_obj.Ex, cc_obj.Ex_starters, cc_obj.Ex_int_dim, cc_obj.Ex_dim = init_blocks(rec, ['OCC','VRT'])
	cc_obj.ExEx, cc_obj.ExEx_starters, cc_obj.ExEx_int_dim, cc_obj.ExEx_dim = init_blocks(rec, ['OCC','VRT','OCC','VRT'])


	return cc_obj


def H_operator(rec):
	cc_obj = cc_operator()

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
	cc_obj.DxDx, cc_obj.DxDx_starters, cc_obj.DxDx_int_dim, cc_obj.DxDx_dim = init_blocks(rec, block_dict['DxDx'])
	return cc_obj
