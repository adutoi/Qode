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
# import numpy  as np
import ctypes as ct

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



def print_contents(list_double_star):
	try:
		print('---------- List of double* Pointers ------------------')
		for ptr_obj in list_double_star:
			print(ptr_obj.contents,  hex(ct.addressof(ptr_obj.contents)) )
		print('------------------------------------------------------')
	except:
		print('------------------------------------------------------')



