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



class nested_operator(object):
	def __init__(self, num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt):
		self._bounds = num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt	# immutable
		self._Ex_amps   = None
		self._Fv_amps   = None
		self._Fo_amps   = None
		self._Dx_amps   = None
		self._Ex_subops = None
		self._Fv_subops = None
		self._Fo_subops = None
		self._Dx_subops = None

	def __iter__(self):
	def __call__(self,the_state):
	def scale(self,c):
	def dot(self,other):
	def increment(self,other,c=1):
	def add_term(self,term):
	def right_multiply(self, transition):
	def left_multiply(self, transition, shallow=False):
		# Start here and fill out the rest as needed
	def _subtarget(self, transition):



	def _get_Ex_amps(self):
		if self._Ex_amps is None:
			num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt = self._bounds
			self._Ex_amps = numpy.array( [[0]*num_OccDst]*num_VrtCrt )
		return self._Ex_amps
	def _get_Fv_amps(self):
		if self._Fv_amps is None:
			num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt = self._bounds
			self._Fv_amps = numpy.array( [[0]*num_VrtDst]*num_VrtCrt )
		return self._Fv_amps
	def _get_Fo_amps(self):
		if self._Fo_amps is None:
			num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt = self._bounds
			self._Fo_amps = numpy.array( [[0]*num_OccCrt]*num_OccDst )
		return self._Fo_amps
	def _get_Dx_amps(self):
		if self._Dx_amps is None:
			num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt = self._bounds
			self._Dx_amps = numpy.array( [[0]*num_OccCrt]*num_VrtDst )
		return self._Dx_amps



	def _get_Ex_subops(self):
		if self._Ex_subops is None:
			num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt = self._bounds
			self._Ex_subops = [[nested_operator(num_VrtCrt-1, num_VrtDst, num_OccDst-1, num_OccCrt) for i in range(num_OccDst)] for a in range(num_VrtCrt)]
		return self._Ex_subops
	def _get_Fv_subops(self):
		if self._Fv_subops is None:
			num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt = self._bounds
			self._Fv_subops = [[nested_operator(num_VrtCrt-1, num_VrtDst-1,      0,     num_OccCrt) for b in range(num_VrtDst)] for a in range(num_VrtCrt)]
		return self._Fv_subops
	def _get_Fo_subops(self):
		if self._Fo_subops is None:
			num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt = self._bounds
			self._Fo_subops = [[nested_operator(      0,    num_VrtDst, num_OccDst-1, num_OccCrt-1) for j in range(num_OccCrt)] for i in range(num_OccDst)]
		return self._Fo_subops
	def _get_Dx_subops(self):
		if self._Dx_subops is None:
			num_VrtCrt, num_VrtDst, num_OccDst, num_OccCrt = self._bounds
			self._Dx_subops = [[nested_operator(      0,    num_VrtDst-1,      0,     num_OccCrt-1) for i in range(num_OccCrt)] for a in range(num_VrtDst)]
		return self._Dx_subops
