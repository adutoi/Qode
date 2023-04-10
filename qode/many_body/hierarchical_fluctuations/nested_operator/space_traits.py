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



class space_traits_class(object):
	def __init__(self):
		self.field = numpy.float64
	@staticmethod
	def dot(v,w):
		return v.dot(w)
	@staticmethod
	def add_to(v,w,n=1):
		v.increment(w,n)
	@staticmethod
	def scale(n,v):
		v.scale(n)
	@staticmethod
	def copy(v):
		return v.copy()
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



# for now, this acts like a singleton, but, eventually, I will want to put in some kind of dimension checking, in which case this will be deleted, the above will have "_class" removed, and that will be used

space_traits = space_traits_class()
