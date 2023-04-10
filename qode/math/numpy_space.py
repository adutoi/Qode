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


# Is "diagonal" too specific?  How to make things generally expandable?  What about an "empty_vector" function? ... or should that replace
# the use of explicit 0 (and be done at a higher level)?


class _generic(object):
    """ Since numpy mostly treats complex and real on the same footing, this code is general; for example .conj() just make a copy of real arrays """
    def __init__(self,dim):
        self.dim   = dim
    def check_member(self,v):
        if not isinstance(v,numpy.ndarray):   raise error("numpy array expected")
        if len(v.shape)!=1:                   raise error("1D array expected")
        if v.shape[0]!=self.dim:              raise error("vector has wrong dimension")
        if type(v[0]) is not self.field:      raise error("vector contains scalars of wrong field type")	# Just a "spot check"? ... numpy insists array is homogeneous?
    def check_lin_op(self,op):
        if not isinstance(op,numpy.ndarray):  raise error("numpy array expected")
        if len(op.shape)==2:
            if op.shape[0]!=op.shape[1]:         raise error("square matrix expected")
            if op.shape[0]!=self.dim:            raise error("matrix has wrong dimension")
            if type(op[0,0]) is not self.field:  raise error("matrix contains scalars of wrong field type")
            functionable = False
        elif len(op.shape)==1:
            if op.shape[0]!=self.dim:            raise error("diagonal array has wrong dimension")
            if type(op[0]) is not self.field:    raise error("diagonal array contains scalars of wrong field type")
            functionable = True
        else:  raise error("1D or 2D array expected")
        return functionable
    @staticmethod
    def dot(v,w):
        return numpy.dot(v.conj(),w)
    @staticmethod
    def add_to(v,w,c=1):
        if c==1:  v +=   w
        else:     v += c*w
    @staticmethod
    def scale(n,v):
        v *= n
    @staticmethod
    def copy(v):
        return v.copy()
    @staticmethod
    def function_on_diags(func,op):
        func = numpy.vectorize(func)
        return func(op)
    @staticmethod
    def act_on_vec(op,v):
        if len(op.shape)==2:  return numpy.dot(op,v)
        else:                 return numpy.multiply(op,v)
    @staticmethod
    def act_on_vec_block(op,v_block):		# Part of an incomplete progression towards also allowing blocks of vectors (no efficiency gains here, but other traits classes might be smarter)
        if len(op.shape)==2:  return [ numpy.dot(op,v)      for v in v_block ]
        else:                 return [ numpy.multiply(op,v) for v in v_block ]
    @staticmethod
    def back_act_on_vec(v,op):
        if   len(op.shape)==2:  return numpy.dot(v.conj(),op).conj()
        else:                   return numpy.multiply(op.conj(),v)
    @staticmethod
    def diagonal(op):
        if   len(op.shape)==2:
            diags = numpy.array([self.field(0)]*self.dim)
            for i in range(self.dim):  diags[i] = op[i,i]
            return diags
        else: return self	# am I being lazy, or is it really ok not to make a copy, since any function on diags makes a copy

class real_traits(_generic):
    def __init__(self,dim):
        _generic.__init__(self,dim)
        self.field = numpy.float64

class cplx_traits(_generic):
    def __init__(self,dim):
        _generic.__init__(self,dim)
        self.field = numpy.complex128



