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
from .integrals import overlap, orbital_operator

"""\
This is good enough for now, but it is actually an abuse of the space class.  This notation should also be implemented directly for basis functions.
See note in integrals.py.
"""

class traits(object):
    """ Since numpy mostly treats complex and real on the same footing, this code is general; for example .conj() just make a copy of real arrays """
    def __init__(self,basis):
        self.field = numpy.float64
        self.basis = basis
        self.dim   = len(basis)
        self.olap  = numpy.zeros((self.dim,self.dim))
        for m,Bm in enumerate(basis):
            for n,Bn in enumerate(basis):
                self.olap[m,n] = overlap(Bm,Bn)
        self.olap_inv = numpy.linalg.inv(self.olap)
    def check_member(self,v):
        if not isinstance(v,numpy.ndarray):   raise Exception("numpy array expected")
        if len(v.shape)!=1:                   raise Exception("1D array expected")
        if v.shape[0]!=self.dim:              raise Exception("vector has wrong dimension")
        if type(v[0]) is not self.field:      raise Exception("vector contains scalars of wrong field type")	# Just a "spot check"? ... numpy insists array is homogeneous?
    def check_lin_op(self,op):
        if not isinstance(op,orbital_operator):  raise Exception("class orbital_operator expected")
        return False
    def dot(self,v,w):
        return numpy.dot(v, numpy.dot(self.olap,w))
    def add_to(self,v,w,c=1):
        if c==1:  v +=   w
        else:     v += c*w
    def scale(self,n,v):
        v *= n
    def copy(self,v):
        return v.copy()
    def act_on_vec(self,op,v):
        new = numpy.zeros(self.dim)
        for m,Bm in enumerate(self.basis):
            for n,Bn in enumerate(self.basis):
                new[m] += op(Bm,Bn) * v[n]
        return numpy.dot(self.olap_inv, new)
    def back_act_on_vec(self,v,op):
        return self.act_on_vec(op,v)
