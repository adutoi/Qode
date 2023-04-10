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
from . import  space
from          .space import conjugate

class vector(space.vector):
    def __init__(self,vec):
        if not isinstance(vec,numpy.matrix):  raise error("numpy matrix expected")
        if vec.shape[1]!=1:  raise error("column matrix expected")
        self.dim = vec.shape[0]
        self.field = type(vec[0,0])
        self.vec = vec
    @staticmethod
    def vecs_are_compatible(v,w):
        if v.dim==w.dim:  return True
        else:             return False
    @staticmethod
    def dot(v,w):  return (v.vec.T*w.vec)[0,0]
    @staticmethod
    def sum(v,w):  return vector(v.vec+w.vec)
    @staticmethod
    def scaled(c,v):  return vector(c*v.vec)

class diagonal_operator(space.functionable_op):
    def __init__(self,diags):
        if not isinstance(diags,numpy.ndarray):  raise error("numpy array expected")
        if len(diags.shape)!=1:  raise error("1D array expected")
        self.dim = diags.shape[0]
        self.field = type(diags[0])
        self.diags = diags
    def function_on_diags(self,func):
        func = numpy.vectorize(func)
        self.diags = func(self.diags)
    def copy(self):  return diagonal_operator(self.diags)
    @staticmethod
    def ops_are_compatible(A,B):
        if A.dim==B.dim:  return True
        else:             return False
    @staticmethod
    def op_and_vec_are_compatible(A,v):
        if A.dim==v.dim:  return True
        else:             return False
    @staticmethod
    def act_on_vec(A,v):
        dim = A.dim
        vec = numpy.array([A.field(0)]*dim)
        for i in range(dim):  vec[i] = v.vec[i,0]
        new_vecA = numpy.multiply(A.diags,vec)
        new_vecM = numpy.matrix([[A.field(0)]]*dim)
        for i in range(dim):  new_vecM[i,0] = new_vecA[i]
        return vector(new_vecM)
    @staticmethod
    def back_act_on_vec(v,A):  return conjugate(A).act_on_vec(A,v)

class operator(space.operator):
    def __init__(self,mat):
        if not isinstance(mat,numpy.matrix):  raise error("numpy matrix expected")
        if mat.shape[0]!=mat.shape[1]:  raise error("square matrix expected")
        self.dim = mat.shape[0]
        self.field = type(mat[0,0])
        self.mat = mat
    def diagonal(self):
        dim = self.dim
        diags = numpy.array([self.field(0)]*dim)
        for i in range(dim):  diags[i] = self.mat[i,i]
        return diagonal_operator(diags)
    @staticmethod
    def ops_are_compatible(A,B):
        if A.dim==B.dim:  return True
        else:             return False
    @staticmethod
    def op_and_vec_are_compatible(A,v):
        if A.dim==v.dim:  return True
        else:             return False
    @staticmethod
    def act_on_vec(A,v):  return vector(A.mat*v.vec)
    @staticmethod
    def back_act_on_vec(v,A):  return  vector((v.vec.H*A.mat).H)
