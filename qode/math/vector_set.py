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
from collections.abc import MutableSequence
import numpy
from . import field_traits
from .space import abs, conjugate, is_int_zero
from .space import _member, _operator_base	# for type-checking only



# To do list
#
# 1. This is actually a prototype for the concept of using syntax like (set|opr|set) to return a 2D list of values (castable to a matrix), but none of 
#    that has yet been implemented. ... actually that syntax is more general. there should be a generic kind of rectilinear list (a tensor) that
#    understands various bilinear forms and passes them on to the elements during either outer product or contraction taking of the exterior indices.
# 2. Seems weird to have space as a constructor argument.  It should it be the space that is in charge of issuing vector sets (using a member function).
# 3. At one point, I thought that I should let the user provide some back-end storage for blocks of vectors so that they can be contiguous; however,
#    It turns out to be kind of a nightmare to figure out how to let a user enforce exactly where and how they want the vectors stored when some
#    algorithms generate new vectors one at a time that are used to build vector_sets.  So the conclusion is that this is the correct level to be 
#    implementing this at.  If the user wants to get really picky, they can implement a vector class that is really just a shell, pointing to physical
#    information in some kind of more sophisticated container.  At least then, the user could enforce things like having vectors that are generated
#    sequentially in time being stored next to one another and could then test that this is true before any operations that rely on that.  I think
#    that is about as good as it gets.  It does seem like I want to embed this a bit deeper into space though to get the notation mentioned in point 1
#    to work.  As it is now, for example in lanczos, I have to make a direct call to the space.traits member functions to perform the action on a block.
# 4. higher level tasks, like computing orthonormal bases or managing biorthogonal complements (in general, already implied a bit here)
#    should be left to the sub_basis class (see prototype in tensor.py)


class vector_set(MutableSequence):
    """\
    This mostly just acts like a list of vectors with only a couple extra bells and whistles:
    It knows the scalar field type over which the vector space is defined, accessible as the member field.
    """
    def __init__(self,space,vector_list=None,orthonormal=False):	# orthonormal=True is an unchecked assertion, specifically to prevent superfluous computation of the overlap matrix
        if vector_list is None:  vector_list = []	# mutable default arguments are dangerous.  See:  http://python-guide-pt-br.readthedocs.io/en/latest/writing/gotchas/
        self._vector_list = vector_list
        self.space = space
        self.field = self.space.field
        for vector in vector_list:
            if not is_int_zero(vector) and vector.space is not self.space:  raise Exception('Vectors in a set are expected to be from the same space.')
        self.orthonormal = orthonormal
        self._S_cache    = (None, None)
    # These are the functions needed to emulate a list
    def __getitem__(self,index):
        if isinstance(index,slice):  return vector_set(self.space,self._vector_list[index])
        else:                        return                       self._vector_list[index]
    def __setitem__(self,index,value):
        if not is_int_zero(value) and value.space is not self.space:  raise Exception('Vectors in a set are expected to be from the same space.')
        self._S_cache = (None, None)
        self._vector_list[index] = value
        return value
    def __delitem__(self,index):
        self._S_cache = (None, None)
        del self._vector_list[index]
    def __len__(self):
        return len(self._vector_list)
    def insert(self,index,value):
        if not is_int_zero(value) and value.space is not self.space:  raise Exception('Vectors in a set are expected to be from the same space.')
        self._S_cache = (None, None)
        return self._vector_list.insert(index,value)
    # These are the "extra" functions
    def overlaps(self):
        if self.orthonormal:  raise LogicError("Set was explicitly asserted to be orthonormal.")
        S, Sinv = self._S_cache
        if S is None:
            S = numpy.array(self.space.dot_vec_blocks(self,self))
            for i in range(len(self)):  S[i,i] = field_traits.metric[self.field](S[i,i])      # Should I check before discarding what should be only a round-off error in the cast?
            self._S_cache = (S, None)
        return S
    def metric_matrix(self):
        S, Sinv = self._S_cache
        if S    is None:  S = self.overlaps()
        if Sinv is None:
            try:
                Sinv = numpy.linalg.inv(S)
            except LinAlgError:
                print("Linear dependency in vector set.  Projection is ambiguous.")
                raise
            self._S_cache = (S, Sinv)
        return Sinv
    def orthonormality(self):
        one = self.field(1)
        dim = len(self)
        S = self.overlaps()
        biggest = field_traits.metric[self.field](0)
        for i in range(dim):
            if abs(one-S[i,i])>biggest:  biggest = abs(one-S[i,i])
            for j in range(i+1,dim):
                if abs(S[i,j])>biggest:  biggest = abs(S[i,j])
        return biggest
    def is_orthonormal(self,thresh=1e-10):  return self.orthonormality() < field_traits.metric[self.field](thresh)
    def projections(self,obj):
        # obj might be a vector or it might be an operator.  This returns a list of raw projections onto the members of the vector set.
        # I think eventually this should also be called via the syntax vector_set|object, so that (vector_set|operator|vector_set) gives a matrix-valued (2D list) return
        # should be identical to single line of code:
        #   return [ v|obj for v in self ]
        if   isinstance(obj,_member):
            return [ prj[0] for prj in self.space.dot_vec_blocks(self,[obj]) ]		# w/o prj[0], it comes out as list of list with single elements
        elif isinstance(obj,_operator_base):
            return obj.back_act_on_vec_block(self)
        else:  
            return [ v|obj for v in self ]	# handles case where obj is a scalar ... maybe someday we'll want block scale too?
    def project(self,obj):
        # When executed on a vector, this returns a resolution of that vector as a linear combination of the vectors in this set.
        # When executed on an operator, the result is a list of vectors, the projection of any vector onto which is identical to the resolution of the operator acting on that vector.
        if len(self)>0:
            projections = self.projections(obj)
            if self.orthonormal:
                return projections
            else:
                Sinv = self.metric_matrix()
                return numpy.array([ sum(Sinv_ij*a_j for Sinv_ij,a_j in zip(Sinv_ix,projections)) for Sinv_ix in Sinv.tolist() ])
        else:  return []
    def deproject(self,projection):
        # If the projection is a list of coefficients, then this builds a vector as a linear combination of the vectors in this set.
        # If the projection is a list of vectors (see project() above), then a sum of dyads is returned, suitable for use as an operator projected into a subspace.
        if len(projection)!=len(self):  raise "projection must have same length as vector set for deprojection"
        deprojection = 0
        for vec_i,proj_i in zip(self,projection):
            if proj_i!=0:  deprojection += vec_i * proj_i      # test only affects efficiency, so ok if not stable ... really important for resolving basis vectors (with explicitly set zeros)
        return deprojection
    def projection_opr(self):
        # Projecting and deprojecting the unit operator gives the operator that corresponds to a projection into the space spanned by this set.
        # Right now, this does not work for linearly dependent set, but the outcome is well-defined . . . fix that?
        # There is some danger here in that this makes an operator whose provenance might be forgotten.  Then if you modify this set, that operator is implicitly changed.
        return self.deproject(self.project(1))
