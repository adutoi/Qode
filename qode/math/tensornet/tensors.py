#    (C) Copyright 2023, 2025 Anthony D. Dutoi
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
from textwrap import indent
from qode.util import struct
from .base import resolve, evaluate, raw, scalar_value, increment, resolve_ellipsis
from . import network_logic

_backend_contract_path = False    # if True, let backend handle finding the optimal contraction path upon evaluate() call
def backend_contract_path(TrueFalse):
    global _backend_contract_path
    _backend_contract_path = TrueFalse



# We have 4 different types of top-level tensors
#   contraction_expression
#   tensor_sum
#   tensor_network
#   primitive_tensor
# All of these are built off of tensor_base.  The latter three are built off of resolved_tensor (which
# is based on tensor_base) indicating that they are completely well-formed.  contraction_expression
# behaves like a tensor so long as it is well formed, but it could be incomplete.  The basic problem
# is that, when connecting a bunch of tensors with @, it is not possible to automatically discern
# when the string is finished (ie, when all intended target free indices are present).  Only the user can decide
# when to "close" the string, but we don't want to make the user do this explicitly.  Therefore,
# we allow this incomplete class, whose *resolution* is a tensor_network (it always resolves to a
# tensor_network). This is triggered by taking another action on it that is not @ or *.  A tensor_sum
# can only contain types primitive_tensor and tensor_network.  An attempt to include another
# tensor_sum will cause it to be expanded out as a single sum (and contraction_expression will 
# be automatically resolved to tensor_network).  A tensor_network can only contain
# type primitive_tensor, and including one network (or contraction_expression) into another network causes explicit expansion
# into a single network.  Taking this together, this means that a tensor_network cannot contain a 
# tensor_sum.  An attempt to put a tensor_sum into a tensor_network will result in the hierarchy
# being swapped to a sum of networks.  The reason for this is purely practical.  At the moment,
# I cannot think of a way of finding the best path to evaluate a network that contains sums.



# Expected to have a backend and a shape (if requested), and be able to participate in tensor arithmetic.
# Everything else is implementation specific.
class tensor_base(object):
    #
    # These functions take care of the adminstration of the base-class data and do not demand anything of the child classes
    #
    def __init__(self, backend, shape=None):
        self._backend = backend
        self._shape = shape    # initialization can be delayed if determined by as-of-yet unknown constituent terms or factors
    def _assert_same_backend(self, other):
        if self._backend is not other._backend:
            raise TypeError(f"attempted operation between tensors with different backends: {self._backend.name} and {other._backend.name}")
    def _attempt_set_shape(self, shape):
        if self._shape is None:    # ok if this is the first ...
            self._shape = shape    # ... or only attempt
        if self._shape!=shape:
            raise ValueError(f"failed to set shape to {shape} because previously set to {self._shape}")
    def __getattr__(self, shape):
        # __getattr__ is only called for non hard-coded attributes.  Here I want it only for .shape, which is *defined*
        # for all tensors and *valid* for all properly initialized tensors.  Doing it like this prevents our only 
        # public attribute from being set from the outside (I think?).
        if shape=="shape":
            return self._shape
        else:
            raise AttributeError(f"{type(self)} object has no attribute {shape}")   # By now "shape" is a misnomer, hence the exception
    #
    # As long as __mul__ is defined for the child classes, these methods can be implemented purely abstractly, without needing
    # to first resolve any contraction_expression first (because its product form can contain a scalar).  We also let __imul__ and
    # __itruediv__ default to the python implicit defintion in terms of these (eg, x *= y  ->  x = x * y).
    #
    def __rmul__(self, x):    # needed for leading scalars
        return self * x
    def __truediv__(self, x):
        return self * (1./x)
    def __neg__(self):
        return self * -1
    #
    # In addition to __neg__ (and therefore __mul__) above, these additionally require that any contraction_expression first be resolved,
    # but not evaluated.  The resolve() function (via .resolve()) is implemented specifically for contraction_expression and is just
    # a pass-through for the others. Along with (implicit) __iadd__ and __isub__, this completes all the basic mathematical operations outside
    # of contraction, all of which are lazy in terms of evaluation.
    #
    def __sub__(self, other):
        return self + (-other)
    def __add__(self, other):
        return tensor_sum(self._backend, terms=[resolve(self), resolve(other)])    # shape implied by terms (which must match)
    #
    # The final layer of (lazy) mathematics is contraction.  When __call__ is applied to a tensor, it prepares that tensor for
    # contraction using @, which is defined only for contraction_expression because one *must* first specify free and contraction
    # indices in order to contract.  If no contraction indices are specified, the free indices could just specify permutation.
    # Contracting a contraction_expression with another tensor first triggers resolution of any existing existing expression
    # into a closed and error-checked tensor_network (which has no memory of the "letters" used to specify the internal contractions).
    # Also, the otherwise useless act of using __call__ with an empty argument string (otherwise only valid for scalar tensors)
    # simply forces resolution to a error-checked tensor_network.  This is *logically* equivalent (though the return type is different)
    # to giving the trivial (non)permutation of the correct length, which is *mathematically* equivalent to doing nothing (assuming the
    # contractions string is, in fact, valid.
    #
    def __call__(self, *indices):
        if indices==():
            return resolve(self)
        else:
            return contraction_expression(self._backend, [(resolve(self), indices)])
    #
    # These are the non-mathematical operations (rather informational) which also require any contraction_expression expression to
    # be resolved.  After resolution, the action is implemented separately for each child of resolved_tensor, which overrides __getitem__
    # (ie, the default implementation is only for contraction_expression, but this documents also that it is an expected tensor behavior).
    # __setitem__ is never overridden.  If the indices specify a specific element, it could also then result in evaluation.
    #
    def __getitem__(self, indices):    # called only by contraction_expression (overridden by resolved_tensor types)
        indices = resolve_ellipsis(indices, len(self._shape))
        return resolve(self)[indices]
    def __setitem__(self, indices):
        raise TypeError("slices and/or elements of tensornet tensors are not assignable")
    #
    # The methods below are the tensornet-specific stacke manipulations (ie not the usual mathematical or informational operators expected
    # for all tensors), which are accessible by unbound functions.  The function _resolve has a different implementation for
    # both contraction_expression and the children of resolved_tensor and is implemented only in the child classes.
    # The function evaluate() automatically triggers resolve() for contraction_expression, but after that, it is handled by the child classes
    # of resolved_tensor that override it separately (ie, this implementation is only called by contraction_expression). In turn,
    # raw(), scalar_value(), and increment() each trigger evaluate(), which must result in a primitive_tensor, which further handles
    # the requested task on the evaluated data via overridden functions.  Therefore, tensor_sum and tensor_network need not implement these
    # three functions.
    #
    def _resolve(self):     # overridden by all child classes
        raise NotImplementedError
    def _evaluate(self):    # called only by contraction_expression (overridden by resolved_tensor types)
        return evaluate(resolve(self))
    def _raw(self):                              # All of these are overridden by primitive_tensor (but only primitive_tensor).
        #return raw(evaluate(self))              # <- Better code aesthetics, but raw() on primitive_tensor result of evaluate() forces a copy.
        return evaluate(self)._raw_tensor        # <- More efficient, since we know that evaluate() called on a tensor_sum or tensor_network ...
    def _scalar_value(self):                     #    ... returns freshly generated data that would be transient (never bound to a variable) ...
        return scalar_value(evaluate(self))      #    ... if immediately copied via raw().  Since the ._scalar member equals 1, .raw_tensor ...
    def _increment(self, raw_result):            #    ... is all that is needed, and copying is unnecessary.
        increment(raw_result, evaluate(self))    #



class contraction_expression(tensor_base):
    def __init__(self, backend, tensors_indices, _scalar=1):
        tensor_base.__init__(self, backend)
        self._scalar = _scalar
        self._tensors_indices = tensors_indices
        shape = network_logic.logic.get_shape(tensors_indices, self._assert_same_backend)    # checks backend and contraction lengths
        self._attempt_set_shape(tuple(shape))
    def __mul__(self, x):
        try:
            new_scalar = x * self._scalar
        except:
            raise TypeError(f"* not defined between contraction_expression and {type(scalar)}.\nUse @ or @= to join tensors via contraction or outer product.")
        return contraction_expression(self._backend, self._tensors_indices, new_scalar)
    def __matmul__(self, other):
        try:
            new_tensors_indices = self._tensors_indices + other._tensors_indices
        except:
            raise TypeError(f"@ not defined between contraction_expression and {type(other)}.\nUse * or *= for multiplication by a scalar.")
        return contraction_expression(self._backend, new_tensors_indices, self._scalar)
    def _resolve(self):
        return network_logic.logic.resolve(self._scalar, self._tensors_indices, tensor_network)



class resolved_tensor(tensor_base):
    # __init__ is inherited
    def _resolve(self):
        return self



# The tensornet type for the primitive tensors that the user sees and uses and builds networks from.
# The tensor importantly knows its backend, via a provided module (implemented by the user if not
# already provided for that backend type).
# Only this class uses _scalar.  Use * or *= from outside the class.
class primitive_tensor(resolved_tensor):
    def __init__(self, backend, raw_tensor, _scalar=1):
        resolved_tensor.__init__(self, backend, backend.shape(raw_tensor))
        self._raw_tensor = raw_tensor
        self._scalar     = _scalar    # so that we can define * without changing raw tensor data
    def __mul__(self, x):
        new = primitive_tensor(self._backend, self._raw_tensor, _scalar=x*self._scalar)
        return new
    def __getitem__(self, indices):
        indexed_tensor = self._backend.subscript(self._raw_tensor, indices)
        if any(isinstance(index,slice) for index in indices):
            new = primitive_tensor(self._backend, indexed_tensor, _scalar=self._scalar)
        else:
            new = self._scalar * indexed_tensor    # this is a scalar if we get here
        return new
    def __str__(self):
        data = indent(self._backend.str(self._raw_tensor), "    ")
        return f"tensornet.primitive_tensor(backend = {self._backend.name}\n{self._scalar} *\n{data}\n)"
    def _evaluate(self):
        return self
    def _raw(self):
        return self._backend.mult(self._scalar, self._raw_tensor)    # force copy even if self._scalar==1 in case value modified (by user or tensor_sum)
    def _scalar_value(self):
        if len(self.shape)>0:
            raise ValueError("cannot take the scalar value of a tensornet tensor with >0 free indices")
        return self._scalar * self._backend.scalar_value(self._raw_tensor)
    def _increment(self, raw_result):
        if self._scalar!=1:
            self._backend.increment(raw_result, self._backend.mult(self._scalar, self._raw_tensor))
        else:
            self._backend.increment(raw_result, self._raw_tensor)    # avoid making unnecessary copy just to use as increment



# in an informational sense, this is the most basic tensor, since all it does is store a list
# of tensor_networks and primitive_tensors.
class tensor_sum(resolved_tensor):
    def __init__(self, backend, terms=None, shape=None):
        resolved_tensor.__init__(self, backend)    # assume empty instaniated as empty/blank ...
        self._attempt_set_shape(shape)             # ... (though can specify shape, and ok if again set to None) ...
        self._terms = []                           # ... as accumulator ...
        if terms is not None:                      # ... unless some terms are given
            for term in terms:
                self._assert_same_backend(term)        # <- otherwise raises exception
                self._attempt_set_shape(term.shape)    # <- exception if mismatch, but will initiate if starts undefined
                try:
                    subterms = term._terms
                except AttributeError:
                    subterms = [term]
                self._terms += subterms    # no copies because never modified in-place
    def __mul__(self, x):
        return tensor_sum(self._backend, terms=[x*term for term in self._terms])    # changes scalar prefactors, not raw tensors.  forces copies
    def __getitem__(self, indices):
        indexed_tensors = [term[indices] for term in self._terms]
        if any(isinstance(index,slice) for index in indices):
            new = tensor_sum(self._backend, terms=indexed_tensors)
        else:
            new = sum(indexed_tensors)    # should be a list of scalars if we get here
        return new
    def __str__(self):
        str_terms = "\n+\n".join([str(term) for term in self._terms])
        return f"tensornet.tensor_sum(backend = {self._backend.name}\n{str_terms}\n)"
    def _evaluate(self):
        raw_result = None
        for term in self._terms:
            if term._scalar!=0:
                if raw_result is None:  raw_result = raw(term)    # will make copy of internal data of a primitive_tensor
                else:                   increment(raw_result, term)
        if raw_result is None:
            if self._shape is not None:    # might have been given explicitly or taken from a tensor that was multiplied by zero
                raw_result = self._backend.zeros(self._shape)
            else:
                raise ValueError("Cannot evaluate empty tensor_sum because no dimension information.\nPerhaps use primitive_tensor_factory.zeros(shape).")
        return primitive_tensor(self._backend, raw_result)



# Barebones theory (written much later after a forensic debug battle).  A tensor_network object contains
# two fundamental pieces of information (and other incidental info).  One is a list of contractions.
# Each contraction itself is a list of two-tuples; each two-tuple identifies a tensor and the index
# of that tensor involved in the contraction.  So if a contraction contains two two-tuples, then 
# two indices are contracted with one another, but there might be more for unusual contractions.
# The ordering of the list of contractions is irrelevant, as is the ordering of the list of two-tuples
# that defines a given contraction. The other piece of information is a list of free indices of the result,
# which correspond to uncontracted indices of the tensors in the network.  This list is ordered, and each
# free index is itself a list of two-tuples with the same tensor-index structure as the two-tuples that
# define a contraction.  In the most usual case, each such list corresponding to a free index will be of
# length one, but if there is more than one tensor index that corresponds to a single free index it is because
# those indices are set equal to each other and reduced to a single free index.

class tensor_network(resolved_tensor):
    def __init__(self, backend, scalar, contractions, free_indices):
        values = network_logic.logic.init(backend, scalar, contractions, free_indices)
        resolved_tensor.__init__(self, values.backend, values.shape)
        self._scalar = values.scalar
        self._contractions = values.contractions
        self._free_indices = values.free_indices
        self._hashable = values.hashable
        self._result_hash = values.result_hash
    def __mul__(self, x):
        new = tensor_network(self._backend, self._scalar, self._contractions, self._free_indices)
        new._scalar *= x
        return new
    def __getitem__(self, indices):
        return network_logic.logic.subscript(self, indices, tensor_network)
    def _evaluate(self):
        return network_logic.logic.evaluate(self, tensor_network, primitive_tensor, _backend_contract_path)
