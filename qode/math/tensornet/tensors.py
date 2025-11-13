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

from copy     import copy
from textwrap import indent
from .base    import increment, raw, resolve, resolve_ellipsis, to_contract
from . import network_logic

_backend_contract_path = False    # if True, let backend handle finding the optimal contraction path upon evaluate() call
def backend_contract_path(TrueFalse):
    global _backend_contract_path
    _backend_contract_path = TrueFalse



# expected to have a backend and a shape (and a scalar if not a sum), everything else is specific to single class
class tensor_base(object):
    def __init__(self, shape, backend):
        self.shape     = shape
        self._backend  = backend
    def __setitem__(self, item):
        raise RuntimeError("elements of tensornet tensors are not assignable")
    def __call__(self, *indices):
        return to_contract(self, indices, tensor_network)
    def _resolve(self):
        return self
    def _raw(self):
        return self._evaluate()._raw_tensor
    def _scalar_value(self):
        if len(self.shape)>0:
            raise RuntimeError("cannot take the scalar value of a tensornet tensor with >0 free indices")
        return self._backend.scalar_value(self._raw())
    # __mul__, __rmul__, __neg__, and __sub__ use __imul__ from child; multiplication here is always with a scalar
    def __itruediv__(self, x):
        self *= (1./x)
        return self
    def __mul__(self, x):
        new = copy(self)
        new *= x
        return new
    def __truediv__(self, x):
        return self * (1./x)
    def __rmul__(self, x):
        return self * x
    def __neg__(self):
        return self * -1
    def __add__(self, other):
        return tensor_sum([self, other])
    def __sub__(self, other):
        return self + (-other)


# in an informational sense, this is the most basic tensor, since all it does is store a list
# of tensor_networks and primitive_tensors.
class tensor_sum(tensor_base):
    def __init__(self, tensor_terms=None):
        tensor_base.__init__(self, None, None)
        if tensor_terms==None:  tensor_terms = []    # for instantiation of empty sum as accumulator
        self._tensor_terms = []
        for term in tensor_terms:
            term_ = resolve(term)
            if self.shape is None:
                try:
                    self.shape     = term_.shape        # A little dirty to
                    self._backend  = term_._backend     # mess with these
                except:
                    raise TypeError("only tensornet tensors can be summed")
            if (term_._backend is not self._backend) or (term_.shape!=self.shape):
                raise ValueError("only tensornet tensors with the same backend and shape can be summed")
            try:
                tensor_subterms = term_._tensor_terms
            except AttributeError:
                new_terms = [copy(term_)]                                       # we want copies, ...
            else:
                new_terms = [copy(sub_term) for sub_term in tensor_subterms]    # ... in case we use *=
            self._tensor_terms += new_terms
    def _increment(self, result):
        self._backend.increment(result, raw(self))
        return
    def _evaluate(self):
        result = None
        for term in self._tensor_terms:
            if term._scalar!=0:
                if result is None:  result = raw(term)
                else:               increment(result, term)    # move actual math out of here and let child classes decided how to add
        if result is None:
            try:
                result = raw(self._tensor_terms[0])                # will produce zero tensor of correct dimensions (should I used primitive_tensor.zeros here?)
            except IndexError:
                raise ValueError("cannot evaluate empty tensor_sum because no dimension information.  perhaps use primitive_tensor.zeros(shape).")
        return primitive_tensor(result, self._backend)
    def __copy__(self):
        return tensor_sum(self._tensor_terms)    # makes a copy of list with copies of terms (bc both modified by += and *=)
    def __getitem__(self, indices):
        indices = resolve_ellipsis(indices, self.shape)
        indexed_tensors = [tens[indices] for tens in self._tensor_terms]
        if any(isinstance(index,slice) for index in indices):
            new = tensor_sum(indexed_tensors)
        else:
            new = sum(indexed_tensors)    # should be a list of scalars if we get here
        return new
    def __imul__(self, x):                # enables __mul__, __rmul__, __neg__, and therefore also __sub__
        for term in self._tensor_terms:
            term *= x                     # will change only scalar prefactors, not raw tensors inside of primitive_tensors
        return self
    # extra functionality just for tensor_sum
    def __iadd__(self, other):
        # duplicates some code in __init__.  how to combine?
        other = resolve(other)
        try:
            other_shape    = other.shape
            other_backend  = other._backend
        except:
            raise TypeError("only tensornet tensors can be added to a tensornet tensor_sum")
        if len(self._tensor_terms)==0 and self._backend is None:    # must have started as an empty accumulator
            self.shape     = other_shape
            self._backend  = other_backend
        if other_backend is not self._backend:
            raise ValueError("only tensornet tensors with the same backend can be added")
        if other_shape!=self.shape:
            raise ValueError("only tensors with equivalent shapes can be summed")
        try:
            other_tensor_terms = other._tensor_terms
        except AttributeError:
            new_terms = [copy(other)]                                              # we want copies, ...
        else:
            new_terms = [copy(other_term) for other_term in other_tensor_terms]    # ... in case we use *=
        self._tensor_terms += new_terms
        return self
    def __isub__(self, other):            # enabled by __iadd__ and (indirectly) __imul__
        self += -other
        return self



# The tensornet type for the primitive tensors that the user sees and uses and builds networks from.
# The tensor importantly knows its backend, via a provided module (implemented by the user if not
# already provided for that backend type).
# Only this class uses _scalar.  Use * or *= from outside the class.
class primitive_tensor(tensor_base):
    def __init__(self, raw_tensor, backend, _scalar=1):
        tensor_base.__init__(self, backend.shape(raw_tensor), backend)
        self._scalar  = _scalar    # here so that we can define *= without changing original data
        self._raw_tensor = raw_tensor
    @staticmethod
    def zeros(shape, backend):
        return primitive_tensor(backend.zeros(shape), backend)
    #@staticmethod
    #def scalar_tensor(scalar, backend):
    #    return primitive_tensor(backend.scalar_tensor(scalar), backend)
    #def _raw(self):
    #    return self._evaluate()._raw_tensor
    def _increment(self, result):
        if self._scalar==1:
            self._backend.increment(result, self._raw_tensor)    # do not make a copy just to use as an increment, but we want to ...
        else:
            self._backend.increment(result, self._backend.mult(self._scalar, self._raw_tensor))
        return
    def _evaluate(self):
        return primitive_tensor(self._backend.mult(self._scalar, self._raw_tensor), self._backend)    # ... copy the data in case someone (like tensor_sum) modifies the result of raw()
    def __getitem__(self, indices):
        indices = resolve_ellipsis(indices, self.shape)
        indexed_tensor = self._backend.element(self._raw_tensor, indices)
        if any(isinstance(index,slice) for index in indices):
            new = primitive_tensor(indexed_tensor, self._backend, _scalar=self._scalar)
        else:
            new = self._scalar * indexed_tensor    # this is a scalar if we get here
        return new
    # __iadd__ and __isub__ would be confusing since __add__ and __sub__ make a tensor_sum (*)
    # but the increment operators should be of input type.  (* this is for the best because it
    # is more flexible; the user can choose to do a hard data-level add outside of tensornet)
    def __imul__(self, x):    # enables __mul__, __rmul__, __neg__, and therefore also __sub__
        self._scalar *= x
        return self
    def __iadd__(self, _):
        raise NotImplementedError("+= and -= do not exist for primitive_tensor because value of tensornet is lazy evaluation via tensor_sum.\nUse + and - or empty/zero tensor_sum for accumulation.")
    def __isub__(self, _):
        self += None    # just raises NotImplementedError
    # extra functionality just for primitive_tensor
    def __str__(self):
        return "tensornet.primitive_tensor(\n{}\n)".format(indent(self._backend.str(self._raw_tensor), "    "))






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

class tensor_network(tensor_base):
    def __init__(self, scalar, contractions, free_indices, backend):
        values = network_logic.logic.init(scalar, contractions, free_indices, backend)
        tensor_base.__init__(self, values.shape, values.backend)
        self._scalar = values.scalar
        self._contractions = values.contractions
        self._free_indices = values.free_indices
        self._hashable = values.hashable
        self._result_hash = values.result_hash
    @staticmethod
    def build(*tensor_factors):
        values = network_logic.logic.build(*tensor_factors)
        return tensor_network(**values)
    def _increment(self, result):
        self._backend.increment(result, raw(self))
        return
    def _evaluate(self):
        return network_logic.logic.evaluate(self, tensor_network, primitive_tensor, _backend_contract_path)
    def __getitem__(self, indices):
        return network_logic.logic.element(self, indices, tensor_network)
    def __imul__(self, x):
        self._scalar *= x
        return self
