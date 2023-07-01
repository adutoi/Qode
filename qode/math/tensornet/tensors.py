#    (C) Copyright 2023 Anthony D. Dutoi
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
from .base    import tensor_base, shape, evaluate, increment, extract, _scalar, _raw, _backend



# adds general addition to any class above the base class (needs tensor_sum here for its definition)
class summable_tensor(tensor_base):
    def __add__(self, other):
        return tensor_sum([self, other])

# in an informational sense, this is the most basic tensor, since all it does is store a list
# of tensor_networks and primitive_tensors.
class tensor_sum(summable_tensor):
    def __init__(self, tensor_terms=None):
        self.tensornet_backend = None
        self.tensornet_shape   = None
        self.tensor_terms      = []
        if tensor_terms==None:  tensor_terms = []    # for instantiation of empty sum as accumulator
        for term in tensor_terms:
            if self.tensornet_backend is None:
                try:
                    self.tensornet_backend = _backend(term)
                    self.tensornet_shape   = shape(term)
                except:
                    raise TypeError("only tensornet tensors can be summed")
            if _backend(term) is not self.tensornet_backend or shape(term)!=self.tensornet_shape:
                raise ValueError("only tensornet tensors with the same backend and shape can be summed")
            try:
                tensor_subterms = term.tensor_terms
            except AttributeError:
                new_terms = [copy(term)]                                        # we want copies, ...
            else:
                new_terms = [copy(sub_term) for sub_term in tensor_subterms]    # ... in case we use *=
            self.tensor_terms += new_terms
    def tensornet_evaluate(self):
        result = extract(self.tensor_terms[0])
        for term in self.tensor_terms[1:]:
            increment(result, term)    # move actual math out of here and let child classes decided how to add
        return primitive_tensor(result, _backend(self))
    def __getitem__(self, indices):
        indexed_tensors = [tens[indices] for tens in self.tensor_terms]
        if any(isinstance(index,slice) for index in indices):
            new = tensor_sum(indexed_tensors)
        else:
            new = sum(indexed_tensors)    # should be a list of scalars if we get here
        return new
    def __imul__(self, x):             # enables __mul__, __rmul__, __neg__, and therefore also __sub__
        for term in self.tensor_terms:
            term *= x                  # will change only scalar prefactors, not raw tensors inside of primitive_tensors
        return self
    # extra functionality just for tensor_sum
    def __iadd__(self, other):
        try:
            other_backend = _backend(other)
        except:
            raise TypeError("only tensornet tensors can be added to a tensornet tensor_sum")
        if other_backend is not _backend(self):
            raise ValueError("only tensornet tensors with the same backend can be added")
        try:
            other_tensor_terms = other.tensor_terms
        except AttributeError:
            new_terms = [copy(other)]                                              # we want copies, ...
        else:
            new_terms = [copy(other_term) for other_term in other_tensor_terms]    # ... in case we use *=
        self.tensor_terms += new_terms
        return self
    def __isub__(self, other):         # enabled by __iadd__ and (indirectly) __imul__
        self += -other
        return self



def primitive_tensor_wrapper(backend, copy_data=False):
    def wrapper(raw_tensor):
        return primitive_tensor(raw_tensor, backend, copy_data)
    return wrapper

# The tensornet type for the primitive tensors that the user sees and uses and builds networks from.
# The tensor importantly knows its backend, via a provided module (implemented by the user if not
# already provided for that backend type).
# only users ever use copy_data because we never modify here that which we wrap in a primitive_tensor
# only this class uses _scalar.  Users should use * or *=
class primitive_tensor(summable_tensor):
    def __init__(self, raw_tensor, backend, copy_data=False, _scalar=1):
        self.tensornet_backend    = backend
        self.tensornet_shape      = backend.shape(raw_tensor)
        self.tensornet_scalar     = _scalar    # here so that we can define *= without changing original data
        if copy_data:
            self.tensornet_raw_tensor = backend.copy_data(raw_tensor)
        else:
            self.tensornet_raw_tensor = raw_tensor
    def tensornet_evaluate(self):
        if _scalar(self)==1:
            return primitive_tensor(_raw(self), _backend(self))                  # this might get modified and so should be independent ...
        else:
            return primitive_tensor(_scalar(self)*_raw(self), _backend(self))
    def tensornet_increment(self, result):
        if _scalar(self)==1:
            result += _raw(self)                                                 # ... but do not make a copy just to use as an increment
        else:
            result += _scalar(self) * _raw(self)
        return
    def __getitem__(self, indices):
        indexed_tensor = self.tensornet_raw_tensor[indices]
        if any(isinstance(index,slice) for index in indices):
            new = primitive_tensor(indexed_tensor, self.tensornet_backend, _scalar=self.tensornet_scalar)
        else:
            new = self.tensornet_scalar * indexed_tensor    # this is a scalar if we get here
        return new
    # __iadd__ and __isub__ would be confusing since __add__ and __sub__ make a tensor_sum (*)
    # but the increment operators should be of input type.  (* this is for the best because it
    # is more flexible; the user can choose to do a hard data-level add outside of tensornet)
    def __imul__(self, x):    # enables __mul__, __rmul__, __neg__, and therefore also __sub__
        self.tensornet_scalar *= x
        return self
    # extra functionality just for primitive_tensor
    def __str__(self):
        return "tensornet.primitive_tensor(\n{}\n)".format(indent(str(_raw(self)), "    "))
