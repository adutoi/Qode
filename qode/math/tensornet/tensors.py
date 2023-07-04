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
from .base    import tensor_base, increment, extract, _resolve_contract_ops



# adds general addition to any class above the base class (needs tensor_sum here for its definition)
class summable_tensor(tensor_base):
    def __init__(self, shape, backend, contract):
        tensor_base.__init__(self, shape, backend, contract)
    def __add__(self, other):
        return tensor_sum([self, other])

# in an informational sense, this is the most basic tensor, since all it does is store a list
# of tensor_networks and primitive_tensors.
class tensor_sum(summable_tensor):
    def __init__(self, tensor_terms=None):
        summable_tensor.__init__(self, None, None, None)
        self._tensor_terms = []
        if tensor_terms==None:  tensor_terms = []    # for instantiation of empty sum as accumulator
        for term in tensor_terms:
            term_ = _resolve_contract_ops(term)
            if self.shape is None:
                try:
                    self.shape     = term_.shape        # A little dirty to
                    self._backend  = term_._backend     # mess with these
                    self._contract = term_._contract    # directly like this (rather than via base class)
                except:
                    raise TypeError("only tensornet tensors can be summed")
            if term_._backend is not self._backend or term_.shape!=self.shape:
                raise ValueError("only tensornet tensors with the same backend and shape can be summed")
            try:
                tensor_subterms = term_._tensor_terms
            except AttributeError:
                new_terms = [copy(term_)]                                        # we want copies, ...
            else:
                new_terms = [copy(sub_term) for sub_term in tensor_subterms]    # ... in case we use *=
            self._tensor_terms += new_terms
    def _increment(self, result):
        result += extract(self)
        return
    def _evaluate(self):
        result = extract(self._tensor_terms[0])
        for term in self._tensor_terms[1:]:
            increment(result, term)       # move actual math out of here and let child classes decided how to add
        return primitive_tensor(result, self._backend, self._contract)
    def __getitem__(self, indices):
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
        other = _resolve_contract_ops(other)
        try:
            other_backend = other._backend
        except:
            raise TypeError("only tensornet tensors can be added to a tensornet tensor_sum")
        if len(self._tensor_terms)==0 and self._backend is None:
            self._backend = other_backend    # must have started as an empty accumulator
        if other_backend is not self._backend:
            raise ValueError("only tensornet tensors with the same backend can be added")
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
# only users ever use copy_data because we never modify here that which we wrap in a primitive_tensor
# only this class uses _scalar.  Users should use * or *=
class primitive_tensor(summable_tensor):
    def __init__(self, raw_tensor, backend, contract, copy_data=False, _scalar=1):
        summable_tensor.__init__(self, backend.shape(raw_tensor), backend, contract)
        self._scalar  = _scalar    # here so that we can define *= without changing original data
        if copy_data:
            self._raw_tensor = backend.copy_data(raw_tensor)
        else:
            self._raw_tensor = raw_tensor
    def _increment(self, result):
        if self._scalar==1:
            result += self._raw_tensor                                                           # do not make a copy just to use as an increment, but we want to ...
        else:
            result += self._scalar * self._raw_tensor
        return
    def _evaluate(self):
        return primitive_tensor(self._scalar*self._raw_tensor, self._backend, self._contract)    # ... copy the data in case someone (like tensor_sum) modifies the result of extract()
    def __getitem__(self, indices):
        indexed_tensor = self._raw_tensor[indices]
        if any(isinstance(index,slice) for index in indices):
            new = primitive_tensor(indexed_tensor, self._backend, self._contract, _scalar=self._scalar)
        else:
            new = self._scalar * indexed_tensor    # this is a scalar if we get here
        return new
    # __iadd__ and __isub__ would be confusing since __add__ and __sub__ make a tensor_sum (*)
    # but the increment operators should be of input type.  (* this is for the best because it
    # is more flexible; the user can choose to do a hard data-level add outside of tensornet)
    def __imul__(self, x):    # enables __mul__, __rmul__, __neg__, and therefore also __sub__
        self._scalar *= x
        return self
    # extra functionality just for primitive_tensor
    def __str__(self):
        return "tensornet.primitive_tensor(\n{}\n)".format(indent(str(self._raw_tensor), "    "))
