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
from ...util import timer



def shape(tensor):    # just in case a user expects such an unbound function to exist
    return tensor.shape

def resolve(tensor):
    return tensor._resolve()

def evaluate(tensor):
    return tensor._evaluate()

def raw(tensor):
    return tensor._raw()

def scalar_value(tensor):
    return tensor._scalar_value()

def increment(raw_result, tensor):    # for incrementing raw tensors of the same shape (implicitly evaluates)
    return tensor._increment(raw_result)




def resolve_ellipsis(indices, shape):
    ellipsis_resolution = [slice(None)] * (len(shape) - len(indices) + 1)
    new_indices = []
    found_ellipsis = False
    for index in indices:
        if index is Ellipsis:
            if found_ellipsis:  raise IndexError("indices can only have a single ellipsis")
            new_indices += ellipsis_resolution
            found_ellipsis = True
        else:
            new_indices += [index]
    return tuple(new_indices)




class to_contract(object):
    def __init__(self, tensor, indices, tensor_network, _from_list=None):
        if _from_list is None:
            self._tensors = [(tensor, indices)]
        else:
            self._tensors = list(_from_list)    # for internal use only.  ignores first two args
        self._tensor_network = tensor_network
    def divulge(self):    # logically only called from _tensor_network when self._tensors is of length 1
        return self._tensors[0]
    def _evaluate(self):
        return self._resolve()._evaluate()
    def _raw(self):
        return self._resolve()._raw()
    def _scalar_value(self):
        return self._resolve()._scalar_value()
    def _increment(self, result):
        return self._resolve()._increment(result)
    def _resolve(self):
        def _resolve_sum(tensor_factors):
            # This function turns a contractions of sums into sums of contractions and passes the terms off to self._tensor_network
            # Implicit type checking is essentially delegated to self._tensor_network.
            outer_terms = [[]]
            for factor in tensor_factors:
                try:
                    tens, indices = factor.divulge()
                    inner_factors = tens._tensor_terms
                except:
                    for outer_term in outer_terms:
                        outer_term += [factor]
                else:
                    new_outer_terms = []
                    for outer_term in outer_terms:
                        for inner_factor in inner_factors:
                            new_outer_terms += [outer_term + [to_contract(inner_factor, indices, self._tensor_network)]]
                    outer_terms = new_outer_terms
            tensor_terms = []
            for outer_term in outer_terms:
                tensor_terms += [self._tensor_network.build(*outer_term)]
            #
            #result_hashes = sorted((term._result_hash, i) for i,term in enumerate(tensor_terms))
            #tensor_terms, tensor_terms_ = [], tensor_terms
            #previous = None
            #for result_hash,i in result_hashes:
            #    if result_hash==previous:
            #        tensor_terms[-1]._scalar += tensor_terms_[i]._scalar
            #    else:
            #        tensor_terms += [tensor_terms_[i]]
            #    previous = result_hash
            #
            if len(tensor_terms)==1:
                return tensor_terms[0]
            else:
                the_sum = tensor_terms[0] + tensor_terms[1]
                for term in tensor_terms[2:]:
                    the_sum += term
                return the_sum
        args = (to_contract(*tensor, self._tensor_network) for tensor in self._tensors)
        return _resolve_sum(args)
    def __getattr__(self, attr):    # only called for non hard-coded attributes
        # doing it this way instead of defining method named shape lets call be tens.shape instead of tens.shape() [just set it in __init__?]
        if attr=="shape":
            return self._resolve().shape
        if attr=="_backend":
            return self._resolve()._backend
        else:
            raise AttributeError("'to_contract' object has no attribute '{}'".format(attr))
    def __call__(self, *indices):
        return self._resolve()(*indices)
    def __setitem__(self, item):
        raise RuntimeError("elements of tensornet tensors are not assignable")
    def __getitem__(self, indices):
        return self._resolve()[indices]
    def __imul__(self, other):
        try:
            other_tensors = other._tensors
        except AttributeError:
            self._tensors += [(other, None)]    # assume it is a scalar.  means pure outer pdt must be written as A() @ B()
        else:
            raise TypeError("use @ or @= to join tensors via contraction or outer product")
        return self
    def __imatmul__(self, other):
        try:
            other_tensors = other._tensors
        except AttributeError:
            raise TypeError("use * or *= for multiplication by a scalar")
        else:
            self._tensors += other_tensors 
        return self
    def __itruediv__(self, x):
        self *= (1./x)
        return self
    def __mul__(self, other):
        new = to_contract(None, None, self._tensor_network, _from_list=self._tensors)
        new *= other
        return new
    def __matmul__(self, other):
        new = to_contract(None, None, self._tensor_network, _from_list=self._tensors)
        new @= other
        return new
    def __truediv__(self, x):
        return self * (1./x)
    def __rmul__(self, x):          # only needed for leading scalars
        return self * x
    def __neg__(self):
        return self * -1
    def __add__(self, other):
        return self._resolve() + other._resolve()
    def __sub__(self, other):
        return self + (-other)
    def __iadd__(self, _):
        raise NotImplementedError("+= and -= do not exist for to_contract")
    def __isub__(self, _):
        self += None    # just raises NotImplementedError






_timings = None

def initialize_timer():    # calling more than once just clears out the old timer
    global _timings
    _timings = timer()

def print_timings(header=None):
    global _timings
    if header is None:  header = "tensornet contraction engine"
    _timings.print(header)

def timings_start():
    global _timings
    if _timings is not None:  _timings.start()

def timings_record(label):
    global _timings
    if _timings is not None:  _timings.record(label)
