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

from copy import copy
from ...util import timer



def evaluate(tensor):
    return _resolve_contract_ops(tensor)._evaluate()

def increment(result, tensor):    # for incrementing raw tensors of the same shape (implicitly evaluates)
    return _resolve_contract_ops(tensor)._increment(result)

def raw(tensor):
    return evaluate(tensor)._raw_tensor

def scalar_value(tensor):
    tensor = _resolve_contract_ops(tensor)
    if len(tensor.shape)>0:
        raise RuntimeError("cannot take the scalar value of a tensornet tensor with >0 free indices")
    return tensor._backend.scalar_value(raw(tensor))



def _resolve_contract_ops(tensor):
    try:
        temp = tensor._call_contract()
    except AttributeError:
        temp = tensor
    return temp

class to_contract(object):
    def __init__(self, tensor, indices, contract, _from_list=None):
        if _from_list is None:
            self._tensors = [(tensor, indices)]
        else:
            self._tensors = list(_from_list)    # for internal use only.  ignores first two args
        self._contract = contract
    def divulge(self):    # logically only called from _contract when self._tensors is of length 1
        return self._tensors[0]
    def _call_contract(self):
        return self._contract(*(to_contract(*tensor, self._contract) for tensor in self._tensors))
    def __getattr__(self, attr):
        if attr=="shape":    # should only be called for non hard-coded attributes
            return self._call_contract().shape
        else:
            raise AttributeError("'to_contract' object has no attribute '{}'".format(attr))
    def __call__(self, *indices):
        return self._call_contract()(*indices)
    def __setitem__(self, item):
        raise RuntimeError("elements of tensornet tensors are not assignable")
    def __getitem__(self, indices):
        return self._call_contract()[indices]
    def __imul__(self, other):
        try:
            other_tensors = other._tensors
        except AttributeError:
            self._tensors += [(other, None)]    # assume it is a scalar.  means pure outer pdt must be written as A() * B()
        else:
            raise TypeError("use @ operator to join tensors via contraction or outer product")
        return self
    def __imatmul__(self, other):
        try:
            other_tensors = other._tensors
        except AttributeError:
            raise TypeError("use * operator for multiplication by a scalar")
        else:
            self._tensors += other_tensors 
        return self
    def __itruediv__(self, x):
        self *= (1./x)
        return self
    def __mul__(self, other):
        new = to_contract(None, None, self._contract, _from_list=self._tensors)
        new *= other
        return new
    def __matmul__(self, other):
        new = to_contract(None, None, self._contract, _from_list=self._tensors)
        new @= other
        return new
    def __truediv__(self, x):
        return self * (1./x)
    def __rmul__(self, x):          # only needed for leading scalars
        return self * x
    def __neg__(self):
        return self * -1
    def __add__(self, other):
        return self._call_contract() + other._call_contract()
    def __sub__(self, other):
        return self + (-other)



# expected to have a backend and a shape (and a scalar if not a sum), everything else is specific to single class
class tensor_base(object):
    def __init__(self, shape, backend, contract):
        self.shape     = shape
        self._backend  = backend
        self._contract = contract    # inject the top-most function in the module so that they can contract "themselves"
    def __setitem__(self, item):
        raise RuntimeError("elements of tensornet tensors are not assignable")
    def __call__(self, *indices):
        return to_contract(self, indices, self._contract)
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
    # the following needs __add__ but cannot define here because no knowledge of tensor_sum
    def __sub__(self, other):
        return self + (-other)



timings = None

def initialize_timer():    # calling more than once just clears out the old timer
    global timings
    timings = timer()

def print_timings(header=None):
    global timings
    if header is None:  header = "tensornet contraction engine"
    timings.print(header)

def timings_start():
    global timings
    if timings is not None:  timings.start()

def timings_record(label):
    global timings
    if timings is not None:  timings.record(label)
