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




def evaluate(tensor):
    return tensor.tensornet_evaluate()

def increment(result, tensor):
    return tensor.tensornet_increment(result)

def raw(tensor):   # user beware, does not return scalar
    try:
        raw_tensor = tensor.tensornet_raw_tensor
    except AttributeError:
        raise TypeError("can only return raw tensor from a tensornet primitive")
    return raw_tensor

def scalar(tensor):    # user beware, does not return raw
    return tensor.tensornet_scalar

def backend(tensor):
    return tensor.tensornet_backend

def scalar_value(tensor):
    try:
        value = tensor.tensornet_scalar_value()
    except AttributeError:
        raise TypeError("can only attempt scalar evaluation of a tensornet primitive")
    return value

def shape(tensor):
    return tensor.tensornet_shape


# expected to have a backend and a shape (and a scalar if not a sum), everything else is specific to single class
class tensor_base(object):
    # __mul__, __rmul__, __neg__, and __sub__ use __imul__ from child
    # multiplication here is always by a scalar
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
    # needs __add__ but cannot define here because no knowledge of tensor_sum
    def __sub__(self, other):
        return self + (-other)
