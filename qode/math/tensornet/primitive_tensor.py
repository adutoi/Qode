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

from textwrap import indent



# The tensornet type for the primitive tensors that the user sees and uses and builds networks from.
# It should mostly act like the raw tensor type of the chosen backend (ie, forwarding calls), but with
# some extra tensornet-specific functionality that makes the names of the info we need uniform
# across backends.  In order to accomplish this, the tensor importantly knows its backend (via a
# provided module, some of which is implemented by the user if not already provided for that backend type).  
class primitive_tensor(object):
    def __init__(self, raw_tensor, backend):
        # long member names to stay out of way of __getattr__ (never seen by user)
        self.tensornet_raw_tensor = raw_tensor             # only ever accessed by wrapper below
        self.tensornet_backend    = backend
        self.tensornet_shape      = backend.shape(raw_tensor)
    def tensornet_scalar_value(self):                      # only ever accessed by wrapper below
        return self.tensornet_backend.scalar_value(self.tensornet_raw_tensor)
    def tensornet_evaluate(self):
        return self
    def __getattr__(self, name):
        return getattr(self.tensornet_raw_tensor, name)    # forward all other calls for user convenience
    # calls to __XXX__ functions are not passed through by __getattr__ ... for good reason (would return wrong types)
    def __str__(self):
        return "tensornet.primitive_tensor(\n{}\n)".format(indent(str(self.tensornet_raw_tensor), "    "))
    def _check_args(self, other, op):
        try:
            other_raw     = other.tensornet_raw_tensor
            other_backend = other.tensornet_backend
        except:
            raise TypeError("unsupported operand type(s) for {}: \'tensornet.primitive_tensor\' and \'{}\'".format(op, type(other)))
        if other_backend is not self.tensornet_backend:
            raise ValueError("primitive tensors with differing backends cannot be added")
    def __add__(self, other):
        self._check_args(other, "+")
        return primitive_tensor(self.tensornet_raw_tensor + other_raw, self.tensornet_backend)
    def __sub__(self, other):
        self._check_args(other, "-")
        return primitive_tensor(self.tensornet_raw_tensor - other_raw, self.tensornet_backend)


# This would be unnecessary except that we will need a hashable value to keep track of the tensors in 
# a tensor network.  For this we use the builtin id() function, and this would mean that multiple
# occurances of the same tensor in a contraction network would erroneously get the same hash.
# So we temporarily wrap every instance of a primitive_tensor that is given to make it a distinct object.
# While we are at it, we make expose the (very limited) member information that we will need
# about this tensor in a much shortened syntax.
class thin_wrapper(object):
    def __init__(self, prim_tensor):
        self._prim_tensor = prim_tensor
        self.raw          = prim_tensor.tensornet_raw_tensor    #
        self.shape        = prim_tensor.tensornet_shape         # only accessed when we are sure
    def scalar_value(self):                                     # that it is a wrapped primitive
        return self._prim_tensor.tensornet_scalar_value()       #



# This is the only function in this module available to users.  With a syntax like
#    tensor = primitive_tensor_wrapper(tensornet.numpy_backend)
#    A = tensor(numpy.array(...))
# the user can easily wrap all the raw tensors they like in their chosen backend.
def primitive_tensor_wrapper(backend):
    def wrapper(raw_tensor):
        return primitive_tensor(raw_tensor, backend)
    return wrapper
