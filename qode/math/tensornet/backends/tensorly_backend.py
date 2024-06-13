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

import tensorly
from ....util import timer



timings = None

def initialize_timer():    # calling more than once just clears out the old timer
    global timings
    timings = timer()

def print_timings(header=None):
    global timings
    if header is None:  header = "tensorly_backend module"
    timings.print(header)



def copy_data(tensor):
    return tensorly.copy(tensor)    # Badly implemented in tensorly; raises warning for pytorch.

def scalar_value(tensor):
    return tensor.item()            # Object-bound method works for numpy and pytorch, but cannot find generic tensorly wrapper.

def scalar_tensor(scalar):
    return tensorly.tensor(scalar)

def shape(tensor):
    return tensorly.shape(tensor)

def contract(*tensor_factors):
    global timings
    if timings is not None:  timings.start()
    ####
    # args = []
    # for factor in tensor_factors:
    #     try:
    #         tensor,*indices = factor
    #     except:
    #         args += [factor]
    #     else:
    #         args += [(id(tensor),*indices)]
    # print("backend.contract called with", *args)
    ####
    def letters(excluded):
        i = 0
        def next_letter():
            nonlocal i
            candidate = 'abcdefghijklmnopqrstuvwxyz'[i]
            i += 1
            if candidate in excluded:
                return next_letter()
            else:
                return candidate
        return next_letter
    tensors = []
    index_strings = []
    all_indices = set()
    max_int = -1
    scalar = 1
    for factor in tensor_factors:
        try:
            tens, *indices = factor
        except:
            scalar *= factor
        else:
            all_indices |= set(indices)
            tensors += [tens]
            index_strings += [indices]
            for index in indices:
                if isinstance(index,int):
                    if index>max_int:
                        max_int = index
    free_indices = []
    if max_int>=0:
        next_letter = letters(all_indices)
        free_indices = [next_letter() for _ in range(max_int+1)]
    instructions = []
    for indices in index_strings:
        for i in range(len(indices)):
            if isinstance(indices[i],int):
                indices[i] = free_indices[indices[i]]
        instructions += ["".join(indices)]
    instructions = ",".join(instructions)
    instructions += "->"
    instructions += "".join(free_indices)
    if timings is not None:  timings.record("admin")
    # print("einsum called with", instructions)
    if timings is not None:  timings.start()
    value = scalar * tensorly.einsum(instructions, *tensors)    # works with numpy and pytorch because einsum is available ... write more generally if needed
    if timings is not None:  timings.record("einsum")
    # print("value id", id(value))
    return value



# Had to do this because a module cannot be pickled so we need our classes to carry around
# objects that can be pickled.  The ID() thing is because we used to compare id of module
# using the "is" operator, but that does not work for objects of this class which have been
# pickled because it unpickles as a distinct object even though it is the same backend!
# However, the module does not get pickled!  So when the object is unpickled it is still of
# this precise class (a class is also an object), so we can compare this id (could probably
# just use a direct reference to the class and still use "is" in calling code ...).
class _functions(object):
    def ID(self):
        return id(_functions)
    def copy_data(self, tensor):
        return copy_data(tensor)
    def scalar_value(self, tensor):
        return scalar_value(tensor)
    def scalar_tensor(self, scalar):
        return scalar_tensor(scalar)
    def shape(self, tensor):
        return shape(tensor)
    def contract(self, *tensor_factors):
        return contract(*tensor_factors)

functions = _functions()
