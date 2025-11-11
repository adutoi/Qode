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
import tensorly
from ....util import timer



_timings = None

def initialize_timer():    # calling more than once just clears out the old timer
    global _timings
    _timings = timer()

def print_timings(header=None):
    global _timings
    if header is None:  header = "tensorly_backend module"
    _timings.print(header)



# Need to place these in an object (a class is an object) because modules cannot be pickled.
# On a similar note, they need to be in a class and not a class instance because instances
# that are the same upon pickling may be unpickled as distinct objects, which interferes
# with checking if two backends are the same (best to used "is" for absolute assurance that
# they are exactly the same).  That said, since modules are not pickled (referenced by fully
# qualified path and loaded only once), we can check that two backends are the same using the
# "is" operator between the class defined here.

class functions(object):
    @staticmethod
    def copy_data(tensor):
        return tensorly.copy(tensor)    # Badly implemented in tensorly; raises warning for pytorch.
    @staticmethod
    def scalar_value(tensor):
        return tensor.item()            # Object-bound method works for numpy and pytorch, but cannot find generic tensorly wrapper.
    @staticmethod
    def scalar_tensor(scalar):
        return tensorly.tensor(scalar)
    @staticmethod
    def shape(tensor):
        return tensorly.shape(tensor)
    @staticmethod
    def zeros(shape):
        return tensorly.zeros(shape)
    @staticmethod
    def increment(tensor, delta):
        tensor += delta
        return
    @staticmethod
    def mult(scalar, tensor):
        return scalar * tensor
    @staticmethod
    def element(tensor, indices):
        return tensor[indices]
    @staticmethod
    def str(tensor):
        return str(tensor)
    @staticmethod
    def contract(*tensor_factors):
        global _timings
        if _timings is not None:  _timings.start()
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
        if _timings is not None:  _timings.record("admin")
        # print("einsum called with", instructions)
        if _timings is not None:  _timings.start()
        value = scalar * tensorly.einsum(instructions, *tensors)    # works with numpy and pytorch because einsum is available ... write more generally if needed
        if _timings is not None:  _timings.record("einsum")
        # print("value id", id(value))
        return value
