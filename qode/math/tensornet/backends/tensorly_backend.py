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

def copy_data(tensor):
    return tensorly.copy(tensor)    # Badly implemented in tensorly; raises warning for pytorch.

def scalar_value(tensor):
    return tensor.item()            # Object-bound method works for numpy and pytorch, but cannot find generic tensorly wrapper.

def scalar_tensor(scalar):
    return tensorly.tensor(scalar)

def shape(tensor):
    return tensorly.shape(tensor)

def contract(*tensor_factors):
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
    # print("einsum called with", instructions)
    value = scalar * tensorly.einsum(instructions, *tensors)    # works with numpy and pytorch because einsum is available ... write more generally if needed
    # print("value id", id(value))
    return value
