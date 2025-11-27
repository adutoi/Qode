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

_timings = None

class ContractionError(ValueError):    # this is mostly for debugging
    pass

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

def increment(raw_result, tensor):    # for incrementing raw tensors of the same shape (for efficiency)
    return tensor._increment(raw_result)

def subscript(tensor, indices):
    return tensor._subscript(indices)

def resolve_ellipsis(indices, n_dim):
    ellipsis_resolution = [slice(None)] * (n_dim - len(indices) + 1)
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
