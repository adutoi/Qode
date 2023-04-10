#    (C) Copyright 2018 Anthony D. Dutoi
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
# Any time one is writing low-level utilities in python, something is probably wrong.
# These will likely be deprecated, if the below is correct:
#   If one wants the indices of the smallest/largest elements, it is probably because
#   there are two lists with the same ordering, in which case, zip, sort/slice, and built-in
#   min are probably the tools of choice.  If one does not want the index, then sort/slice
#   and built-in min are probably all one needs.
# Consider the above options before using these in the future.

def idx_val(float_list, num=None):
    N = len(float_list)
    n = 1 if num is None else num
    n = min(n,N)
    indices  = []
    smallest = []
    for _ in range(n):
        index = None
        small = float('inf')
        for i in range(N):
            if i not in indices and float_list[i]<small:
                index = i
                small = float_list[i]
        indices  += [index]
        smallest += [small]
    if num is None:
        indices  = indices[0]
        smallest = smallest[0]
    return indices,smallest

def value(float_list, num=None):		# should have same result at built-in min() function for num=None
    _,val = idx_val(float_list,num)
    return val

def index(float_list, num=None):
    idx,_ = idx_val(float_list,num)
    return idx
