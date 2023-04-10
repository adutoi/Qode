#  (C) Copyright 2023 Anthony D. Dutoi
#
#  This file is part of TonyUtil.
#
#  TonyUtil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Check that I am not duplicating things in itertools

def recursive_looper(loops, kernel, order_aware=False, _arguments=None, _level=0):
    """ loop over the arguments of the function 'kernel' in the order and ranges given in 'loops' """
    # loops  is a list of 2-tuples, where the first element of the list is a pair that
    #        specifies (1) the position of the argument of kernel to iterate over in the
    #        outermost loop and (2) an iterable (like a range object) yielding the respective
    #        argument values to iterate over, and so on (for further loops elements).
    # kernel is a function (probably a wrapped version of what the user is actually
    #        intetested in) that is called in the inner-most loop, with the argument-value
    #        combintation established by the looping.
    if _arguments is None:  _arguments = [None]*len(loops)    # initialize at top level
    if len(loops)>0:  # use >1 and move innermost loop to 'else' to minimize function calls?
        position,values = loops[0]
        loops = loops[1:]
        for v in values:
            if order_aware:
                _arguments[position] = (_level, v)
            else:
                _arguments[position] = v
            recursive_looper(loops, kernel, order_aware=order_aware, _arguments=_arguments, _level=_level+1)
    else:
        kernel(*_arguments)

def compound_range(dims, inactive=None):
    """ returns a list, in which each element is itself a list representing a successive value of a compound index """
    # dims     is a list whose elements are each an iterable (like a range object) that yields the values that
    #          that each respective component of the compound index can take on.
    # inactive is a list that positionally indicates components of the returned compound indices that are not
    #          iterated over.  Their values are replaced by None in the returned list of compound values.
    if inactive is None:  inactive = []
    values = [[]]
    for i,d in reversed(list(enumerate(dims))):
        if i in inactive:
            values = [[None]+v for v in values]
        else:
            new_values = []
            for j in d:
                new_values += [[j]+v for v in values]
            values = new_values
    return values
