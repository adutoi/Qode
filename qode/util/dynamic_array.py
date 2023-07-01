#    (C) Copyright 2019 Anthony D. Dutoi
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

"""\
The class dynamic_array lives somewhere between a dictionary, a tensor, and a function.

For a given rule, given upon creation of a dynamic_array object (or list of rules, which are nested into a single rule),
the object accepts an index in square brackets [] (more generally, indices, given in a tuple) and returns
a value that it determines on the fly, by calling the rule(s) with those same indices (unpacked, if a tuple) as
arguments (so the rule is just a function that accepts a fixed number of arguments ... or a variable
number via the func(*args) syntax for generic wrapper rules).

The most important differnce between dynamic_array and simple calls to a function (namely the rule function) is that 
a dynamic_array has explicit awareness of, and can report on, the values allowed for each of its arguments,
and the anticipated use cases are those for which these values are discrete (and also usually finite in 
number).  For example, it could then be used to generate blocks for block-tensor contractions.

To this end, the object behaves a little bit like a tensor, except that the values are dynamically
generated by a backing function.  Also different from a typical tensor is that the index values 
do not need to be integers, or even have an inherent ordering.  To this extent it is a bit like a
multi-dimensional dynamic dictionary, in which access to *any combination* of the allowed keys is expected to
be *logically* meaningful (which is also a distinction between a dynamic_array and a generic function),
though the values might be zero (or divergent, I guess).

A user must specify at the time of creation the allowed ranges for each of the indices as a list containing
one range for each index.  These ranges may be lists of values, native python range objects or abstractions of a users
choice.  If the user wants a range to be completely unrestricted (no bounds or type checking), then None
may be given.  A suitable range object r should do something sensible for "if i in r" at least, and
len(r), str(r) and indexing/iteration over r are also good things to implement. 
In other words, a user defined option should have __contains__ function defined (and perhaps __len__, __str__, __iter__, __getitem__).
(If ranges are not restricted via giving None, the range is replaced here by a suitable such object).
As a concrete example, of the object represents all integers, __len__, __iter__, and __getitem__ should raise exceptions, __str__ should
just describe that it is the infinite set of integers, and __contains__ should check that the object is of integer type.  A user is free to embellish
this object with other member functions that may, e.g., order subgroups of indices (if they have an inherent order).
Some algorithms that use dynamic_arrays may not function if the ranges are not finite, but there is no
trouble specifying the abstract existance of such an infinite array.
One important point is that this way of doing things does not allow the user to escape at least specifying
the number of array axes, which is the length of the list given as the ranges argument.

Non-trivial advantages are to be had by the layering of rules.  For example, two rules that are implemented
here are caching and masking.  By wrapping the base rule in the cached rule (made available to the user),
previously computed rules are stored.  This could be useful when only some elements of a larger tensor will
be needed repeatedly, but some will not be needed and it is tedious to figure out in advance which 
ones.  The mask function wraps the dynamic_tensor given by the user (without modification) to allow a restriction
of its dimensionality and the ranges of the remaining dimensions (similar to slices ... and might someday be
implemented as such, except that our indices need not be integers).  A final use case for which this was designed
is the wrapping of a base rule in a layer that enforces some symmetry (left to the user to implement).  The ability to give all the rules in
a list (outermost first) is just a convenience to avoid tedious counting of parentheses.  Generally, a rule
wrapper is a function that takes a function (another rule) as an argument and returns a function that takes
indices as an argument and does its thing to either the indices or return value of the base rule and returns the (modified?) value.

It is worth mentioning that nothing prevents the rule from returning mutable object.  This can be a tool (if used with caching)
and a danger.

There is a bit of a dark corner as respects using tuples as indices for *single* tensor axes, such as:
val = mytensor[(4,6),(1,0)]
This could be valueable for packing sparse matrices when not all value combinations are used, or 
for indexing "catalogs" where the length of the tuple is variable even.  However, allowing tuples
as single indices could give ambiguity as to whether tns[1,2] has two indices or one tuple index.  To relieve this issue, the 
forgoing is considered to have two integer indices, and tns[((1,2),)] has one tuple index.  In other words, the user must explicitly pack
tuple indices as 1-tuples if there is only one.  Since this means that a dynamic_array indexed by a single 1-tuple
is read as having one index, then this means that tns[5] and tns[(5,)] need to be semantically the same.  In fact, they
are and one should think if tns[5] as just a convenient short form that works when one only has one index which
itself is not an indexable type.  In fact, "5" is immediately promoted internally to "(5,)".  The same is true with ranges
given for masking; if a single value is given (which is not iteself indexable), then the range is read as consisting of a single
value (a 1-tuple) ... this is not allowed for the initial ranges argument when instantiating an object, because it makes no sense
as that single index would be sliced out, per the behavior of mask, which subsequently ignores axes of length one,
makeing the new tensor appear to be of lower dimensionality.
"""

# Possible improvement:
# We could allow __setitem__ if the rules took an optional final keyword argument, which when passed makes an assignment.
# Rules for which this does not makes sense can simply decide not to accept such an argument and the code will crash if it is attempted.




def _tpl(indices):
    """ A utility that promotes single values to 1-tuples for uniform processing ... means that tuples used as single-axis indices must always be wrapped as 1-tuples explicitly by user"""
    try:
        indices[0]
    except:
        return (indices,)	# it was not a tuple, so make it a 1-tuple
    else:
        return indices		# it was a tuple, so return it (assumes users wrap single-tuple indices inside a 1-tuple)

class _range_not_restricted_class(object):
    """ The default range does not complain no matter what you look for inside it """
    def __init__(self):  pass
    def __contains__(self, item):  return True		# essentially means that anything is "in" this range
    def __str__(self):  return "[unrestricted range]"
    def __len__(self):  raise LogicError("[unrestricted range] does not have finite length")
_range_not_restricted = _range_not_restricted_class()	# singleton


# The main event
class dynamic_array(object):
    """ A class that acts a bit like a function that only takes discrete values on a rectilinear grid, and which is aware of such dimensionality information, and allows for value caching, etc"""
    def __init__(self, rules, ranges):
        self.ranges = [(_range_not_restricted if r is None else r) for r in ranges]	# 'None' is how a user indicates no range restriction on an index
        try:
            rules_rev = list(reversed(rules))	# let rules be a list with the latter rules being inside the former ... so that common modifications (like caching) may be applied to any base rule
        except:
            self.rule = rules			# there was only one rule
        else:
            self.rule = rules_rev[0]
            for rule in rules_rev[1:]:  self.rule = rule(self.rule)
        self.keys_cache = None			# generating all keys could be pretty intense and might be repeated
    def __getitem__(self, indices):
        indices = _tpl(indices)					# promote any single indices to 1-tuples (single indices that *are* tuples must be wrapped as 1-tuples by the user to avoid ambiguity)
        for idx,rng in zip(indices,self.ranges):		# single indices were promoted to 1-tuples
            if idx not in rng:
                print("dyn_array recieved key:", idx)
                print("Value not present in:\n", rng)
                raise KeyError
        return self.rule(*indices)



def wrapper_rule(dyn_arr, sliced_out=None):
    """ Mostly just grabs an element from an underlying array, but 'sliced_out' can indicate some indices of the underlying array that have fixed value """
    if sliced_out is None:  sliced_out = {}
    def get_val(*indices):
        if len(sliced_out)>0:
            dim = len(indices) + len(sliced_out)	# single indices were promoted to 1-tuples
            new_indices = ()
            i = 0
            for j in range(dim):
                if j in sliced_out:
                    new_indices += (sliced_out[j],)
                else:
                    new_indices += (indices[i],)
                    i += 1
            indices = new_indices
        return dyn_arr[indices]
    return get_val

def mask(dyn_arr, ranges):
    """ restricts the values that can be used to index the underlying array (essentially a view), and, if one of those ranges only allows one value, the view behaves as lower dimensional ... can also be used to promote a dict to a dynamic_array"""
    sliced_out = {}
    new_ranges = []
    for i,r in enumerate(ranges):
        if r is None:			# 'None' is how a user indicates no range restriction on an index
            new_ranges += [None]
        else:
            rr = _tpl(r)				# promote any single indices to 1-tuples (single indices that *are* tuples must be wrapped as 1-tuples by the user to avoid ambiguity)
            if len(rr)==1:  sliced_out[i] = rr[0]	# only single possibility given for this axis
            else:           new_ranges += [rr]
    return dynamic_array(wrapper_rule(dyn_arr, sliced_out), ranges=new_ranges)

def wrap(dyn_arr, rules):
    """ allows the creation of a new dynamic_array that effectively has more rules added (outside existing ones) """
    return dynamic_array(rules+[wrapper_rule(dyn_arr)], ranges=dyn_arr.ranges)


class cached(object):
    def __init__(self, rule):
        self.rule = rule
        self.storage = {}
    def __call__(self, *indices):
        if indices not in self.storage:  self.storage[indices] = self.rule(*indices)
        return self.storage[indices]