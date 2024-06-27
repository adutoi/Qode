#  (C) Copyright 2020, 2024 Anthony D. Dutoi
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

# This is one step more explicit than passing around tuples, where you have to keep track of
# the order, but has a lighter syntax than dictionaries.  Variants of unpacking, splatting,
# and double-splatting are supported, also for subsets of the stored values.  In principle,
# for unpacking/splatting, you could rely on the dict order, but it is a bad idea.
#
# The basic  usage syntax is:
# >>> a = struct(attr1=arg1, attr2=arg2, attr3=arg3, ...)
# >>> a.later = another
# >>> (a.attr1 is arg1) and (a.later is another)    # True
# >>> x, y = as_tuple(a("attr1", "attr2"))    # user determines order, can be a subset
# >>> x, y = as_tuple(a("attr1 attr2"))       # same thing
# >>> f(*as_tuple(a("attr1 attr2")))          # trivial use- ...
# >>> f(**as_dict(a("attr1 attr2")))          # ... case extension
# >>> f(**as_dict(a))                         # This is completely safe/predictable if want all members
#
# Another use is for one-off classes where field names (accessible via dot syntax) are predetermined:
# >>> class data(struct):
# >>>     def __init__(self, attr1, attr2):
# >>>         struct.__init__(self, attr1=attr1, attr2=attr2)
# >>>     def __str__(self):
# >>>         return f"My contents are: {self.attr1}, {self.attr2}"
#
# What is different than collections.namedtuple?
# - instances are mutable (can be advantage)
# - unpacking/splatting syntax is heavier (disadvantage) but more flexible (advantage)
# - double splatting exists (and light-weight syntax if everything extracted)
#
# The names (attr1, attr2, etc) can be any valid identifier name (like what would be valid
# for the name of a variable, function, class, etc), but should not begin with an underscore
# (not enforced though, and you can get a way with it for all but a small number of cases).
# See implementation notes at the end of the file (to understand some of the terse comments)

import itertools



def as_tuple(to_tuple):
    return to_tuple._as_tuple()

def as_dict(to_dict):
    return to_dict._as_dict()

class struct(object):
    def __init__(self, **kwargs):                         # implicitly restricts allowed keys of internal dict
        object.__setattr__(self, "_data_dict", kwargs)    # avoids circular recursion
    def __setattr__(self, attr, value):
        self._data_dict[attr] = value
        return value
    def __getattr__(self, attr):
        try:
            _data_dict = object.__getattribute__(self, "_data_dict")    # avoids circular recursion
        except:
            return object.__getattribute__(self, attr)
        else:
            try:
                return _data_dict[attr]
            except KeyError:
                raise AttributeError(repr(attr))
    def __delattr__(self, attr):
        try:
            del self._data_dict[attr]
        except KeyError:
            raise AttributeError(repr(attr))
    def __call__(self, *attrs):    # allows the taking of a sub-namespace
        attrs = [attr for attr in itertools.chain.from_iterable([attr.split(" ") for attr in attrs]) if attr!=""]
        kwargs = {attr:self._data_dict[attr] for attr in attrs}
        return struct(**kwargs)
    def __repr__(self):
        arguments = ", ".join("{}={}".format(k,repr(v)) for k,v in self._data_dict.items())
        return "{}({})".format(type(self).__name__, arguments)
    def _as_tuple(self):
        return tuple(self._data_dict.values())
    def _as_dict(self):
        return dict(self._data_dict)



# Implementation notes:
#
# Some basic theory:  First, in general, __getattr__, if defined, is only called if dot syntax is
# otherwise unresolved, whereas __getattribute__, if defined, will take over any dot-syntax call.
# Defining both __getattr__ and __getattribute__ makes no sense.  The built-in object class only
# has __getattribute__ defined.  Conversely, __setattr__ and __setattribute__ cannot meaningfully
# be different from each other, and so only __setattr__ is ever defined.  Second, when class
# instances (henceforth, obj or objs, to distinguish them from the ojbect class) are pickled, the
# class definitions (ie code for member functions) are not stored.  The fully qualified class name
# (ie, including module) is stored, and the modules containing the definitions must be importable
# later (this is so that old pickles can be loaded with updated code).  Only the data for the obj is 
# stored, and this is restored upon unpickling after a new obj of that class has been created using
# the __new__ function (which, if not defined by the user, is taken from the object base class).
# For the most part we do not want to mess with __new__ and all the other machinery that lives down
# at that level, but it is important to know that we can, because this can be modified, and other
# special __xxx__ functions can be defined for a class, to define how an obj of that class is 
# unpickled, and we need to understand that such functions might be called during unpickling.
# When an obj is unpickled, its __init__ function is not called!  This makes sense because the obj
# is not being initialized, pe se, but rather restored.  Even more detail is found here:
# https://docs.python.org/3.10/library/pickle.html#what-can-be-pickled-and-unpickled
#
# Now the application of the above:  First, we want an object that is blank in a sense, having no
# names reserved and unavailable to label user data (that do not begin with an underscore, at least).
# We also do not want labels for user data co-mingling with the administrative members of the object.
# So we put the user data in a separate data dictionary (._data_dict member), which is itself one of
# the administrative members.  Since we do not know what attributes the user might want to set, we
# use __setattr__ to redirect *all* dot-syntax assignments to the data dictionary.  But this means
# that we lose any ability to assign administrative members through the dot-syntax mechanism; just
# trying to assign "self._data_dict" would lead to a circular recursion.  This is circumvented by
# using the __setattr__ method of the parent object class.  Under normal use, it would not be
# necessary to use the methods of the object base class to *get* an attribute because __getattr__
# does not redirect *all* dot-syntax recalls to the the data dictionary.  Assuming __init__ has
# been run, then the data dictionary will exist and the dot-syntax for it will be resolved.  However,
# if the obj is being reconstructed from a pickle, __init__ is not run and it can be in such a state
# that ._data_dict does not yet exist and optional functions from the class that we did not implement
# (for guiding the unpickling) are being requested using the dot syntax.  These would get redirected
# to the data dictionary, which does not yet exist, leading to a circular recursion.  Therefore,
# we first used the __getattribute__ method of the parent object class to test if ._data_dict exists
# if not, we use the same object.__getattribute__ function to access the member that is being 
# requested, which presumably raises AttributeError, causing the unpickling mechanism to move on
# (if the attribute did exist, __getattr__ would never be called).  Even if ._data_dict does exist,
# as long as it does not contain a key with one of these special function names, it will anyway
# raise AttributeError, when the key is not found.  So, double leading underscores are unsafe for
# data labels, but too lazy to enforce it now.
#
# as_tuple and as_dict are implemented only as external functions (accessing leading-underscore
# member functions) just to adhere to the philosophy that any labels that do not begin with
# underscores are fair game.
