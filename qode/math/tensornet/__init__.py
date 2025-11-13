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

from .base           import evaluate, increment, raw, scalar_value, shape, initialize_timer, print_timings
from .tensors        import tensor_sum, backend_contract_path    # tensor_sum() can initialize an empty accumulator for += use
from .tensors        import primitive_tensor as _primitive_tensor
from .backends       import dummy_backend, numpy_backend, tensorly_backend

print("""
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
""")

# Philosopy:  A backend (in the sense of the value of the argument backend below) is a class that
# contains static methods that perform basic operations on objects that represent tensors of a specific
# type. These will henceforth be called backend tensors.  Tensornet then uses the backend functions to
# manipulate the backend tensors according to what is encoded by its high-level syntax.  Objects of the
# class defined below are factories for the primitive tensornet tensors which each wrap a backend tensor
# (and never modifies it).  An important point is that a given backend can have multiple factories
# depending on how what a user is doing.
#
# The creation of a tensornet primitive tensor should only take place via a factory (which is why
# that class name is "hidden" here with the underscore convention) because this is the mechanism by
# which a backend tensor is associated with its backend functions and the utility for contraction creation.
# In the most basic case, this creation is done by the .init() method (not .__init__, which creates the 
# factory).  The method .zeros() is a convenience function for creation of a tensornet primitive storing
# a new backend tensor (via a backend method) so the user need not create such a zero tensor first and
# then "wrap" it (.scalar_tensor() is similar and was added for completion, but commenting it out essentially
# documents that we have not yet found a use case ... the necessity for the backend to be able to create
# such a thing however has to do with obeying expectations about type uniformity).
#
# Finally we come to the arguments data_wrap and data_unwrap, the method .data(), and the idea of multiple
# factories for the same backend.  It is important that tensornet works only with a notion of backend
# tensors and the functions it is given for them without knowing anything further about their properties.
# But the user may be working with some abstracted backend tensors.  Although they contain, or
# somehow represent the existence of, a specific concrete tensor (eg, numpy, tensorly, pytorch), they 
# do not have the usual expected behavior (e.g., [] might not even get an element).  We will call such 
# things meta tensors, it might be very inconvenient for the user to make their in-house meta tensor behave 
# exactly like a concrete tensor, right down to letting a numpy.ndarray (etc) be instantiated from it.  It
# will also be very inconvenient if the user has to always use the unbound function raw() and then afterward
# extract the concrete data for cooperation with other linear algebra libraries.  The converse is also true
# that it will be annoyingly repetitive to always create the meta tensor from a concrete tensor from elsewhere,
# and then pass that to the factory.  So, therefore, a factory can be created that takes the functions
# data_wrap and data_unwrap.  Although their actions are not constrained as such, the idea is that data_wrap
# is a function that takes a concrete tensor and returns an object of the meta tensor type with which
# backend works, and data_unwrap does the opposite, extracting and returning (by reference) the concrete
# tensor represented (eg, contained within) the meta tensor.  The convenience function .data() called
# on a tensornet tensor causes evaluation of the (potential?) network (sum?) and returns directly the concrete
# tensor associated with the result.  By default, data_wrap and data_unwrap are just pass-throughs,
# which means that .data() is just an alias for the unbound function raw(), which forces evaluation and
# returns a backend tensor (which is typically already a concrete tensor type).  Since this automatic
# wrapping and unwrapping might not always be the desired behavior, more than one factory can be defined,
# per need.
#
# Addendum 1:  There is no requirement that a meta tensor is not itself also a concrete type, and that
# this concrete type might be initialized from multiple other concrete types, which might not even be
# different than its own type.  So if we interpret the concrete and meta layers to both be tensorly,
# where the meta layer can be initialized also from numpy or pytorch (etc), then there it is a perfectly
# clean for the user to provide as data_wrap a function that makes a tensorly tensor from some input
# which may or may not be of tensorly type), and to leave data_unwrap as the default pass-through.
#
# Addendum2:  One question that can be asked about the design is: why am I giving in a separated
# backend with functions that perform manipulations on a tensor, rather than just insisting that backend
# tensors have certain properties, like += (which would also make the insides of tensornet prettier).
# The answer is that one might have enough serious work to do in the abstracted/meta backend already,
# without having to define yet more operations that make the thing look almost like a tensor, but not
# completely.  One can imagine various inheritance and virtual method schemes, but one is always going
# to end up in one of two places: 1. Forcing the user to make a new class for their backend tensors
# (and then maybe always, even when using "normal" tensors), or 2. automating the creation of the 
# backend tensors with the expected operators in terms of some other/divorced set of instructions.
# In case 2, we are back where we started from the user point of view, though the insides of tensornet
# look nicer (for the price of another file/class and layer of abstraction).
#
# Addendum 3: The .__init__() method used to take a copy_data argument which meant that the factory should
# copy the provided primitive before storing it as an argument for further use.  Probably unusual, but
# a user might do this if they intended on modifying that primitive for some other use while it was still
# in use by the tensornet code (you never know ...).  This can now be done by the user with the data_wrap
# option.

class primitive_tensor_factory(object):
    def __init__(self, backend, data_wrap=lambda d: d, data_unwrap=lambda d: d):
        self._backend     = backend
        self._data_wrap   = data_wrap
        self._data_unwrap = data_unwrap
    def init(self, data):
        raw_tensor = self._data_wrap(data)
        return _primitive_tensor(raw_tensor, self._backend)    # injects 'contract' into tensor_base so it they can contract "themselves"
    def zeros(self, shape):
        return _primitive_tensor.zeros(shape, self._backend)
    #def scalar_tensor(self, scalar):
    #    return _primitive_tensor.scalar_tensor(scalar, self._backend)
    def data(tensor):
        return self._data_unwrap(raw(tensor))

dummy_tensor = primitive_tensor_factory(dummy_backend.functions)
np_tensor    = primitive_tensor_factory(numpy_backend.functions)
tl_tensor    = primitive_tensor_factory(tensorly_backend.functions)
