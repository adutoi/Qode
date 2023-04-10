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

import numpy
import ctypes
from .PyC_types import PyCtypes, py_types



def native2PyCtype(obj):
    """ Associates native python types with PyCtypes, according to PyC_types.py_types """
    for pytype,PyCtype in py_types.items():
        if obj is pytype or isinstance(obj,pytype):  return PyCtype
    return None    # returned if obj is not (and not instance of) one of the supported native types

def dtype2PyCtype(array):
    """ Associates numpy data types with PyCtypes, according to PyC_types.PyCtypes """
    PyCtype = None
    for test in PyCtypes:
        if array.dtype==test.numpy:  PyCtype = test
    if PyCtype is None:  raise AssertionError
    return PyCtype

def numpy_Cptr(array, PyCtype):
    """ Returns a ctypes pointer of the specified type for the data in the given numpy array """
    return array.ctypes.data_as(ctypes.POINTER(PyCtype.ctypes))



def C_arg(obj):
    """ Returns an object suitable for passing to a library loaded with ctypes.cdll.LoadLibrary """
    # This will convert native python types according to the py_types dictionary above, and it will give the ctypes data pointers for numpy arrays.
    # It also understands simple (non-nested) indexable objects (e.g. lists) of native types and numpy arrays and will return a C pointer to an array (of pointers).
    # There is beauty in simplicity, and this should cover the majority of usage cases (simple numbers and arrays of simple numbers and pointers),
    # hiding the low level from the user.  If the user wants more, they should be able to translate their array of pointers to a flat list, etc.
    # (or build another function on top of this one).
    PyCtype = native2PyCtype(obj)
    if PyCtype:
        return PyCtype.ctypes(obj)
    elif isinstance(obj, numpy.ndarray):
        PyCtype = dtype2PyCtype(obj)
        return numpy_Cptr(obj, PyCtype)
    else:
        # assume it is a homogeneous indexable object of one of the above supported types
        try:     obj0 = obj[0]
        except:  raise AssertionError
        PyCtype = native2PyCtype(obj0)
        if PyCtype:
            type0 = type(obj0)
            array = (PyCtype.ctypes * len(obj))()
            for i,element in enumerate(obj):
                if type(element) is not type0:  raise AssertionError
                array[i] = PyCtype.ctypes(element)
            return array
        elif isinstance(obj0, numpy.ndarray):
            PyCtype = dtype2PyCtype(obj0)
            array = (ctypes.POINTER(PyCtype.ctypes) * len(obj))()
            for i,np_array in enumerate(obj):
                if not isinstance(np_array, numpy.ndarray) or np_array.dtype!=PyCtype.numpy:  raise AssertionError
                array[i] = numpy_Cptr(np_array, PyCtype)
            return array
        else:
            raise AssertionError

def C_args(*py_args):
    """ Converts multiple objects in one function call to types suitable for passing to a library loaded with ctypes.cdll.LoadLibrary """
    return [C_arg(py_arg) for py_arg in py_args]



def pythonize(C_function, return_pytype=None):
    """ Given a function in a library loaded with ctypes.cdll.LoadLibrary, wrap and hide the python/numpy<->ctypes conversions of arguments and the return value """
    # The signature of the C function must use the types conversion specified PyC_types.py/.h and in the py_types dictionary above
    PyCtype = native2PyCtype(return_pytype)
    if PyCtype:  C_function.restype = PyCtype.ctypes
    def py_function(*py_args):
        return C_function(*C_args(*py_args))
    return py_function
