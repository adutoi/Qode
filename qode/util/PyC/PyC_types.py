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

# Here we instantiate some python objects that share the same names as
# some C typedefs in the sister file PyC_types.h, in order that the
# information that they convey may be used to create and communicate
# information that is compatible with C routines.
#
# Below we document the names of some numpy types that are defined to be
# equivalent to certain types in the local C implementation and also give
# name that ctypes uses for this.  Other safe types could be, added but
# these are the most useful to use now.
#
# An important theoretical aspect is that C and numpy really specify the
# data type itself.  The ctypes name only specifies how existing data is
# passed (I think).
#
# C type          ctypes name          numpy dtype          PyC name    
# ------          -----------          -----------          --------
# int             c_int                intc                 Int
# ssize_t         c_ssize_t            intp                 BigInt
# double          c_double             float64*             Double
#   * In theory this one could break(?), but only on weird systems
#
# The source of the information here is largely from these two pages:
# https://docs.scipy.org/doc/numpy-1.13.0/user/basics.types.html
# https://docs.python.org/3/library/ctypes.html

import ctypes
import numpy

class PyC_type(object):
    def __init__(self, numpy_dtype, ctypes_name):
        self.numpy  = numpy_dtype
        self.ctypes = ctypes_name

# These lines attach numpy and ctypes naming to our PyC name.
# The PyC_types.h connects this information to the C types.

Int    = PyC_type(numpy.intc,    ctypes.c_int)
BigInt = PyC_type(numpy.intp,    ctypes.c_ssize_t)
Double = PyC_type(numpy.float64, ctypes.c_double)

# Just a list of all the low-level types that are mapped

PyCtypes = [Int, BigInt, Double]

# Associates native python types with PyCtypes defined here.
# Unlike the PyCtype<->ctypes<->numpy translations, which are guaranteed safe (but low-level) specifications,
# the suitability of this mapping depends on the user, but it is likely fine for most. """
# Should I make this somehow dynamic (user defined)? ... but how to make sure it gets communicated to C code without being too clunky?

py_types = { int: BigInt, float: Double }
