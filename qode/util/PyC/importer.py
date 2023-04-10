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

import os
import inspect
from   .loader import load_C
from   .args   import pythonize

# Usage, assuming you have a file named code.c with a function func in it that returns an int
#
# from qode.util.PyC import import_C
# fast_module = import_C("code", ... other args ...)
#
# # Now this
# result = fast_module.func.return_type(int)( ... arguments ...)
#
# # Or this
# fast_module.func.return_type(int)
# result1 = fast_module.func( ... arguments1 ...)
# result2 = fast_module.func( ... arguments2 ...)
#
# # Or this
# fast = fast_module.func.return_type(int)
# result1 = fast( ... arguments1 ...)
# result2 = fast( ... arguments2 ...)

# The biggest weakness here is that this does not handle connections between multiple C-source files ...
# ... until I think out the best way to do that, we will just #include and recompile code we want to reuse.
# Actually, if I ever need anything more than this (which has the advantage of being super lightweight,
# automatically present, and conforming to the idea of "as little C as possible"), then I should look
# into real tools for this job.  For a good discussion of Cython, SWIG, etc, see:
# https://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html 

# Also, as of now there is no recompile if the arguments to the import (such as compilation) change.
# This can be fixed by hashing the arguments to a string that gets stored in a file and checked.



class _function_wrapper(object):
    def __init__(self, C_function):
        self.C_function    = C_function
        self.function      = None
        self.return_pytype = None
    def __call__(self, *args):
        if self.function is None:  self.function = pythonize(self.C_function, self.return_pytype)
        return self.function(*args)
    def return_type(self, pytype):
        if self.function is None:		# Can "change our mind" up until the time it is first called ... after that, it *really* makes no sense
            self.return_pytype = pytype
            return self
        elif pytype is self.return_pytype:
            return self
        else:
            raise "Cannot change return type of C function after function has been called"


class import_C(object):
    def __init__(self, C_filestem, flags="", include=None, directory=None, cc=None, ld=None):
        if directory is None:  directory = os.path.dirname(os.path.abspath(inspect.stack()[1][1]))   # assume .c code is in the same directory as the calling .py module.  See:  https://stackoverflow.com/questions/13699283/how-to-get-the-callers-filename-method-name-in-python
        self._C_module = load_C(C_filestem, flags, include, directory, cc, ld)
        self._wrapped = {}
    def __getattr__(self, name):
        if name not in self._wrapped:  self._wrapped[name] = _function_wrapper(getattr(self._C_module, name))
        return self._wrapped[name]
