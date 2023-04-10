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
import builtins
import math
import cmath
import numpy

metric = {\
  float            : float         \
, complex          : float         \
, numpy.float64    : numpy.float64 \
, numpy.complex128 : numpy.float64 \
}

# even though the ints do not build a field, some of these are defined to work on them because we should not restrict what can be done

def abs(number):
    """\
    The values in metric above are defined by the return type of this function.
    I don't include integer classes, because they are not fields, and, outside of this module,
    this abs is qualified by the field_traits namespace.
    """
    if   isinstance(number,(int,float,complex)):               return builtins.abs(number)
    elif isinstance(number,(numpy.float64,numpy.complex128)):  return    numpy.abs(number)
    else:  raise error("field_traits.abs executed on unavailable type")

def conjugate(number):
    if   isinstance(number,(int,float,numpy.float64)):   return number
    elif isinstance(number,(complex,numpy.complex128)):  return number.conjugate()
    else:  raise error("field_traits.conjugate executed on unavailable type")

def exp(number):
    if   isinstance(number,(int,float)):                       return  math.exp(number)
    elif isinstance(number,complex):                           return cmath.exp(number)
    elif isinstance(number,(numpy.float64,numpy.complex128)):  return numpy.exp(number)
    else:  raise error("field_traits.exp executed on unavailable type")

def sqrt(number):
    if   isinstance(number,(int,float)):                       return  math.sqrt(number)
    elif isinstance(number,complex):                           return cmath.sqrt(number)
    elif isinstance(number,(numpy.float64,numpy.complex128)):  return numpy.sqrt(number)
    else:  raise error("field_traits.sqrt executed on unavailable type")

# Keep going with more functions as needed ...
