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
from .contraction import contraction

data = {\
  "H 1s" : contraction(
    exponents    = [5.4471780, 0.8245470],
    coefficients = [0.1562850, 0.9046910]
    ),
  "H 2s" : contraction(
    exponents    = [0.1831920],
    coefficients = [1.0000000]
    ),
  "He 1s" : contraction(
    exponents    = [13.6267000, 1.9993500],
    coefficients = [ 0.1752300, 0.8934830]
    ),
  "He 2s" : contraction(
    exponents    = [0.3829930],
    coefficients = [1.0000000]
    )
}
