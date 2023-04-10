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
  "H" : contraction(
    exponents    = [3.42525091, 0.62391373, 0.16885540],
    coefficients = [0.15432897, 0.53532814, 0.44463454]
    ),
  "He" : contraction(
    exponents    = [6.36242139, 1.15892300, 0.31364979],
    coefficients = [0.15432897, 0.53532814, 0.44463454]
    )
}
