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
from math import sqrt, pi, sin, cos

def Y_00(th,ph):  return sqrt(1/(4*pi))

def Y_10(th,ph):  return sqrt(3/(4*pi)) * cos(th)

rot_eigen = {}
rot_eigen['0 0'] = Y_00
rot_eigen['1 0'] = Y_10

def wf():  return rot_eigen
