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
from math import sqrt, exp, pi

def wf0(x,a):  return ((a/pi)**(1/4.)) * exp(-(a/2.)*(x**2))

def wf1(x,a):  return  sqrt(2*a)   *  x * wf0(x,a)

def wf2(x,a):  return  sqrt(1/2.)  * (2*a*(x**2)    - 1) * wf0(x,a)

def wf3(x,a):  return  sqrt(3*a)   * (2*a*(x**3)/3. - x) * wf0(x,a)

def wf4(x,a):  return  sqrt(1/6.)  * (2*(a**2)*(x**4) -  6*a*(x**2) +    3/2.) * wf0(x,a)

def wf5(x,a):  return  sqrt(a/15.) * (2*(a**2)*(x**5) - 10*a*(x**3) + 15*x/2.) * wf0(x,a)
