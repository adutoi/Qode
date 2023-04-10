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
from math import sqrt, exp, acos, atan2
import spherical_harmonics

Y = spherical_harmonics.wf()

def pol_coord(x,y,z):
	r = sqrt(x**2 + y**2 + z**2)
	th = acos(z/r)
	ph = atan2(y,x)
	return (r,th,ph)

def R_10(r,Q,a):  return  2 * ((Q/a)**(3/2.)) * exp(-Q*r/a)

def R_21(r,Q,a):  return  (2/sqrt(3)) * ((Q/(2*a))**(5/2.)) * r * exp(-Q*r/(2*a))

def wf1S(x,y,z,Q,a):
	r,th,ph = pol_coord(x,y,z)
	return R_10(r,Q,a) * Y['0 0'](th,ph)

def wf2Pz(x,y,z,Q,a):
	r,th,ph = pol_coord(x,y,z)
	return R_21(r,Q,a) * Y['1 0'](th,ph)

Hatom_eigen = {}
Hatom_eigen['1s']   = wf1S
Hatom_eigen['2p_z'] = wf2Pz

def eigen(i,Q,a):
	def wavefunction(x,y,z):  return Hatom_eigen[i](x,y,z,Q,a)
	return wavefunction

def wf(Q,a):
	wavefunctions = {}
	for i in ['1s','2p_z']:  wavefunctions[i] = eigen(i,Q,a)
	return wavefunctions
