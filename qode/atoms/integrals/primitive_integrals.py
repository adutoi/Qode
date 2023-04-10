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
from math import *

"""\
About the most primitive integrals engine one can imagine.  It only does Gaussian s orbitals.
"""



class gaussian(object):
	def __init__(self,prefactor,exponent):
		self.prefactor = float(prefactor)
		self.exponent  = float(exponent)	# subject to interpretation by calling code (is it normalized before multiplication by prefactor?)



#
# Special mathematical primitives
#

# The combined exponent of a gaussian product
def alpha(g1,g2):
	return g1.exponent + g2.exponent

# The reduced exponent of a gaussian product (multiplies the distance dependence of the norm)
def mu(g1,g2):
	return g1.exponent * g2.exponent / alpha(g1,g2)

# The displacement vector D=r2-r1, and its norm
def displacement(r1,r2):
	x1,y1,z1 = r1
	x2,y2,z2 = r2
	Dx,Dy,Dz = x2-x1, y2-y1, z2-z1
	R = sqrt(Dx**2 + Dy**2 + Dz**2)
	return (Dx,Dy,Dz), R

# If a gaussian at the origin is multiplied with another centered a distance +R away, on a line, this gives
# the exponent and center, c, of the resulting distribution.  The prefactor is for a distribution
# which has a classical norm (not squared norm) of unity.
def chg_dist_line(g1,g2,R):
	a = alpha(g1,g2)
	m = mu(g1,g2)
	s = g1.prefactor * g2.prefactor * exp(-m*(R**2)) * sqrt((pi/a)**3)
	c = R * g2.exponent / a
	return ( gaussian(s,a) , c )

# If a gaussian at 3-vector position r1 is multiplied with another centered at r2, this gives
# the exponent and center, r, of the resulting distribution.  The prefactor is for a distribution
# which has a classical norm (not squared norm) of unity.
def chg_dist(g1,r1, g2,r2):
	(Dx,Dy,Dz), R = displacement(r1,r2)
	g,c = chg_dist_line(g1,g2,R)
	if R==0:
		r = r1
	else:
		x1,y1,z1 = r1
		r = x1+Dx*(c/R) , y1+Dy*(c/R), z1+Dz*(c/R)
	return g,r



#
# Coulomb integral kernel
#

# This can be thought of as just erf(sqrt(a)*r)/r, but with safe evaluation near r=0
def repulsion_primitive(a,r):
	ar = sqrt(a) * r
	answer = 0
	if ar<0.1:
		T = ar*ar
		term = 1
		for i in range(1000):
			answer = answer + term
			term = term * (2*T) / (2*(i+1)+1)
		answer = 2 * sqrt(1/pi) * exp(-T) * answer
	else:
		answer = erf(ar) / ar
	return sqrt(a) * answer



#
# The integrals.
# Gaussians g1,g2,g3,g4 have no implicit norm [= prefactor*exp(-exponent*x*x)]
# r1,r2,r3,r4 are their vector centers.
#

def overlap(g1,r1, g2,r2):
	D,R = displacement(r1,r2)
	a = alpha(g1,g2)
	m = mu(g1,g2)
	return g1.prefactor * g2.prefactor * sqrt((pi/a)**3) * exp(-m*(R**2))

def kinetic(g1,r1, g2,r2):
	D,R = displacement(r1,r2)
	m = mu(g1,g2)
	return ( (1./2) * m * (6 - 4*m*(R**2)) * overlap(g1,r1,g2,r2) )

def charge_op(g1,r1, g2,r2, Z,rZ):	# Z need not be an integer (or even positive)
	g,r = chg_dist(g1,r1, g2,r2)	# classicaly normalized, then with a prefactor!
	D,R = displacement(r,rZ)
	return ( -Z * g.prefactor * repulsion_primitive(g.exponent,R) )

def e_repulsion(g1,r1, g2,r2, g3,r3, g4,r4):
	g12,r12 = chg_dist(g1,r1, g2,r2)	# classicaly normalized, then with a prefactor!
	g34,r34 = chg_dist(g3,r3, g4,r4)	# classicaly normalized, then with a prefactor!
	D,R = displacement(r12,r34)
	a = 1./( 1./g12.exponent + 1./g34.exponent )
	return ( g12.prefactor * g34.prefactor * repulsion_primitive(a,R) )
