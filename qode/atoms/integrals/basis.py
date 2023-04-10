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
import copy
from math import *
import numpy
from ..static_data.Gaussian_basis_sets import data
from . import primitive_integrals

"""\
Right now, this is really over simplified.  It has no concept of having more than one function per atom even.

The grand plan would be, using strings as either dictionary keys or function arguments, to have those dictionaries
or functions return functions that take as their sole argument a position vector, which then returns a completed basis
function from a standard set like STO-3G.

That said, it should be easy to construct custom basis functions.  First, think about what constructors of basis functions
would look like if there were no such things as standard sets, and then think about how to smoothly congjure up functions 
from standard sets.
"""



def Snorm(exponent):
	return (2*exponent/pi)**(3./4)

contractions = {}
for basis_set,atoms in data.items():
	new = {}
	for atom,raw in atoms.items():
		new[atom] = [primitive_integrals.gaussian(c*Snorm(a),a) for c,a in zip(raw.coefficients,raw.exponents)]
	contractions[basis_set] = new



class basis_function(object):
	def __init__(self,g,r):
		self.gaussians = g	# list of Gaussians (including prefactors) to sum up ...
		self.center    = r	# ... that all share the same center
	def __or__(self,other):
		if isinstance(other,basis_function):
			B1 = other
			B2 = self
			r1 = B1.center
			r2 = B2.center
			result = 0
			for g1 in B1.gaussians:
				for g2 in B2.gaussians:
					result += primitive_integrals.overlap(g1,r1, g2,r2)
			return result
		else:
			return other.__ror__(self)
	def _copy(self):
		return copy.deepcopy(self)
	def _scale(self,c):
		for g in self.gaussians:  g.prefactor *= c
	def __mul__(self,c):
		new = self._copy()
		new._scale(c)
		return new
	def __rmul__(self,c):
		return self * c		# calls __mul__
	def __imul__(self,c):
		self._scale(c)
		return self
	def __truediv__(self,c):
		return self * (1./c)	# calls __mul__
	def __itruediv(self,c):
		self *= (1./c)		# calls __imul__ above
		return self



def symmON_transformation(basis_functions):
	""" returns a numpy 2D array where the columns are coefficients of orthonormal orbitals in the input basis """
	n_bas = len(basis_functions)

	overlaps = numpy.zeros((n_bas,n_bas))
	for m,Bm in enumerate(basis_functions):
		for n,Bn in enumerate(basis_functions):
			overlaps[m,n] = (Bm|Bn)

	evals, evecs = numpy.linalg.eigh(overlaps)

	invrtD_evecsT = numpy.zeros((n_bas,n_bas))
	for row in range(n_bas):
		for column in range(n_bas):
			invrtD_evecsT[row,column] = evecs[column,row] / sqrt(evals[row])

	U = numpy.tensordot(evecs,invrtD_evecsT, (1,0))

	return U


def transform2(T_0, U):
	T_1 = numpy.tensordot(U, T_0, (0,1))
	T_2 = numpy.tensordot(U, T_1, (0,1))
	return T_2

def transform4(T_0, U):
	T_1 = numpy.tensordot(U, T_0, (0,3))
	T_2 = numpy.tensordot(U, T_1, (0,3))
	T_3 = numpy.tensordot(U, T_2, (0,3))
	T_4 = numpy.tensordot(U, T_3, (0,3))
	return T_4




