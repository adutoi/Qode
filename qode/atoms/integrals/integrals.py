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
from . import primitive_integrals

"""\
see comment down near e_repulsion integrals
"""


def overlap(B1,B2):
	r1 = B1.center
	r2 = B2.center
	result = 0
	for g1 in B1.gaussians:
		for g2 in B2.gaussians:
			result += primitive_integrals.overlap(g1,r1, g2,r2)
	return result

class orbital_operator(object):
	pass

class kinetic_class(orbital_operator):
	def __init__(self):
		pass
	def __call__(self, B1, B2):
		r1 = B1.center
		r2 = B2.center
		result = 0
		for g1 in B1.gaussians:
			for g2 in B2.gaussians:
				result += primitive_integrals.kinetic(g1,r1, g2,r2)
		return result
kinetic = kinetic_class()

class charge(orbital_operator):
	def __init__(self, Z, r):
		self.charge = Z
		self.center = r
	def __call__(self, B1, B2):
		r1 = B1.center
		r2 = B2.center
		result = 0
		for g1 in B1.gaussians:
			for g2 in B2.gaussians:
				result += primitive_integrals.charge_op(g1,r1, g2,r2, self.charge,self.center)
		return result

class charges(orbital_operator):
	def __init__(self, charge_ops):
		self.charge_ops = charge_ops
		self.classical_interaction_energy = 0
		for charge_op1 in charge_ops:
			Z1,r1 = charge_op1.charge, charge_op1.center
			for charge_op2 in charge_ops:
				Z2,r2 = charge_op2.charge, charge_op2.center
				if charge_op1 is not charge_op2:
					D,R = primitive_integrals.displacement(r1,r2)
					try:
						self.classical_interaction_energy += Z1 * Z2 / R
					except ZeroDivisionError:
						print("Two charges on top of each other.  Try just making it one summed charge?")
						raise
		self.classical_interaction_energy /= 2
	def __call__(self,B1,B2):
		result = 0
		for charge_op in self.charge_ops:  result += charge_op(B1,B2)
		return result



def e_repulsion(D12,D34):
	B1,B2 = D12
	B3,B4 = D34
	r1 = B1.center
	r2 = B2.center
	r3 = B3.center
	r4 = B4.center
	result = 0
	for g1 in B1.gaussians:
		for g2 in B2.gaussians:
			for g3 in B3.gaussians:
				for g4 in B4.gaussians:
					result += primitive_integrals.e_repulsion(g1,r1, g2,r2, g3,r3, g4,r4)
	return result


"""\
At first, what is below would seem like a bit of a hack, making thes objects look like they were in a "space" when they are not.
For sure, the implementation here is a hack, and the interface is not quite right, but here is where this is going.

In fact, all of the above operators should be "upgraded" to have the behavior below.  In this way, if one just wants integrals
and has no intention of doing transformations, then the space-like notation works anyway.  Of course, space traits is free to
use this notation in its own implementation.  The advantage is the tensor engine which is coming.  outer products of tensors
whose elements are basis functions with 1x1 tensors of operators can be used to build integrals tensors.  Contractions with
tensors that contain numbers can do transformations quickly ("spaces" is for a whole other purpose, not meant for
manipulations of big data).

Remember that the contraction engine is supposed to define different kinds of outer products (basis|basis) is an outer
product where "multiplication" (general bilinear op) is defined as overlap an this builds the S matrix, but basis*basis
builds the pair-distribution basis.
"""



class V_on_pair(object):
	def __init__(self,left_arg,basis):
		self.left_arg = left_arg
		self.basis    = basis
	def __or__(self,right_arg):
		v1,v2 = self.left_arg
		v3,v4 = right_arg
		result = 0
		for c1,B1 in zip(v1.v,self.basis):
			for c2,B2 in zip(v2.v,self.basis):
				for c3,B3 in zip(v3.v,self.basis):
					for c4,B4 in zip(v4.v,self.basis):
						result += c1*c2*c3*c4 * e_repulsion((B1,B2),(B3,B4))
		return result

class repulsion_ints(object):
	def __init__(self,basis):
		self.basis = basis
	def __ror__(self,left_arg):
		return V_on_pair(left_arg,self.basis)
