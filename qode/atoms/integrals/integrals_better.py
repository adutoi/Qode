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
import numpy
from . import primitive_integrals



def resolve_arrays(integral, arguments, _processed=[]):
	""" each of the arguments might be single number or an array, interpret as an outer product of inputs and return array of outputs with dimensionality equal to number of array inputs """
	if not arguments:
		return integral._contracted_integral(*_processed)
	else:
		if isinstance(arguments[0],list):
			result = []
			for arg in arguments[0]:
				new_arguments = [arg] + arguments[1:]
				result += [resolve_arrays(integral, new_arguments, _processed)]
			return result
		else:
			new_processed = _processed + [arguments[0]]
			return resolve_arrays(integral, arguments[1:], new_processed)



class bra_e_repulsion(object):
	def __init__(self,left_arg):
		self.left_arg = left_arg
	def __or__(self,right_arg):
		B1,B2 = self.left_arg
		B3,B4 = right_arg
		return numpy.array(resolve_arrays(self, [B1,B2,B3,B4]))
	def _contracted_integral(self,B1,B2,B3,B4):
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

class e_repulsion(object):
	def __init__(self):
		pass
	def __ror__(self,left_arg):
		return bra_e_repulsion(left_arg)



class bra_kinetic(object):
	def __init__(self,left_arg):
		self.left_arg = left_arg
	def __or__(self,right_arg):
		return numpy.array(resolve_arrays(self, [self.left_arg, right_arg]))
	def _contracted_integral(self,B1,B2):
		r1 = B1.center
		r2 = B2.center
		result = 0
		for g1 in B1.gaussians:
			for g2 in B2.gaussians:
				result += primitive_integrals.kinetic(g1,r1, g2,r2)
		return result

class kinetic(object):
	def __init__(self):
		pass
	def __ror__(self,left_arg):
		return bra_kinetic(left_arg)



class bra_charge(object):
	def __init__(self, Z, r, left_arg):
		self.charge = Z
		self.center = r
		self.left_arg = left_arg
	def __or__(self,right_arg):
		return numpy.array(resolve_arrays(self, [self.left_arg, right_arg]))
	def _contracted_integral(self,B1,B2):
		r1 = B1.center
		r2 = B2.center
		result = 0
		for g1 in B1.gaussians:
			for g2 in B2.gaussians:
				result += primitive_integrals.charge_op(g1,r1, g2,r2, self.charge, self.center)
		return result

class charge(object):
	def __init__(self, Z, r):
		self.charge = Z
		self.center = r
	def __ror__(self,left_arg):
		return bra_charge(self.charge, self.center, left_arg)



class bra_charges(object):
	def __init__(self, charge_ops, left_arg):
		self.charge_ops = charge_ops
		self.left_arg = left_arg
	def __or__(self,right_arg):
		return numpy.array(resolve_arrays(self, [self.left_arg, right_arg]))
	def _contracted_integral(self,B1,B2):
		result = 0
		for charge_op in self.charge_ops:  result += (B1|charge_op|B2)
		return result

class charges(object):
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
	def __ror__(self,left_arg):
		return bra_charges(self.charge_ops, left_arg)

# overlap is implemented in basis.py
