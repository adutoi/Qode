#    (C) Copyright 2018, 2019 Yuhong Liu and Anthony D. Dutoi
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
from ....util.PyC import import_C, Double, BigInt
from ....math     import space
hamiltonian    = import_C("hamiltonian",    include=["partition_storage.h"], flags="-O3 -lm")
preconditioner = import_C("preconditioner", include=["partition_storage.h"], flags="-O3")



# These three function could be coded more compactly in a recursive manner
def _initialize_0():
	return numpy.zeros(1, dtype=Double.numpy)	# just a single numer accessible by "reference"
def _initialize_1(Nfrag, Nvrt, transitions):
	transitions, = transitions
	size_tot = 0
	for M in range(Nfrag):
		size = 1
		if transitions[0]=="X":  size  = Nvrt[M]
		if transitions[1]=="X":  size *= Nvrt[M]
		size_tot += size
	return numpy.zeros(size_tot, dtype=Double.numpy)
def _initialize_2(Nfrag, Nvrt, transitions):
	transitions1, transitions2 = transitions
	size_tot = 0
	for M1 in range(Nfrag):
		size1 = 1
		if transitions1[0]=="X":  size1  = Nvrt[M1]
		if transitions1[1]=="X":  size1 *= Nvrt[M1]
		for M2 in range(Nfrag):
			size2 = 1
			if transitions2[0]=="X":  size2  = Nvrt[M2]
			if transitions2[1]=="X":  size2 *= Nvrt[M2]
			if M1!=M2:  size_tot += size1*size2	# Do not allocate space for unreferenced diagonal blocks
	return numpy.zeros(size_tot, dtype=Double.numpy)

_block_dims = {\
"K":     [], \
"Ex":    [("O","X")], \
"Fo":    [("O","O")], \
"Fv":    [("X","X")], \
"Dx":    [("X","O")], \
"ExEx":  [("O","X"),("O","X")], \
"ExFo":  [("O","X"),("O","O")], \
"ExFv":  [("O","X"),("X","X")], \
"ExDx":  [("O","X"),("X","O")], \
"FoFo":  [("O","O"),("O","O")], \
"FoFv":  [("O","O"),("X","X")], \
"FoDx":  [("O","O"),("X","O")], \
"FvFv":  [("X","X"),("X","X")], \
"FvDx":  [("X","X"),("X","O")], \
"DxDx":  [("X","O"),("X","O")]  \
}

class operator(object):
	"""operator class holding H, T, or Omega"""
	def __init__(self):	# not meant to be called from the outside world
		self.blocks = {}
	def scale(self, factor):
		if factor!=1:
			for block in self.blocks:  self.blocks[block] *= factor
	def dot(self, other):
		result = 0
		for block in self.blocks:
			factor = 1
			if block in ("ExEx", "FoFo", "FvFv", "DxDx"):  factor = 1/2.	# takes care of redundancy due to symmetry (in the general case, since factor is a factorial, could test if 'ExEx' (etc) is in string and multiply by 1/2 (etc)
			if block in other.blocks:
				result += factor * numpy.dot(self.blocks[block], other.blocks[block])
		return result
	def add(self, increment, scalar=1):
		for block in increment.blocks:
			if block in self.blocks:
				self.blocks[block] += scalar * increment.blocks[block]
			else:
				self.blocks[block]  = numpy.copy(increment.blocks[block])
				if scalar!=1:  self.blocks[block] *= scalar
	def copy(self):
		the_copy = operator()
		the_copy.add(self)
		return the_copy
	@staticmethod
	def initialize(states_per_frag, blocks):
		Nfrag = len(states_per_frag)
		Nvrt  = [n-1 for n in states_per_frag]
		op = operator()
		for block in blocks:	# when _initialize functions are recusively coded, there are no 'if's and len(...) is an argument
			if len(_block_dims[block])==0:
				op.blocks[block] = _initialize_0()
			if len(_block_dims[block])==1:
				op.blocks[block] = _initialize_1(Nfrag, Nvrt, _block_dims[block])
			if len(_block_dims[block])==2:
				op.blocks[block] = _initialize_2(Nfrag, Nvrt, _block_dims[block])
		return op



def SD_operator(states_per_frag):
	""" Sets up a cluster operator with single and double substitutions """
	return operator.initialize(states_per_frag, ("Ex", "ExEx"))

def RSD_operator(states_per_frag):
	""" Sets up a cluster operator with single and double substitutions, and an reference amplitude """
	return operator.initialize(states_per_frag, ("K", "Ex", "ExEx"))

def H_012_operator(states_per_frag, monomer_Hamiltonians, dimer_Couplings, thresh=1e-12):
	""" Sets up XR-CC Hamiltonian with up to dimer couplings """
	dummy = numpy.zeros(1, dtype=Double.numpy)	# need an actual pointer so as to not confuse C interface
	dimer_Couplings_flatlist = []
	for row in dimer_Couplings:
		for matrix in row:  dimer_Couplings_flatlist += [dummy if matrix is None else matrix]
	Nfrag = len(states_per_frag)
	Nvrt  = numpy.array([n-1 for n in states_per_frag], dtype=BigInt.numpy)
	H = operator.initialize(states_per_frag, ("K", "Ex", "Fo", "Fv", "Dx", "ExEx", "ExFo", "ExFv", "ExDx", "FoFo", "FoFv", "FoDx", "FvFv", "FvDx", "DxDx"))
	Hb = H.blocks
	hamiltonian.load(Nfrag, Nvrt, thresh, monomer_Hamiltonians, dimer_Couplings_flatlist, Hb["K"], Hb["Ex"], Hb["Fo"], Hb["Fv"], Hb["Dx"], Hb["ExEx"], Hb["ExFo"], Hb["ExFv"], Hb["ExDx"], Hb["FoFo"], Hb["FoFv"], Hb["FoDx"], Hb["FvFv"], Hb["FvDx"], Hb["DxDx"])
	return Hb["K"][0], H	# Return constant part of H separately as well, since it is ignored in the BCH implementation



class inv_diagE(object):
	"""Class for preconditioning Omega by diving through by energy differences """
	def __init__(self, diag_energies, states_per_frag):
		self.diag_energies   = [numpy.array(Em,dtype=Double.numpy) for Em in diag_energies]
		self.states_per_frag = states_per_frag
	def __call__(self, omega):
		Nfrag = len(self.states_per_frag)
		Nvrt  = numpy.array([n-1 for n in self.states_per_frag], dtype=BigInt.numpy)
		preconditioned_omega = omega.copy()
		preconditioned_omega.blocks["K"][0] = 0
		Ob = preconditioned_omega.blocks
		preconditioner.divide_by_diagE(Nfrag, Nvrt, Ob["Ex"], Ob["ExEx"], self.diag_energies)
		return preconditioned_omega



def _dot(v,w):          return v.dot(w)
def _add_to(v,w,n=1):   v.add(w,n)
def _scale(n,v):        v.scale(n)
def _copy(v):           return v.copy()
def _act_on_vec(op,v):  return op(v)	# This is a single tasker, only here for acting energy denominators
def _check_member(v):   pass		# skip for now (but still needs to be implemented to avoid a crash)
def _check_lin_op(op):  return False	# skip for now (but still needs to be implemented to avoid a crash)
space_traits = space.traits(field=Double.numpy, dot=_dot, add_to=_add_to, scale=_scale, copy=_copy, act_on_vec=_act_on_vec, check_member=_check_member, check_lin_op=_check_lin_op)
