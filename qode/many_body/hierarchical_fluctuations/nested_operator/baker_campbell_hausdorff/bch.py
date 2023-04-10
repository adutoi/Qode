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
import pickle 	# for debugging
from .BCHtruncation import BCHtruncation
from ..             import      operator, zero_mask, multiply, commute



def _commute(X, T, resources):
	""" evaluate X = [X,T] using wisdom about not computing things that will never matter, keeping excitations """
	direct = False	# direct is slower !?
	if direct:
		return commute(X,T,resources)
	else:
		commutator =         multiply(X, T, resources=resources)
		commutator.increment(multiply(T, X, resources=resources), -1)
		return commutator

def _partition_and_strip(target, X, mask, dump):
	""" strip X down to excitations and parts that matter going forward and *move* excitations from source to target (incrementing target and zeroing out excitations in source) """
	dump, order = dump
	if dump:
		to_dump = mask.to_keep(X).copy()
		to_dump.increment(mask.excitations(X))
		pickle.dump(to_dump.dump(),open("X_{}.pickled".format(order),"wb"))
	#
	new_X = mask.to_keep(X).copy()
	target.increment(mask.excitations(X))
	debug = False
	if debug:
		mask.excitations(X).scale(0)
		mask.to_keep(X).scale(0)
		mask.to_trash(X).scale(0)
		check = X.dot(X)
		if check>1e-10:  raise Exception("Somewhere in the BCH understanding on which this is based, something went wrong.  Error = {}".format(check))
	return new_X

def _BCH_term(to_increment, order, to_be_commuted, wisdom, resources, textlog, dump):
	X,T  = to_be_commuted
	mask = wisdom.masks(order)
	#
	if order==0:
		textlog("... order 0 commutator ...")
	else:
		textlog("... order {} commutator ... running on {} cores ...".format(order, resources.n_cores))
		X = _commute(X, T, resources)
		X.scale(1./order)	# factors accumulate, generating factorial
	#
	X = _partition_and_strip(to_increment, X, mask, (dump,order))
	to_be_commuted[0] = X		# in place modification of list element accessible to caller



class BakerCampbellHausdorff(object):
	""" This class is a wrapper for the nested-operator-specific dirty-work.  It has an interface needed by the coupled_cluster module, but hides the internal nature of the operator.  """
	def __init__(self, fluctuation_classes, resources, H_order=2, CC_order=2):
		self.fluctuations = fluctuation_classes
		self.wisdom       = BCHtruncation(CC_order, H_order)
		self.resources    = resources
		self.n_calls      = 0			# here for debugging
	def computeOmega(self, H, T, textlog):	# call must have this signature!   Computes excitation part (0th, 1st and 2nd order) of HBar
		self.n_calls += 1
		dump = (self.n_calls==4)
		dump = False		# manual override.  Comment out this line if you want this to dump a cycle for debugging
		if dump:
			print("DUMPING OPERATORS! ... Did you remember to set direct=True in _commute? ... it's slower but produces cleaner dumps")
			pickle.dump(H.dump(),open("H.pickled","wb"))
			pickle.dump(T.dump(),open("T.pickled","wb"))
		#
		Omega = operator(self.fluctuations)				# initialized to zero
		X = H.copy()							# zeroth-order commutator is H, but copy because will be modified in place
		commutator = [X,T]						# X in this list will be modified in place, just a happy accident that it looks like commutator notation
		for n in range(self.wisdom.max_commutator+1):			# Assumes than n-body Hamiltonian includes n-th-order deexcitations (otherwise can stop sooner)
			_BCH_term(Omega, n, commutator, self.wisdom, self.resources, textlog, dump)
		if dump:  pickle.dump(Omega.dump(),open("Omega.pickled","wb"))
		return Omega
	@staticmethod
	def Energy(Omega, textlog):		# call must have this signature!   Retrieves energy from excitation part (0th order component) of HBar
		return Omega._C[0]
