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
# the next biggest upgrade here is that there is some crossover between when I use an index and when i use a transition
# a transition was really just meant to be a container for an index that could identify the type, but now, I find myself
# needing an index that can identify the type in at least one place.
# a transition has also been upgraded to operating on a
# state, making it "heavier" . . . think this through more clearly.
# It seems like an index should just identify the type, if that is not obvious from some other kind of examination (am I going in circles now?)


from ....util import parallel
from .nested_operator_iterators import iter_single, common, latter


class ZeroCommutator(Exception):
	def __init__(self):  pass
	def __str__(self):  return "You should not be seeing this because it should have been handled."


#def print_nested(encoding, prestring=""):
#	K, E, F, D = range(4)
#	chars = ["K", "E", "F", "D"]
#	if encoding is True:  print(prestring+"...",end="")
#	elif encoding is False:  pass
#	else:
#		if    encoding[K] and prestring=="":  print(chars[K], "  ", end="")
#		elif  encoding[K]:  print(prestring, "  ", end="")
#		for i in [E,F,D]:
#			if encoding[i] is not False:  print_nested(encoding[i], prestring+chars[i])
#	if prestring=="":  print("\n")



# define/use state.add(state1,state2)
def increment(target,delta):
	""" for use in parallel aggregation of state.  no need for a coefficient """
	target.increment(delta)
	return target

def subcall(suboperator,transition,the_state):
	""" need unbound function for parallelization of call """
	return transition(suboperator(the_state))



class nested_operator(object):
	"""\
	Let an operator O be written as:
		O = C + sum_t( G_t * O_t )
	where O_t can be written
		O_t = C_t + sum_u( G_u * O_tu )
	and so forth.

	The set {G_t} is the set of all single particle transitions (D,E,F below, for deexcitation, excitation and flat, respectively).
	C is a scalar, and the scalars C_t build a vector and the scalars C_tu (implicit above) build a matrix, etc.  Therefore, once
	expanded, we have
		O = C + sum_t( C_t * G_t ) + sum_tu( C_tu * G_t * G_u ) + ...
	which builds an arbitrary number-conserving n-body operator.

	The nested form is the form in which this class stores an operator.  This is done by defining operators of the form similar to
	O and O_t above.  Since O and O_t are fundamentally the same in structure, there is a recursion, where by an object of this type contains more
	instances of this type, each implicitly associated with a given G_t.

	Furthermore, the G_t are blocked by their character as
	excitation, deexcitation, and flat (virtual-virtual and occupied-occuied).  Furthermore, in order that every operator (at the highest level) has a 
	unique representation (eliminate redundancy) there are restrictions on indices as one goes deeper into the structure.
	Using commutation rules, it is possible always to insist on the following:
		1. An excitation never occurs to the right of a deexciation or flat operator (i.e. nested inside at a location associated with a deexciation or flat) 
		2. A flat operator never occurs to the right of a deexcitation operator
		3. A flat virtual operator never appears in a string with a flat occupied (they are reorganized to an excitation and a dexcitation)
			The above 3 conditions imply, for example that, if t is a deexcitation, the index u is only over deexciations.
		4. There are four kinds of indices (occ-create, occ-destroy, vrt-create, vrt-destroy), and if two indices of the same kind appear
		   in a string, then the lower index appears to the right of the higher index (lower indices more deeply nested).

	At each layer of recursion, there are two fundamental possibilities.  Either the recursion bottoms out, and the nested operator is just a number,
	which might be zero (i.e. just C, with the blocks for the further G_t absent), or there are some G_t present (any of which might coincidentally be
	associated with zero coefficient at the next layer).  At a futher level of detail, the G_t are present or not in blocks of D, E, or F type transitions.

	"""
	def __init__(self, fluctuation_classes, C=0, Tx_indices=None):
		"""\
		The arguments are as follows:

		fluctuation_classes:  Contains lists of all possible indices in the order which defines "descending", along with reorder and commutator functions reorder, commutator
		C:                    The scalar that multiplies the operator string that *encloses* this entire operator and its sublayers (it does not multiply the sublayers).
		Tx_indices:           A tuple containing 3 lists of (compound) indices for excitation, flat, and deexcitation operators respectively that may be executed at the present layer of recursion.
			              Indices are expected to be immutable types, suitable for keying a dictionary by value, and where equality (==) has an inherent meaning (probably a tuple!), even though ordering is defined by these lists.
		                      If None, default is all possible indices ... normal for top/user layer
		"""
		self._TxClasses = fluctuation_classes
		self._C  = [C]		# wrapped because wrappers will need access to original
		self._Ex = None		# \
		self._Fx = None		#  |- blocks start as None, meaning they are implicitly zero
		self._Dx = None		# /
		if Tx_indices is None:  self._Tx_indices = fluctuation_classes.all_transition_indices()
		else:                   self._Tx_indices = fluctuation_classes.order_indices(Tx_indices)
		self._mask     = None	# Only the iterator and derived classes ever check this
		self._original = self	# for use with destructive masks (since they don't just modify info in containers but overwrite containers in original structure with empty ones)
	def __iter__(self):
		return iter_single(self)
	def __call__(self,the_state,resources=parallel.resources(1)):
		""" Acts the present layer on the state given by recursively calling lower layers and then operating with the enclosing operator """
		# makes assumptions that state member functions .new and .increment are available
		new_state = the_state.new()
		new_state.increment(the_state,self._C[0])
		n_cores = resources.n_cores
		if n_cores==1:
			for suboperator,transition in self:
				new_state.increment(subcall(suboperator,transition,the_state))		# same as:   new_state.increment(transition(suboperator(the_state)))
		else:
			new_state = parallel.aggregate(subcall, (the_state,), increment, n_cores).run(new_state, iter(self))
		return new_state
	def copy(self, _target=None):
		""" copying is masked, but mask is not copied"""
		# _target is for internal recursion only
		if _target is None:  target = nested_operator(self._TxClasses, 0, self._Tx_indices)
		else:                target = _target
		target._C[0] = self._C[0]
		for target_subop,self_subop in latter(target,self):  self_subop.copy(_target=target_subop)		# loops over indices present in latter, getting (implicitly making) subblocks in former as necessary
		return target
	def scale(self,c):
		""" use this also to numerically zero out an operator.  There is no *zero* function because one might interpret that as shrinking storage size, which is only possible if masking is ignored.  use new to get something small """
		self._C[0] *= c
		for suboperator,_ in self:  suboperator.scale(c)
		# intentionally no return to avoid confusion with .scaled()
	def scaled(self,c):
		""" return a scaled copy """
		new = self.copy()
		new.scale(c)
		return new
	def dot(self,other,indent=""):
		product = 0
		primitive = self._C[0] * other._C[0]
		product += primitive
		for subop1,subop2 in common(self,other,indent):  product += subop1.dot(subop2,indent+"  ")			# loops only over indices in common
		return product
	def increment(self,other,c=1):
		""" Add another nested operator to this one """
		self._C[0] += c * other._C[0]
		for subop1,subop2 in latter(self,other):  subop1.increment(subop2,c)			# loops over indices present in latter, getting (implicitly making) subblocks in former as necessary
	def add_term(self,term):
		"""\
		term comes in as a list of operators, presumed to obey the ordering restrictions, which terminates with a scalar the multiplies that string.
		The task is essentially to find the right place in the recursion scheme to place the scalar at the end, which it adds to whatever is there,
		creating blocks along the way if it is needed (works by recursion).
		"""
		try:
			self._C[0] += term[0]
		except:
			self._subtarget(term[0]).add_term(term[1:])
	def dump(self):
		this_dump = [ ((),self._C[0]) ]
		for suboperator,transition in self:
			index = (transition.type, transition.index)
			sub_dump = suboperator.dump()
			for indices,amplitude in sub_dump:  this_dump += [ ( (index,)+indices , amplitude ) ]
		return this_dump
	@staticmethod
	def load(the_dump, fluctuation_classes):
		new = nested_operator(fluctuation_classes)
		for indices,amplitude in the_dump:
			#print(indices,amplitude)
			new.add_term([fluctuation_classes.make_transition(t_type,index) for t_type,index in indices] + [amplitude])
		return new
	def right_multiply(self, transition):
		""" Return a new nested operator that obeys proper index ordering and is equivalent to this operator multipled from the right by a given 1e- number-conserving transition """
		Tx_indices = self._TxClasses.append_indices(self._Tx_indices, transition)
		target = nested_operator(self._TxClasses, 0, Tx_indices)
		target.add_term([transition, self._C[0]])
		for source,inner_trans in self:
			target.increment( source.right_multiply(transition).left_multiply(inner_trans,shallow=True) )
		return target
	def left_multiply(self, transition, shallow=False):
		""" Return a new nested operator that obeys proper index ordering and is equivalent to this operator multipled from the left by a given 1e- number-conserving transition """
		Tx_indices = self._TxClasses.append_indices(self._Tx_indices, transition)
		target = nested_operator(self._TxClasses, 0, Tx_indices)
		target.add_term([transition, self._C[0]])
		for source,inner_trans in self:
			products = self._TxClasses.reorder(transition, inner_trans)
			for pdt in products:
				if shallow:
					if   len(pdt)==1:
						target.increment(source, pdt[0])
					elif len(pdt)==2:
						target._subtarget(pdt[0]).increment(source, pdt[1])
					elif len(pdt)==3:
						target._subtarget(pdt[0])._subtarget(pdt[1]).increment(source, pdt[2])
				else:
					if   len(pdt)==1:
						target.increment(source, pdt[0])
					elif len(pdt)==2:
						target.increment(source.left_multiply(pdt[0]), pdt[1])
					elif len(pdt)==3:
						target._subtarget(pdt[0]).increment(source.left_multiply(pdt[1]), pdt[2])
		return target
	def right_commute(self, transition):
		""" Return a new nested operator that obeys proper index ordering and is equivalent to this operator commuted with a transition on the right """
		Tx_indices = self._TxClasses.append_indices(self._Tx_indices, transition)
		target = nested_operator(self._TxClasses, 0, Tx_indices)
		for suboperator,inner_trans in self:
			target.increment(suboperator.right_commute(transition).left_multiply(inner_trans))
			try:
				commutator = self._TxClasses.commute(inner_trans, transition)
			except ZeroCommutator:
				pass
			else:
				for term in commutator:
					if isinstance(term,int):  target.increment(suboperator.scaled(term))
					else:                     target.increment(suboperator.left_multiply(term))
		return target
	def _subtarget(self, transition):
		index = transition.index
		if   transition.is_excitation:
			if index not in self._Tx_indices.Excitations:       raise Exception("operator string does not obey index restrictions")
			block = self._get_Ex_block()
		elif transition.is_flat:
			if index not in self._Tx_indices.Flat_transitions:  raise Exception("operator string does not obey index restrictions")
			block = self._get_Fx_block()
		elif transition.is_deexcitation:
			if index not in self._Tx_indices.Deexcitations:     raise Exception("operator string does not obey index restrictions")
			block = self._get_Dx_block()
		return block[index]
	def _get_Ex_block(self):
		""" Return the block of excitations, whether or not it previously existed """
		if self._Ex is None:
			self._Ex = {}
			for index in self._Tx_indices.Excitations:
				transition = self._TxClasses.excitation(index)      
				Tx_indices = self._TxClasses.slice_indices(self._Tx_indices, transition)
				self._Ex[index] = nested_operator(self._TxClasses, 0, Tx_indices)
		return self._Ex
	def _get_Fx_block(self):
		""" Return the block of flat transitions, whether or not it previously existed """
		if self._Fx is None:
			self._Fx = {}
			for index in self._Tx_indices.Flat_transitions:
				transition = self._TxClasses.flat_transition(index)      
				Tx_indices = self._TxClasses.slice_indices(self._Tx_indices, transition)
				self._Fx[index] = nested_operator(self._TxClasses, 0, Tx_indices)
		return self._Fx
	def _get_Dx_block(self):
		""" Return the block of deexcitations, whether or not it previously existed """
		if self._Dx is None:
			self._Dx = {}
			for index in self._Tx_indices.Deexcitations:
				transition = self._TxClasses.deexcitation(index)      
				Tx_indices = self._TxClasses.slice_indices(self._Tx_indices, transition)
				self._Dx[index] = nested_operator(self._TxClasses, 0, Tx_indices)
		return self._Dx



class masked(nested_operator):
	""" mostly just a shallow copy of the parent, but the mask will be passed through the iterator """
	def __init__(self, the_mask, parent):
		self._TxClasses = parent._TxClasses
		self._C  = parent._C  if the_mask.C_on  else [0]
		self._Ex = parent._Ex if the_mask.Ex_on else None
		self._Fx = parent._Fx if the_mask.Fx_on else None
		self._Dx = parent._Dx if the_mask.Dx_on else None
		self._Tx_indices = parent._Tx_indices
		self._mask     = the_mask
		self._original = parent._original
		if the_mask.is_destructive:  self._destroy_rest()
	def _destroy_rest(self):
		raise Exception("This untested code should not be being called yet")
		if not self._mask.C_on:   self._original._C  = [0]
		if not self._mask.Ex_on:  self._original._Ex = None
		if not self._mask.Fx_on:  self._original._Fx = None
		if not self._mask.Dx_on:  self._original._Dx = None
		for suboperator,transition in self:  pass		# The simple act of recursive looping over what is not destroyed applies destructive masks to everything left

def no_mask(operator):  return operator

def mask_it(primitive_mask, destructive=False, inverted=False):
	if (not inverted and primitive_mask is True) or (inverted and primitive_mask is False):
		return no_mask
	else:
		return mask(primitive_mask,destructive,inverted)

class mask(object):
	def __init__(self, primitive_mask, destructive=False, inverted=False):		# primitive_mask = [True, [True, [False, [False, True, False, False], False, False], False, False], [True, False, False, False], False], for example.
		K, E, F, D = range(4)
		self.is_destructive = destructive
		if (not inverted and primitive_mask is True) or (inverted and primitive_mask is False):
			raise Exception("What are you doing?  This is not a mask.")
		elif (not inverted and primitive_mask is False) or (inverted and primitive_mask is True):
			self.C_on  = False	# \
			self.Ex_on = False	#  |- Since these are all "off" the submask should never be asked for
			self.Fx_on = False	#  |
			self.Dx_on = False	# /
		else:
			if inverted:
				raise Exception("This untested code should not be being called yet")
				self.C_on = (not primitive_mask[K])
				if primitive_mask[E] is True or primitive_mask[E] is False:
					self.Ex_on   = (not primitive_mask[E])
					self.Ex_mask = mask_it((not primitive_mask[E]), destructive)
				if primitive_mask[F] is True or primitive_mask[F] is False:
					self.Fx_on   = (not primitive_mask[F])
					self.Fx_mask = mask_it((not primitive_mask[F]), destructive)
				if primitive_mask[D] is True or primitive_mask[D] is False:
					self.Dx_on   = (not primitive_mask[D])
					self.Dx_mask = mask_it((not primitive_mask[D]), destructive)
			else:
				self.C_on  =      primitive_mask[K]		# True or False directly
				self.Ex_on = bool(primitive_mask[E])		# non-empty lists test as True in python (oddly ==True does not work here though)
				self.Fx_on = bool(primitive_mask[F])		# non-empty lists test as True in python
				self.Dx_on = bool(primitive_mask[D])		# non-empty lists test as True in python
				self.Ex_mask = mask(primitive_mask[E], destructive)
				self.Fx_mask = mask(primitive_mask[F], destructive)
				self.Dx_mask = mask(primitive_mask[D], destructive)
		self._arguments = (primitive_mask, destructive, inverted)
	def __call__(self, operator):
		return masked(self, operator)
	def inverted(self):
		primitive_mask, destructive, inverted = self._arguments
		return mask(primitive_mask, destructive, not inverted)
	def destructive(self):
		primitive_mask, destructive, inverted = self._arguments
		destructive = True
		return mask(primitive_mask, destructive, inverted)

zero_mask = mask(False, destructive=True)



# Unbound functions for purposes of parallelization
def mult_func_right(subop2,transition,Op1):  return multiply(Op1.right_multiply(transition),subop2,True)
def mult_func_left( subop1,transition,Op2):  return multiply(subop1,Op2,False).left_multiply(transition)

class aggregator(object):
	""" This aggregator enforces the correct final structure onto the target """
	def __init__(self, target):
		self._arguments = (target._TxClasses, 0, target._Tx_indices)
	def __call__(self, op1, op2):
		target = nested_operator(*self._arguments)
		target.increment(op1)
		target.increment(op2)
		return target

def multiply(Op1, Op2, use_right_multiply=False, resources=parallel.resources(1)):
	if Op1._TxClasses!=Op2._TxClasses:  raise Exception("illegal multiplication of different operator types")
	fluctuation_classes = Op1._TxClasses
	Tx_indices = fluctuation_classes.combine_indices(Op1._Tx_indices, Op2._Tx_indices)
	target = nested_operator(fluctuation_classes, 0, Tx_indices)
	#
	n_cores = resources.n_cores
	if use_right_multiply:
		target.increment(Op1, Op2._C[0])
		if n_cores==1:					# essential so that parallelization only at top layer
			for suboperator,transition in Op2:
				target.increment(mult_func_right(suboperator,transition,Op1))
		else:
			target = parallel.aggregate(mult_func_right, (Op1,), aggregator(target), n_cores).run(target, iter(Op2))
	else:
		target.increment(Op2, Op1._C[0])
		if n_cores==1:					# essential so that parallelization only at top layer
			for suboperator,transition in Op1:
				target.increment(mult_func_left(suboperator,transition,Op2))
		else:
			target = parallel.aggregate(mult_func_left, (Op2,), aggregator(target), n_cores).run(target, iter(Op1))
	return target



def comm_func_1(subop2,transition,Op1):  return multiply(Op1.right_commute(transition),subop2)
def comm_func_2(subop2,transition,Op1):  return commute(Op1,subop2).left_multiply(transition)

# Currently, this will fail unless Op2 contains only excitations
def commute(Op1, Op2, resources=parallel.resources(1)):
	if Op1._TxClasses!=Op2._TxClasses:  raise Exception("illegal commutation of different operator types")
	fluctuation_classes = Op1._TxClasses
	Tx_indices = fluctuation_classes.combine_indices(Op1._Tx_indices, Op2._Tx_indices)
	target = nested_operator(fluctuation_classes, 0, Tx_indices)
	#
	n_cores = resources.n_cores
	if n_cores==1:					# essential so that parallelization only at top layer
		for suboperator,transition in Op2:
			target.increment(comm_func_1(suboperator,transition,Op1))
			target.increment(comm_func_2(suboperator,transition,Op1))
	else:
		target = parallel.aggregate(comm_func_1, (Op1,), aggregator(target), n_cores).run(target, iter(Op2))
		target = parallel.aggregate(comm_func_2, (Op1,), aggregator(target), n_cores).run(target, iter(Op2))
	return target
