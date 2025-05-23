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
from copy import deepcopy
import queue
from ..util                     import parallel
from .state                     import state, dot, create, annihilate, op_string
from .reorder                   import excitation, deexcitation, occ_rearrange, vrt_rearrange, reorder_t_t, commute_t_e, ZeroCommutator
from .nested_operator_iterators import iter_single, common, latter

def print_nested(encoding, prestring=""):
	K, E, F, D = range(4)
	chars = ["K", "E", "F", "D"]
	if encoding is True:  print(prestring+"...",end="")
	elif encoding is False:  pass
	else:
		if    encoding[K] and prestring=="":  print(chars[K], "  ", end="")
		elif  encoding[K]:  print(prestring, "  ", end="")
		for i in [E,F,D]:
			if encoding[i] is not False:  print_nested(encoding[i], prestring+chars[i])
	if prestring=="":  print("\n")



# define/use state.add
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
	def __init__(self,occ_strings,CrtOcc_range,DstOcc_range,CrtVrt_range,DstVrt_range,C=0):
		"""\
		The arguments are as follows:

		occ_strings:	This is tantamount to specifying the space on which this operator is operating in the context
				of the calling code.  Probably this shoudl be deprecated in some way, with the state being operated on
				communicating this.  But for now it stands.
		CrtOcc_range:	A list of all occupied indices on which a creation operator may be executed at the present layer of recursion
				(so excluding higher occupied indices on which a creation is executed in the enclosing layers).  This list must
				be in descending order, but the indices need not form a contiguous block.
		DstOcc_range:	Like CrtOcc_range, but for destruction of an occupied
		CrtVrt_range:	Like CrtOcc_range, but for creation of a virtual
		DstVrt_range:	Like CrtOcc_range, but for destuction of a virtual
				^ Note: relying on user that Occ and Vrt sets are disjoint at highest level, regardless of whether operator otherwise restricted somehow
		C:		The scalar that multiplies the operator string that encloses this entire operator and its sublayes.
		"""
		self._C  = [C]		# wrapped because wrappers will need access to original
		self._Ex = None
		self._Dx = None
		self._Fo = None
		self._Fv = None
		# The sortings are necessary because reorder functions work on raw index values and not on the order they occur in the list ... unless one looks inside, you'd never know I did this though
		self._CrtOcc_range = sorted(CrtOcc_range, reverse=True)
		self._CrtVrt_range = sorted(CrtVrt_range, reverse=True)
		self._DstOcc_range = sorted(DstOcc_range, reverse=True)
		self._DstVrt_range = sorted(DstVrt_range, reverse=True)
		self._occ_strings  = occ_strings
		self._mask     = None		# Only the iterator and derived classes ever check this
		self._original = self		# for use with destructive masks (since they don't just modify info in containers but overwrite containers in original structure with empty ones)
	def __iter__(self):
		return iter_single(self)
	def __call__(self,the_state,resources=parallel.resources(1)):
		""" Acts the present layer on the state given by recursively calling lower layers and then operating with the enclosing operator """
		new_state = state(self._occ_strings)
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
		if _target is None:  target = nested_operator(self._occ_strings, self._CrtOcc_range, self._DstOcc_range, self._CrtVrt_range, self._DstVrt_range)
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
		#print("   ",primitive)
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
	def right_multiply(self, transition):
		""" Return a new nested operator that obeys proper index ordering and is equivalent to this operator multipled from the right by a given 1e- number-conserving transition """
		CrtOcc_range = sorted(list( set(self._CrtOcc_range) | set(transition.CrtOcc) ),reverse=True)
		DstOcc_range = sorted(list( set(self._DstOcc_range) | set(transition.DstOcc) ),reverse=True)
		CrtVrt_range = sorted(list( set(self._CrtVrt_range) | set(transition.CrtVrt) ),reverse=True)
		DstVrt_range = sorted(list( set(self._DstVrt_range) | set(transition.DstVrt) ),reverse=True)
		target = nested_operator(self._occ_strings, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
		target.add_term([transition, self._C[0]])
		for source,inner_trans in self:
			target.increment( source.right_multiply(transition).left_multiply(inner_trans,shallow=True) )
		return target
	def left_multiply(self, transition, shallow=False):
		""" Return a new nested operator that obeys proper index ordering and is equivalent to this operator multipled from the left by a given 1e- number-conserving transition """
		CrtOcc_range = sorted(list( set(self._CrtOcc_range) | set(transition.CrtOcc) ),reverse=True)
		DstOcc_range = sorted(list( set(self._DstOcc_range) | set(transition.DstOcc) ),reverse=True)
		CrtVrt_range = sorted(list( set(self._CrtVrt_range) | set(transition.CrtVrt) ),reverse=True)
		DstVrt_range = sorted(list( set(self._DstVrt_range) | set(transition.DstVrt) ),reverse=True)
		target = nested_operator(self._occ_strings, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
		target.add_term([transition, self._C[0]])
		for source,inner_trans in self:
			products = reorder_t_t(transition, inner_trans)
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
	def right_commute(self, Ex_trans):
		""" Return a new nested operator that obeys proper index ordering and is equivalent to this operator commuted with an excitation on the right """
		if not Ex_trans.is_excitation:  raise Exception("sorry, only commutation with excitations on the right is currently supported")
		CrtOcc_range = sorted(list( set(self._CrtOcc_range) | set(Ex_trans.CrtOcc) ),reverse=True)
		DstOcc_range = sorted(list( set(self._DstOcc_range) | set(Ex_trans.DstOcc) ),reverse=True)
		CrtVrt_range = sorted(list( set(self._CrtVrt_range) | set(Ex_trans.CrtVrt) ),reverse=True)
		DstVrt_range = sorted(list( set(self._DstVrt_range) | set(Ex_trans.DstVrt) ),reverse=True)
		target = nested_operator(self._occ_strings, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
		for suboperator,transition in self:
			target.increment(suboperator.right_commute(Ex_trans).left_multiply(transition))
			try:
				commutator = commute_t_e(transition, Ex_trans)
			except ZeroCommutator:
				pass
			else:
				for term in commutator:
					if isinstance(term,int):  target.increment(suboperator.scaled(term))
					else:                     target.increment(suboperator.left_multiply(term))
		return target
	def _subtarget(self, transition):
		if   transition.is_excitation:
			a,i = transition.E_transition
			if a not in self._CrtVrt_range or i not in self._DstOcc_range:  raise Exception("operator string does not obey index restrictions")
			return self._get_Ex_block()[a][i]
		elif transition.is_deexcitation:
			a,i = transition.D_transition
			if a not in self._DstVrt_range or i not in self._CrtOcc_range:  raise Exception("operator string does not obey index restrictions")
			return self._get_Dx_block()[a][i]
		elif transition.is_occ_rearrange:
			i,j = transition.Fo_transition
			if i not in self._DstOcc_range or j not in self._CrtOcc_range:  raise Exception("operator string does not obey index restrictions")
			return self._get_Fo_block()[i][j]
		elif transition.is_vrt_rearrange:
			a,b = transition.Fv_transition
			if a not in self._CrtVrt_range or b not in self._DstVrt_range:  raise Exception("operator string does not obey index restrictions")
			return self._get_Fv_block()[a][b]
	def _get_Ex_block(self):
		""" Return the block of excitations, whether or not it previously existed """
		if self._Ex is not None:  return self._Ex
		else:
			self._Ex = {}					# Use dictionaries because might not be contiguous.  A better version will used np.array and an index map.
			for na,a in enumerate(self._CrtVrt_range):
				row = {}
				for ni,i in enumerate(self._DstOcc_range):  row[i] = nested_operator(self._occ_strings,self._CrtOcc_range,self._DstOcc_range[ni+1:],self._CrtVrt_range[na+1:],self._DstVrt_range)
				self._Ex[a] = row
			return self._Ex
	def _get_Dx_block(self):
		""" Return the block of deexcitations, whether or not it previously existed """
		if self._Dx is not None:  return self._Dx
		else:
			self._Dx = {}					# Use dictionaries because might not be contiguous.  A better version will used np.array and an index map.
			for na,a in enumerate(self._DstVrt_range):
				row = {}
				for ni,i in enumerate(self._CrtOcc_range):  row[i] = nested_operator(self._occ_strings,self._CrtOcc_range[ni+1:],[],[],self._DstVrt_range[na+1:])
				self._Dx[a] = row
			return self._Dx
	def _get_Fo_block(self):
		""" Return the block of flat occupied transitions, whether or not it previously existed """
		if self._Fo is not None:  return self._Fo
		else:
			self._Fo = {}					# Use dictionaries because might not be contiguous.  A better version will used np.array and an index map.
			for ni,i in enumerate(self._DstOcc_range):
				row = {}
				for nj,j in enumerate(self._CrtOcc_range):  row[j] = nested_operator(self._occ_strings,self._CrtOcc_range[nj+1:],self._DstOcc_range[ni+1:],[],self._DstVrt_range)
				self._Fo[i] = row
			return self._Fo
	def _get_Fv_block(self):
		""" Return the block of flat virtual transitions, whether or not it previously existed """
		if self._Fv is not None:  return self._Fv
		else:
			self._Fv = {}					# Use dictionaries because might not be contiguous.  A better version will used np.array and an index map.
			for na,a in enumerate(self._CrtVrt_range):
				row = {}
				for nb,b in enumerate(self._DstVrt_range):  row[b] = nested_operator(self._occ_strings,self._CrtOcc_range,[],self._CrtVrt_range[na+1:],self._DstVrt_range[nb+1:])
				self._Fv[a] = row
			return self._Fv



class masked(nested_operator):
	""" mostly just a shallow copy of the parent, but the mask will be passed through the iterator """
	def __init__(self, the_mask, parent):
		#print("my mom")
		#print("C:{} Ex:{} Dx:{} Fo:{} Fv:{}".format(parent._C, parent._Ex, parent._Dx, parent._Fo, parent._Fv))
		#print(the_mask.C_on)
		self._C  = parent._C  if the_mask.C_on  else [0]
		self._Ex = parent._Ex if the_mask.Ex_on else None
		self._Dx = parent._Dx if the_mask.Dx_on else None
		self._Fo = parent._Fo if the_mask.Fo_on else None
		self._Fv = parent._Fv if the_mask.Fv_on else None
		#print("me")
		#print("C:{} Ex:{} Dx:{} Fo:{} Fv:{}".format(self._C, self._Ex, self._Dx, self._Fo, self._Fv))
		self._CrtOcc_range = parent._CrtOcc_range
		self._CrtVrt_range = parent._CrtVrt_range
		self._DstOcc_range = parent._DstOcc_range
		self._DstVrt_range = parent._DstVrt_range
		self._occ_strings  = parent._occ_strings
		self._mask     = the_mask
		self._original = parent._original
		if the_mask.is_destructive:  self._destroy_rest()
	def _destroy_rest(self):
		raise Exception("This untested code should not be being called yet")
		if not self._mask.C_on:   self._original._C  = [0]
		if not self._mask.Ex_on:  self._original._Ex = None
		if not self._mask.Dx_on:  self._original._Dx = None
		if not self._mask.Fo_on:  self._original._Fo = None
		if not self._mask.Fv_on:  self._original._Fv = None
		for suboperator,transition in self:  pass		# The simple act of recursive looping over what is not destroyed applies destructive masks to everything left

def no_mask(operator):  return operator

def mask_it(primitive_mask, destructive=False, inverted=False):
	if (not inverted and primitive_mask is True) or (inverted and primitive_mask is False):
		#print("block not masked")
		return no_mask
	else:
		return mask(primitive_mask,destructive,inverted)

class mask(object):
	def __init__(self, primitive_mask, destructive=False, inverted=False):		# primitive_mask = [True, [True, [False, [False, True, False, False], False, False], False, False], [True, False, False, False], False], for example.
		K, E, F, D = range(4)
		#print("making mask:  ", end="")
		#print_nested(primitive_mask)
		#print(primitive_mask)
		self.is_destructive = destructive
		if (not inverted and primitive_mask is True) or (inverted and primitive_mask is False):
			raise Exception("What are you doing?  This is not a mask.")
		elif (not inverted and primitive_mask is False) or (inverted and primitive_mask is True):
			#print("hit False")
			self.C_on  = False	#
			self.Ex_on = False	#
			self.Fo_on = False	# Since these are all "off" the submask should never be asked for
			self.Fv_on = False	#
			self.Dx_on = False	#
		else:
			if inverted:
				raise Exception("This untested code should not be being called yet")
				self.C_on = (not primitive_mask[K])
				if primitive_mask[E] is True or primitive_mask[E] is False:
					self.Ex_on   = (not primitive_mask[E])
					self.Ex_mask = mask_it((not primitive_mask[E]), destructive)
				if primitive_mask[D] is True or primitive_mask[D] is False:
					self.Dx_on   = (not primitive_mask[D])
					self.Dx_mask = mask_it((not primitive_mask[D]), destructive)
				if primitive_mask[F] is True or primitive_mask[F] is False:
					self.Fo_on   = (not primitive_mask[F])
					self.Fv_on   = (not primitive_mask[F])
					self.Fo_mask = mask_it((not primitive_mask[F]), destructive)
					self.Fv_mask = mask_it((not primitive_mask[F]), destructive)
			else:
				#print("parsing block ...")
				#print("deep: primitive_mask[K]=",primitive_mask[K])
				self.C_on  =      primitive_mask[K]		# True or False directly
				#print("???", self.C_on)
				self.Ex_on = bool(primitive_mask[E])		# non-empty lists test as True in python (oddly ==True does not work here though)
				self.Fo_on = bool(primitive_mask[F])		# non-empty lists test as True in python
				self.Fv_on = bool(primitive_mask[F])		# non-empty lists test as True in python
				self.Dx_on = bool(primitive_mask[D])		# non-empty lists test as True in python
				#print("ON??   C:{} Ex:{} Fo:{} Fv:{} Dx:{}".format(self.C_on, self.Ex_on, self.Fo_on, self.Fv_on, self.Dx_on))
				self.Ex_mask = mask(primitive_mask[E], destructive)
				self.Fo_mask = mask(primitive_mask[F], destructive)
				self.Fv_mask = mask(primitive_mask[F], destructive)
				self.Dx_mask = mask(primitive_mask[D], destructive)
				#print("done parsing block")
		self._arguments = (primitive_mask, destructive, inverted)
	def __call__(self, operator):
		#print("???", self.C_on)
		return masked(self, operator)
	def inverted(self):
		primitive_mask, destructive, inverted = self._arguments
		return mask(primitive_mask, destructive, not inverted)
	def destructive(self):
		primitive_mask, destructive, inverted = self._arguments
		destructive = True
		return mask(primitive_mask, destructive, inverted)


zero_mask = mask(False, destructive=True)



def RSD_excitations(reference, the_state, occ_strings, occ_orbs, vrt_orbs):
	occ_orbs = sorted(occ_orbs, reverse=True)
	vrt_orbs = sorted(vrt_orbs, reverse=True)
	X = nested_operator(occ_strings,[],occ_orbs,vrt_orbs,[])
	X.add_term([dot(reference,the_state)])
	for a in vrt_orbs:
		for i in occ_orbs:
			amplitude = dot(op_string(create(a),annihilate(i))(reference),the_state)
			X.add_term([excitation(a,i),amplitude])
	for idxa,a in enumerate(vrt_orbs):
		for idxi,i in enumerate(occ_orbs):
			for b in vrt_orbs[idxa+1:]:
				for j in occ_orbs[idxi+1:]:
					amplitude = dot(op_string(create(a),annihilate(i),create(b),annihilate(j))(reference),the_state)
					X.add_term([excitation(a,i),excitation(b,j),amplitude])
	return X



# Unbound functions for purposes of parallelization
def mult_func_right(subop2,transition,Op1):  return multiply(Op1.right_multiply(transition),subop2,True)
def mult_func_left( subop1,transition,Op2):  return multiply(subop1,Op2,False).left_multiply(transition)

class aggregator(object):
	""" This aggregator enforces the correct final structure onto the target """
	def __init__(self, target):
		occ_strings  = target._occ_strings
		CrtOcc_range = target._CrtOcc_range
		DstOcc_range = target._DstOcc_range
		CrtVrt_range = target._CrtVrt_range
		DstVrt_range = target._DstVrt_range
		self._arguments = (occ_strings, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)	
	def __call__(self, op1, op2):
		target = nested_operator(*self._arguments)
		target.increment(op1)
		target.increment(op2)
		return target

def multiply(Op1, Op2, use_right_multiply=False, resources=parallel.resources(1)):
	CrtOcc_range = sorted(list( set(Op1._CrtOcc_range) | set(Op2._CrtOcc_range) ),reverse=True)
	DstOcc_range = sorted(list( set(Op1._DstOcc_range) | set(Op2._DstOcc_range) ),reverse=True)
	CrtVrt_range = sorted(list( set(Op1._CrtVrt_range) | set(Op2._CrtVrt_range) ),reverse=True)
	DstVrt_range = sorted(list( set(Op1._DstVrt_range) | set(Op2._DstVrt_range) ),reverse=True)
	target = nested_operator(Op1._occ_strings, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
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
	CrtOcc_range = sorted(list( set(Op1._CrtOcc_range) | set(Op2._CrtOcc_range) ),reverse=True)
	DstOcc_range = sorted(list( set(Op1._DstOcc_range) | set(Op2._DstOcc_range) ),reverse=True)
	CrtVrt_range = sorted(list( set(Op1._CrtVrt_range) | set(Op2._CrtVrt_range) ),reverse=True)
	DstVrt_range = sorted(list( set(Op1._DstVrt_range) | set(Op2._DstVrt_range) ),reverse=True)
	target = nested_operator(Op1._occ_strings, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
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









