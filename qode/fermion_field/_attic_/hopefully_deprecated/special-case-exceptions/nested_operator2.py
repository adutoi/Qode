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
from copy import deepcopy #only for debugging
from .state                     import state, dot, create, annihilate, op_string
from .reorder2                  import excitation, deexcitation, occ_rearrange, vrt_rearrange, reorder_t_t, OpPdtIsZero, OpPdtIsOrdered
from .nested_operator_iterators import iter_single, common, latter



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
	def __init__(self,occ_strings,CrtOcc_range,DstOcc_range,CrtVrt_range,DstVrt_range,C=0.):
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
		self._C = C
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
	def __iter__(self):
		return iter_single(self)
	def __call__(self,the_state):
		""" Acts the present layer on the state given by recursively calling lower layers and then operating with the enclosing operator """
		new_state = state(self._occ_strings)
		new_state.increment(the_state,self._C)
		for suboperator,transition in self:  new_state.increment(transition(suboperator(the_state)))
		return new_state
	def scale(self,c):
		self._C *= c
		for suboperator,_ in self:  suboperator.scale(c)
	def dot(self,other):
		product = 0
		product += self._C * other._C
		for subop1,subop2 in common(self,other):  product += subop1.dot(subop2)			# loops only over indices in common
		return product
	def increment(self,other,c=1):
		""" Add another nested operator to this one """
		self._C += c * other._C
		for subop1,subop2 in latter(self,other):  subop1.increment(subop2,c)			# loops over indices present in latter, getting (implicitly making) subblocks in former as necessary
	def add_term(self,term):
		"""\
		term comes in as a list of operators, presumed to obey the ordering restrictions, which terminates with a scalar the multiplies that string.
		The task is essentially to find the right place in the recursion scheme to place the scalar at the end, which it adds to whatever is there,
		creating blocks along the way if it is needed (works by recursion).
		"""
		try:
			self._C += term[0]
		except:
			self._subtarget(term[0]).add_term(term[1:])
	def right_multiply(self, transition):
		""" Return a new nested operator that obeys proper index ordering and is equivalent to this operator multipled from the right by a given 1e- number-conserving transition """
		CrtOcc_range = sorted(list( set(self._CrtOcc_range) | set(transition.CrtOcc) ),reverse=True)
		DstOcc_range = sorted(list( set(self._DstOcc_range) | set(transition.DstOcc) ),reverse=True)
		CrtVrt_range = sorted(list( set(self._CrtVrt_range) | set(transition.CrtVrt) ),reverse=True)
		DstVrt_range = sorted(list( set(self._DstVrt_range) | set(transition.DstVrt) ),reverse=True)
		target = nested_operator(self._occ_strings, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
		target.add_term([transition, self._C])
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
		target.add_term([transition, self._C])
		for source,inner_trans in self:
			try:
				products = reorder_t_t(transition, inner_trans)
			except OpPdtIsZero:
				pass
			except OpPdtIsOrdered:
				target._subtarget(transition)._subtarget(inner_trans).increment(source)
			else:
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





def multiply(Op1, Op2, use_left=True):
	CrtOcc_range = sorted(list( set(Op1._CrtOcc_range) | set(Op2._CrtOcc_range) ),reverse=True)
	DstOcc_range = sorted(list( set(Op1._DstOcc_range) | set(Op2._DstOcc_range) ),reverse=True)
	CrtVrt_range = sorted(list( set(Op1._CrtVrt_range) | set(Op2._CrtVrt_range) ),reverse=True)
	DstVrt_range = sorted(list( set(Op1._DstVrt_range) | set(Op2._DstVrt_range) ),reverse=True)
	target = nested_operator(Op1._occ_strings, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
	#
	target.increment(Op2, Op1._C)
	if use_left:
		for suboperator,transition in zip_ops(suboperators(Op1), transitions(Op1)):
			target.increment(multiply(suboperator,Op2).left_multiply(transition))
	else:
		for suboperator,transition in zip_ops(suboperators(Op2), transitions(Op2)):
			target.increment(multiply(Op1,suboperator).right_multiply(transition))
	return target
