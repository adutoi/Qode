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
import sys
from .reorder import excitation, deexcitation, occ_rearrange, vrt_rearrange



class _iter_done(object):
	def __init__(self):  pass
	def __iter__(self):  return self
	def __next__(self):  raise StopIteration

def _Ex_list(parent):
	composite_indices = []				# This will be small, so readability trumps inefficiency of storing (and may be more time efficient)
	for a in parent._CrtVrt_range:
		for i in parent._DstOcc_range:
			composite_indices += [(a,i)]
	return composite_indices

def _Dx_list(parent):
	composite_indices = []				# This will be small, so readability trumps inefficiency of storing (and may be more time efficient)
	for a in parent._DstVrt_range:
		for i in parent._CrtOcc_range:
			composite_indices += [(a,i)]
	return composite_indices

def _Fo_list(parent):
	composite_indices = []				# This will be small, so readability trumps inefficiency of storing (and may be more time efficient)
	for i in parent._DstOcc_range:
		for j in parent._CrtOcc_range:
			composite_indices += [(i,j)]
	return composite_indices

def _Fv_list(parent):
	composite_indices = []				# This will be small, so readability trumps inefficiency of storing (and may be more time efficient)
	for a in parent._CrtVrt_range:
		for b in parent._DstVrt_range:
			composite_indices += [(a,b)]
	return composite_indices

def _indices(parent):
	all_indices = []
	all_blocks  = []
	#
	if parent._Ex is not None:
		all_indices += [ iter(_Ex_list(parent)) ]
		all_blocks  += [ True ]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	if parent._Dx is not None:
		all_indices += [ iter(_Dx_list(parent)) ]
		all_blocks  += [ True ]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	if parent._Fo is not None:
		all_indices += [ iter(_Fo_list(parent)) ]
		all_blocks  += [ True ]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	if parent._Fv is not None:
		all_indices += [ iter(_Fv_list(parent)) ]
		all_blocks  += [ True ]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	return all_indices, all_blocks



class _index_iter(object):
	def __init__(self, iterator_list):
		self._block = iter(iterator_list)
		self._index = next(self._block)	# Starts with excitations
		self._doing_Ex = True		# might already be "done", consistent with fact that doing_Fv will be true when the whole thing is done
		self._doing_Dx = False
		self._doing_Fo = False
		self._doing_Fv = False
	def __iter__(self):
		return self
	def __next__(self):
		try:
			return next(self._index)
		except StopIteration:
			self._next_block()
			return next(self)
	def _next_block(self):
		self._index = next(self._block)
		if   self._doing_Ex:  self._doing_Ex, self._doing_Dx = False,True
		elif self._doing_Dx:  self._doing_Dx, self._doing_Fo = False,True
		elif self._doing_Fo:  self._doing_Fo, self._doing_Fv = False,True
		# if it were in Fv, then this condition will never be executed, and this leaves the whole thing in a self-consistent, but exhausted state



class no_mask(object):
	def __call__(self, subop):  return subop

class not_masked(object):
	def __init__(self):
		self.Ex_mask = no_mask()
		self.Dx_mask = no_mask()
		self.Fo_mask = no_mask()
		self.Fv_mask = no_mask()

class mask_and_block(object):
	def __init__(self, mask, block):
		self.mask  = mask
		self.block = block



class iter_single(object):
	def __init__(self, parent):
		indices,blocks = _indices(parent)
		self._index = _index_iter(indices)
		mask = parent._mask
		if mask is None:  mask = not_masked()
		self._Ex = mask_and_block(mask.Ex_mask, parent._get_Ex_block()) if blocks[0] else None
		self._Dx = mask_and_block(mask.Dx_mask, parent._get_Dx_block()) if blocks[1] else None
		self._Fo = mask_and_block(mask.Fo_mask, parent._get_Fo_block()) if blocks[2] else None
		self._Fv = mask_and_block(mask.Fv_mask, parent._get_Fv_block()) if blocks[3] else None
	def __iter__(self):
		return self
	def __next__(self):
		p,q = next(self._index)	# may raise StopIteration exception.  No need to handle it, just let it through
		if   self._index._doing_Ex:  return self._Ex.mask(self._Ex.block[p][q]),    excitation(p,q)
		elif self._index._doing_Dx:  return self._Dx.mask(self._Dx.block[p][q]),  deexcitation(p,q)
		elif self._index._doing_Fo:  return self._Fo.mask(self._Fo.block[p][q]), occ_rearrange(p,q)
		elif self._index._doing_Fv:  return self._Fv.mask(self._Fv.block[p][q]), vrt_rearrange(p,q)

class iter_pair(object):
	def __init__(self, parent1, parent2, mask, indent=""):
		indices,blocks = mask(parent1, parent2)
		self._index = _index_iter(indices)
		mask1 = parent1._mask
		if mask1 is None:  mask1 = not_masked()
		mask2 = parent2._mask
		if mask2 is None:  mask2 = not_masked()
		self._Ex1 = mask_and_block(mask1.Ex_mask, parent1._get_Ex_block()) if blocks[0] else None
		self._Ex2 = mask_and_block(mask2.Ex_mask, parent2._get_Ex_block()) if blocks[0] else None
		self._Dx1 = mask_and_block(mask1.Dx_mask, parent1._get_Dx_block()) if blocks[1] else None
		self._Dx2 = mask_and_block(mask2.Dx_mask, parent2._get_Dx_block()) if blocks[1] else None
		self._Fo1 = mask_and_block(mask1.Fo_mask, parent1._get_Fo_block()) if blocks[2] else None
		self._Fo2 = mask_and_block(mask2.Fo_mask, parent2._get_Fo_block()) if blocks[2] else None
		self._Fv1 = mask_and_block(mask1.Fv_mask, parent1._get_Fv_block()) if blocks[3] else None
		self._Fv2 = mask_and_block(mask2.Fv_mask, parent2._get_Fv_block()) if blocks[3] else None
		self.indent=indent
	def __iter__(self):
		return self
	def __next__(self):
		p,q = next(self._index)	# may raise StopIteration exception.  No need to handle it, just let it through
		if   self._index._doing_Ex:
			#print(self.indent+"Ex",p,q)
			return self._Ex1.mask(self._Ex1.block[p][q]), self._Ex2.mask(self._Ex2.block[p][q])
		if   self._index._doing_Dx:
			#print(self.indent+"Dx",p,q)
			return self._Dx1.mask(self._Dx1.block[p][q]), self._Dx2.mask(self._Dx2.block[p][q])
		if   self._index._doing_Fo:
			#print(self.indent+"Fo",p,q)
			return self._Fo1.mask(self._Fo1.block[p][q]), self._Fo2.mask(self._Fo2.block[p][q])
		if   self._index._doing_Fv:
			#print(self.indent+"Fv",p,q)
			return self._Fv1.mask(self._Fv1.block[p][q]), self._Fv2.mask(self._Fv2.block[p][q])








def _take_intersection(parent1, parent2, get_index_list):
	indices1 = get_index_list(parent1)
	indices2 = get_index_list(parent2)
	indices = list( set(indices1) & set(indices2) )		# sets are not inherently ordered!
	indices.sort()						# tuples have inherent python ordering
	if len(indices)==0:  return _iter_done()
	else:                return reversed(indices)		# reversed returns an iterator

def _intersection(parent1, parent2):
	all_indices = []
	all_blocks  = []
	#
	if parent1._Ex is not None and parent2._Ex is not None:
		all_indices += [ _take_intersection(parent1, parent2, _Ex_list) ]
		all_blocks  += [ True ]		# Could result in creation of unnecessary zero block _iter_done() returned above [have _take_intersection raise exception?]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	if parent1._Dx is not None and parent2._Dx is not None:
		all_indices += [ _take_intersection(parent1, parent2, _Dx_list) ]
		all_blocks  += [ True ]		# Could result in creation of unnecessary zero block _iter_done() returned above [have _take_intersection raise exception?]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	if parent1._Fo is not None and parent2._Fo is not None:
		all_indices += [ _take_intersection(parent1, parent2, _Fo_list) ]
		all_blocks  += [ True ]		# Could result in creation of unnecessary zero block _iter_done() returned above [have _take_intersection raise exception?]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	if parent1._Fv is not None and parent2._Fv is not None:
		all_indices += [ _take_intersection(parent1, parent2, _Fv_list) ]
		all_blocks  += [ True ]		# Could result in creation of unnecessary zero block _iter_done() returned above [have _take_intersection raise exception?]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	return all_indices, all_blocks

def _latter_indices(parent1, parent2):
	all_indices, all_blocks = _indices(parent2)
	if parent1._mask is not None:		# Could lead to crazy results, but if the user masks an increment target, this is what it must mean ... (don't increment, even if present in source)
		if not parent1._mask.Ex_on:  all_indices[0], all_blocks[0] = _iter_done(), False
		if not parent1._mask.Dx_on:  all_indices[1], all_blocks[1] = _iter_done(), False
		if not parent1._mask.Fo_on:  all_indices[2], all_blocks[2] = _iter_done(), False
		if not parent1._mask.Fv_on:  all_indices[3], all_blocks[3] = _iter_done(), False
	return all_indices, all_blocks



#primarily for use with dot product
def common(parent1, parent2, indent=""):  return iter_pair(parent1, parent2, _intersection, indent)

# primarily for use with increment
def latter(parent1, parent2):  return iter_pair(parent1, parent2, _latter_indices)







# Deprecated and belongs with older interface ... turned out to be really dangerous to zip two operators because even if they were compatible, they 
# might have much different patterns of "holes" in their index structure
#class zip_ops(object):
#	"""\
#	Like zip, this terminates as soon as the end of any iterator is reached, so it sort of assumes the user has things lined up (operators of same size and masking).
#	I modeled the interface on the builtin zip, but I was a little confused by the combination of return while and yield, so I just hacked it as best I could.
#	"""
#	def __init__(self,*iterable_ops):
#		self._iterators = [iter(it) for it in iterable_ops]
#	def __iter__(self):
#		return self
#	def __next__(self):
#		return tuple([next(it) for it in self._iterators])
