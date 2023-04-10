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
class _iter_done(object):
	def __init__(self):  pass
	def __iter__(self):  return self
	def __next__(self):  raise StopIteration

def _indices(parent):
	all_indices = []
	all_blocks  = []
	#
	if parent._Ex is not None:
		all_indices += [ iter(parent._Tx_indices.Excitations) ]
		all_blocks  += [ True ]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	if parent._Fx is not None:
		all_indices += [ iter(parent._Tx_indices.Flat_transitions) ]
		all_blocks  += [ True ]
	else:
		all_indices += [ _iter_done() ]
		all_blocks  += [ False ]
	if parent._Dx is not None:
		all_indices += [ iter(parent._Tx_indices.Deexcitations) ]
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
		self._doing_Fx = False
		self._doing_Dx = False
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
		if   self._doing_Ex:  self._doing_Ex, self._doing_Fx = False,True
		elif self._doing_Fx:  self._doing_Fx, self._doing_Dx = False,True
		# if it were in Dx, then this condition will never be executed, and this leaves the whole thing in a self-consistent, but exhausted state



class no_mask(object):
	def __call__(self, subop):  return subop

class not_masked(object):
	def __init__(self):
		self.Ex_mask = no_mask()
		self.Fx_mask = no_mask()
		self.Dx_mask = no_mask()

class mask_and_block(object):
	def __init__(self, mask, block):
		self.mask  = mask
		self.block = block



class iter_single(object):
	def __init__(self, parent):
		indices,blocks = _indices(parent)
		self._index = _index_iter(indices)
		self._TxClasses = parent._TxClasses
		mask = parent._mask
		if mask is None:  mask = not_masked()
		self._Ex = mask_and_block(mask.Ex_mask, parent._get_Ex_block()) if blocks[0] else None
		self._Fx = mask_and_block(mask.Fx_mask, parent._get_Fx_block()) if blocks[1] else None
		self._Dx = mask_and_block(mask.Dx_mask, parent._get_Dx_block()) if blocks[2] else None
	def __iter__(self):
		return self
	def __next__(self):
		index = next(self._index)	# may raise StopIteration exception.  No need to handle it, just let it through
		if   self._index._doing_Ex:  return self._Ex.mask(self._Ex.block[index]), self._TxClasses.excitation(index)
		elif self._index._doing_Fx:  return self._Fx.mask(self._Fx.block[index]), self._TxClasses.flat_transition(index)
		elif self._index._doing_Dx:  return self._Dx.mask(self._Dx.block[index]), self._TxClasses.deexcitation(index)

class iter_pair(object):
	def __init__(self, parent1, parent2, mask):
		indices,blocks = mask(parent1, parent2)
		self._index = _index_iter(indices)
		mask1 = parent1._mask
		if mask1 is None:  mask1 = not_masked()
		mask2 = parent2._mask
		if mask2 is None:  mask2 = not_masked()
		self._Ex1 = mask_and_block(mask1.Ex_mask, parent1._get_Ex_block()) if blocks[0] else None
		self._Ex2 = mask_and_block(mask2.Ex_mask, parent2._get_Ex_block()) if blocks[0] else None
		self._Fx1 = mask_and_block(mask1.Fx_mask, parent1._get_Fx_block()) if blocks[1] else None
		self._Fx2 = mask_and_block(mask2.Fx_mask, parent2._get_Fx_block()) if blocks[1] else None
		self._Dx1 = mask_and_block(mask1.Dx_mask, parent1._get_Dx_block()) if blocks[2] else None
		self._Dx2 = mask_and_block(mask2.Dx_mask, parent2._get_Dx_block()) if blocks[2] else None
	def __iter__(self):
		return self
	def __next__(self):
		index = next(self._index)	# may raise StopIteration exception.  No need to handle it, just let it through
		if   self._index._doing_Ex:  return self._Ex1.mask(self._Ex1.block[index]), self._Ex2.mask(self._Ex2.block[index])
		if   self._index._doing_Fx:  return self._Fx1.mask(self._Fx1.block[index]), self._Fx2.mask(self._Fx2.block[index])
		if   self._index._doing_Dx:  return self._Dx1.mask(self._Dx1.block[index]), self._Dx2.mask(self._Dx2.block[index])








def _take_intersection(indices1, indices2, fluctuation_classes):
	indices = list( set(indices1) & set(indices2) )	# sets are not inherently ordered!
	if len(indices)==0:  return _iter_done()
	else:                return iter(fluctuation_classes.order_indices(indices))

def _intersection(parent1, parent2):
	if parent1._TxClasses!=parent2._TxClasses:  raise Exception("illegal operation between different operator types")
	fluctuation_classes = parent1._TxClasses
	all_indices = []
	all_blocks  = []
	#
	if (parent1._Ex is None) or (parent2._Ex is None):  iterator = _iter_done()
	else:                                               iterator = _take_intersection(parent1._Tx_indices.Excitations,      parent2._Tx_indices.Excitations,      fluctuation_classes)	# could return _iter_done instance
	all_indices += [ iterator ]
	if isinstance(iterator,_iter_done):  all_blocks += [ False ]
	else:                                all_blocks += [ True ]
	#
	if (parent1._Fx is None) or (parent2._Fx is None):  iterator = _iter_done()
	else:                                               iterator = _take_intersection(parent1._Tx_indices.Flat_transitions, parent2._Tx_indices.Flat_transitions, fluctuation_classes)	# could return _iter_done instance
	all_indices += [ iterator ]
	if isinstance(iterator,_iter_done):  all_blocks += [ False ]
	else:                                all_blocks += [ True ]
	#
	if (parent1._Dx is None) or (parent2._Dx is None):  iterator = _iter_done()
	else:                                               iterator = _take_intersection(parent1._Tx_indices.Deexcitations,    parent2._Tx_indices.Deexcitations,    fluctuation_classes)	# could return _iter_done instance
	all_indices += [ iterator ]
	if isinstance(iterator,_iter_done):  all_blocks += [ False ]
	else:                                all_blocks += [ True ]
	#		
	return all_indices, all_blocks

def _latter_indices(parent1, parent2):
	all_indices, all_blocks = _indices(parent2)
	if parent1._mask is not None:		# Could lead to crazy results, but if the user masks an increment target, this is what it must mean ... (don't increment, even if present in source)
		if not parent1._mask.Ex_on:  all_indices[0], all_blocks[0] = _iter_done(), False
		if not parent1._mask.Fx_on:  all_indices[1], all_blocks[1] = _iter_done(), False
		if not parent1._mask.Dx_on:  all_indices[2], all_blocks[2] = _iter_done(), False
	return all_indices, all_blocks



#primarily for use with dot product
def common(parent1, parent2, indent=""):  return iter_pair(parent1, parent2, _intersection)

# primarily for use with increment
def latter(parent1, parent2):  return iter_pair(parent1, parent2, _latter_indices)
