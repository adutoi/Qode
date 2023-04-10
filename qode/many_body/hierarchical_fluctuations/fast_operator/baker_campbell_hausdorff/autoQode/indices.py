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
class IndexRelationInconsistency(Exception):
	def __init__(self):  pass
	def __str__(self):  return "You should not be seeing this because it should have been handled."



class Index(object):
	""" The main job of this class is to keep track of relationships to other indices, these should always be reciprocated within a group, but there can relationships to any other ("external") index too """
	def __init__(self, index_set, relationships):
		""" should not be called by an end user, access through ._new_index member function of Indices"""
		self.index_set = index_set			# The set to which it belongs (thought of as attached to a specific summation ... can only belong to one set)
		self.relationships = list(relationships)	# list of tuples, each tuple contains a reference to another index, and a character '<', '>', '=' indicating its relationship to that index (in opposite order)
		for relationship in relationships:  self.reciprocate(*relationship)
		self.letter = ''				# volatile, only meainingful over the course of a single "print" (general output) of an index set
	def compare(self, other):
		result = None
		for relate,index in self.relationships:
			if index is other:	# could be in here multiple times, so double-check consistency
				if   result is None:  result = relate
				elif result!=relate:  raise Exception("conflicting index relationships")
		return result
	def add_relationship(self, relate, other):
		self.relationships += [(relate,other)]
		self.reciprocate(relate,other)
	def reciprocate(self, relate, other):
		if other in self.index_set.indices:
			if relate=='=':   reciprocal = '='
			if relate=='>':   reciprocal = '<'
			if relate=='<':   reciprocal = '>'
			if relate=='>=':  reciprocal = '<='
			if relate=='<=':  reciprocal = '>='
			other.relationships += [(reciprocal,self)]
	### Primitive functions  --- They do not make use of add_relationship and reciprocation because they are utilities for replacing existing indices and relationships
	@staticmethod
	def inherit_relationships(parents):
		""" this sense of inheritance assumes (imminent) death of the parents.  we give all the relatives we are still in touch with our phone number and tell them to call us instead.  Also assumes self is newborn: no existing relationships. """
		index_set = parents[0].index_set
		for parent in parents:
			if parent.index_set is not index_set:  raise Exception("illegal Kronecker")
			for relate,other in parent.relationships:
				if other in parents:
					if relate!="=":  raise IndexRelationInconsistency()	# child is supposed to replace parents, so if parents not equivlent, its asking child not to be itself, so summation over no terms is zero
		child = index_set._new_index([])
		for index in index_set.indices:
			for parent in parents:
				index.replace_in_relationships(parent,child)			# all relatives now have our phone number instead of the parents, intraparent relationships (equal) also point now to child
		for parent in parents:
			for relationship in parent.relationships:
				child.relationships += [relationship]				# finally, copy updated relationships to child (parents now prepared to pass on), might have inherited conflicting relationships
		for parent in parents:  index_set._remove(parent)				# parents now buried and forgotten
		redundant = []
		for n,(relate1,other1) in enumerate(child.relationships):
			if other1 is child:  redundant += [n]					# can happen because parents updated before passing on relationships, but must test as equal (enforced above)
			for relate2,other2 in child.relationships:
				if other1 is other2:
					if relate1!=relate2:  raise IndexRelationInconsistency()	# check logical consistency of multiple relationships (but too hard to eliminate this redundancy here for everyone (and not necessary)
		for n in reversed(redundant):  del child.relationships[n]				# get rid of completely spurious self-relationships though (do backwards so unprocessed indices do not change meaning)
		return child
	def replace_in_relationships(self, old, new):
		""" in place replacement!, assume replacing index already inherited reciprocations from indices it is replacing """
		for n,relationship in enumerate(self.relationships):
			relate,other = relationship
			if other is old:
				self.relationships[n] = (relate,new)
	def __str__(self):
		return self.letter				# volatile, only meainingful over the course of a single "print" (general output) of an index set



class comparable(object):
	" a thin wrapper, because if these are implemented directly in the Index object, they interfere with set logic (needed for resolving multiple deltas) because relationships can be fuzzy (ie, None) "
	def __init__(self, index):
		self.index = index			# go ahead an access this for readout, since it never changes
	def __eq__(self,other):
		relate = self.index.compare(other.index)
		if relate=='=':  return True
		else:            return False		# also if no relationship, must be definitive to be true
	def __neq__(self,other):
		relate = self.index.compare(other.index)
		if relate in ('<','>'):  return True
		else:                    return False		# also if no relationship, must be definitive to be true
	def __lt__(self,other):
		relate = self.index.compare(other.index)
		if relate=='<':  return True
		else:            return False		# also if no relationship, must be definitive to be true
	def __gt__(self,other):
		relate = self.index.compare(other.index)
		if relate=='>':  return True
		else:            return False		# also if no relationship, must be definitive to be true
	def __lte__(self,other):
		relate = self.index.compare(other.index)
		if relate in ('<','=','<='):  return True
		else:                         return False		# also if no relationship, must be definitive to be true
	def __gte__(self,other):
		relate = self.index.compare(other.index)
		if relate in ('>','=','>='):  return True
		else:                         return False		# also if no relationship, must be definitive to be true



class Indices(object):
	""" This just holds a list of indices and a list of letters to use (in order) to label them as they appear in summations """
	def __init__(self, how_many=0, bound="N", letters=["n1", "n2", "n3"]):
		self.indices = []
		self.bound = bound
		self.letters = list(letters)
		if isinstance(how_many,str):
			self.letters = []
			parse_string = "".join([(' ' if c==',' else c) for c in how_many])
			tokens = parse_string.split()
			for token in tokens:
				n = 0
				while n<len(token) and token[n] not in "<=>":  n += 1
				self.letters += [token[0:n]]
				relationships = []
				if len(token[n:])>0:
					conditions = token[n:].split('&')
					for condition in conditions:
						n = 0
						while n<len(condition) and condition[n] in "<=>":  n += 1
						relate = condition[0:n]
						other = self.indices[self.letters.index(condition[n:])]
						relationships += [(relate,other)]
				self._new_index(relationships)		# automatically added to list
		elif isinstance(how_many,list):
			for index in how_many:
				index.index_set = self
				self.indices += [index]
		elif isinstance(how_many,int):
			for _ in range(how_many):  self._new_index()
		else:
			raise ValueError()
	def __call__(self):
		return list(self.indices)	# shallow copy
	def _new_index(self, relationships):
		index = Index(self,relationships)
		self.indices += [index]
		return index
	def _remove(self,index):
		self.indices.remove(index)		# does not remove relationships bc that index might now be external
	def _replace_in_relationships(self, index_pairs):
		""" In all relationships, replace the indicated old index in each pair with the new one """
		for index in self.indices:
			for old,new in index_pairs:
				index.replace_in_relationships(old,new)
	def copy(self, index_replacements):						# purposely not called _copy because not to be confused with things derived from indexed_base
		""" Make a copy of each index and make them have the same relationships among each other as in the original set, keeping relationships to external indices exactly as they were"""
		new = Indices(bound=self.bound, letters=self.letters)
		for index in self.indices:  new._new_index(index.relationships)		# relationships of new indices all point to old set
		new._replace_in_relationships(list(zip(self.indices, new.indices)))	# replace references to those old indices by the new ones in the same order to duplicate relationship pattern
		new._replace_in_relationships(index_replacements)			# this copy might also need to replace references to "external" indices
		return new
	def merge_indices(self,parents):
		merged = self._new_index().inherit_relationships(parents)
		for parent in parents:  self.remove(parent)
		return merged
	def assign_letters(self):
		if logic_check(self.indices):  return self.assign_letters_cased()	# strict ... if it does not check out, it just means indices not fully cased (not that there are no logically allowed cases)
		else:	                       return self.assign_letters_uncased()
	def assign_letters_cased(self):
		""" all indiced will be assigned letters. which letter and which restrictions are printed depends on order, but meaning of result does not, so can reorder indices before this """
		letter = iter(self.letters)			# better than zipping below because would rather run out of letters than miss an index
		char_bank = []
		for index in self.indices:
			index.letter = next(letter)
			chars = [str(index), self.bound]	# bound here just for convenience
			priors = []						# priors holds all indices this one has a relationship to that are already given a letter, including outer external indices
			for n,(relate,other) in enumerate(index.relationships):
				if str(other) is not '':
					keep = True
					for rel,idx in index.relationships[:n]:
						if idx is other:
							if rel==relate:  keep = False
							else:  raise Exception("relationship contradiction encountered")
					if keep:  priors += [other]
			priors = order_indices(priors, strict=False)		# since each other has a relationship to this index, they must have relationships between each other
			n = find_location(index, priors, strict=False)		# where the present index belongs in the ordered list of all others already named
			relate = ''
			if n<len(priors):
				relate = index.compare(priors[n])
				chars += [( relate, str(priors[n]) )]
			if n>0 and relate!='=':
				chars += [( index.compare(priors[n-1]), str(priors[n-1]) )]
			char_bank += [chars]
		return char_bank
	def assign_letters_uncased(self):
		""" all indiced will be assigned letters. which letter and which restrictions are printed depends on order, but meaning of result does not, so can reorder indices before this """
		letter = iter(self.letters)		# better than zipping below because would rather run out of letters than miss an index
		char_bank = []
		for index in self.indices:
			index.letter = next(letter)	# the index hereby remembers the letter by which it is called when nested str calls made on the summand
			chars = [str(index), self.bound]
			for n,(relate,other) in enumerate(index.relationships):
				if str(other) is not '':
					keep = True
					for rel,idx in index.relationships[:n]:
						if idx is other:
							if rel==relate:  keep = False
							else:  raise Exception("relationship contradiction encountered")
					if keep:
						chars += [( relate, str(other) )]
			char_bank += [chars]
		return char_bank



def case_indices(index_set):
	cases = [index_set]
	for n in range(len(index_set.indices)):
		cases = case_index(n, cases)
		logical_cases = []
		for case in cases:
			if logic_check(case.indices[:n+1]):  logical_cases += [case]		# used to do this at the end without slicing, but this makes it much faster
		cases = logical_cases
	return cases

def case_index(n, cases_so_far):
	cases = []
	for index_set in cases_so_far:
		index = index_set.indices[n]				# this call focuses on this index in each case
		list_of_relationship_lists = None
		for other in index_set.indices[:n]:			# search over indices already handled (i.e., casing is order dependent!)
			if index.compare(other) is None:					# only need multiple cases if relationship not well defined
				more_relationships = [('<',other), ('=',other), ('>',other)]	# and there are three possibilities
				if list_of_relationship_lists is None:
					list_of_relationship_lists = [[cond] for cond in more_relationships]		# inner list of simultaneous relationships contains only one relationship so far
				else:
					new_list_of_relationship_lists = []
					for cond1 in list_of_relationship_lists:
						for cond2 in more_relationships:
							new_list_of_relationship_lists += [cond1 + [cond2]]	# tensor product of all previous relationship cases with each new relationship
					list_of_relationship_lists = new_list_of_relationship_lists
		if list_of_relationship_lists is None:
			cases += [index_set]
		else:
			for relationships in list_of_relationship_lists:
				cased_index_set = index_set.copy([])
				old = index_set.indices
				new = cased_index_set.indices
				for relate,other in relationships:  cased_index_set.indices[n].add_relationship(relate, new[old.index(other)])
				cases += [cased_index_set]
	return cases

def find_location_ascending(index, indices, strict=True):			# sometimes not strict if trying to arrange inner indices with outers that may have forgotten about them
	""" return the place to insert index in indices such that it will preserve the ascending order of indices with the new index occuring first among any equivalent indices """
	#print("     placing", index)
	n = 0
	relate1 = ''
	while n<len(indices):
		#print("     relative to:", indices[n])
		relate1 = index.compare(indices[n])
		#print("     relate1=",relate1)
		if not strict and relate1 is None:
			relate1 = indices[n].compare(index)
			if relate1 is None:  raise IndexRelationInconsistency()
			if   relate1=="<":  relate1=">"
			elif relate1==">":  relate1="<"
		#print("     relate1=",relate1)
		if relate1!='>':  break
		n += 1
	m = n		# the money line !! (remainder is consistency check)
	#print("m:", m)
	n = m-1
	#print("            n,m,len=",n, m, len(indices))
	while n>=0 and n<len(indices) and m<len(indices):
		inter_relate = indices[m].compare(indices[n])
		if not strict and inter_relate is None:
			inter_relate = indices[n].compare(indices[m])
			if inter_relate is None:  raise IndexRelationInconsistency()
			if   inter_relate=="<":  inter_relate=">"
			elif inter_relate==">":  inter_relate="<"
		if inter_relate!='=':  break
		relate2 = index.compare(indices[n])
		if not strict and relate2 is None:
			relate2 = indices[n].compare(index)
			if relate2 is None:  raise IndexRelationInconsistency()
			if   relate2=="<":  relate2=">"
			elif relate2==">":  relate2="<"
		if relate2!=relate1:  raise IndexRelationInconsistency()
		n -= 1
	n = m+1
	#print("             n,m,len=",n, m, len(indices))
	while n>=0 and n<len(indices) and m<len(indices):
		inter_relate = indices[m].compare(indices[n])
		if not strict and inter_relate is None:
			inter_relate = indices[n].compare(indices[m])
			if inter_relate is None:  raise IndexRelationInconsistency()
			if   inter_relate=="<":  inter_relate=">"
			elif inter_relate==">":  inter_relate="<"
		if inter_relate!='=':  break
		relate2 = index.compare(indices[n])
		if not strict and relate2 is None:
			relate2 = indices[n].compare(index)
			if relate2 is None:  raise IndexRelationInconsistency()
			if   relate2=="<":  relate2=">"
			elif relate2==">":  relate2="<"
		if relate2!=relate1:  raise IndexRelationInconsistency()
		n += 1
	while n<len(indices) and m<len(indices):
		relate3 = index.compare(indices[n])
		if not strict and relate3 is None:
			relate3 = indices[n].compare(index)
			if relate3 is None:  raise IndexRelationInconsistency()
			if   relate3=="<":  relate3=">"
			elif relate3==">":  relate3="<"
		if relate3!='<':  raise IndexRelationInconsistency()
		n += 1
	return m

def find_location(index, indices, descending=False, strict=True):		# descending means indices already in descending order
	""" if descending new index position puts it last among any equivalent indices, but first for ascending """
	if descending:  indices = list(reversed(indices))
	n = find_location_ascending(index, indices, strict)
	if descending:  return len(indices)-n
	else:           return n

def order_indices(indices, descending=False, strict=True):
	ordered_indices = []
	for index in indices:
		n = find_location(index, ordered_indices, strict=strict)	# do in ascending order, then reverse later if needed
		ordered_indices.insert(n,index)
	if descending:  ordered_indices = list(reversed(ordered_indices))
	#print("     ", end="" )
	#for index in ordered_indices:  print(index,end="")
	#print()
	return ordered_indices

def logic_check(indices):
	try:
		order_indices(indices)		# strict, order is irrelevant
	except IndexRelationInconsistency:
		return False
	else:
		return True



def nest_lists(list_of_lists):
	if len(list_of_lists)>1:
		top_list  = list_of_lists[0]
		sub_lists = nest_lists(list_of_lists[1:])
		#print(len(top_list))
		#print(len(sub_lists))
		if top_list and sub_lists:
			new = []
			for item in top_list:
				for sub_list in sub_lists:
					new += [[item] + sub_list]
			return new
		else:
			raise Exception()
	else:
		return [[item] for item in list_of_lists[0]]
