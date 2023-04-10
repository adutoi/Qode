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
from base import indexed_base



class Coeff(indexed_base):
	def __init__(self, name, bottom=[], top=[]):		# mutable default args ok because read-only
		indexed_base.__init__(self)
		self.name = name
		self.bottom = list(bottom)	# just a list of Index objects (shallow copy of input)
		self.top    = list(top)		# just a list of Index objects (shallow copy of input)
	def _replace_indices(self, index_replacements):		# private helper function.  in place modification.
		if index_replacements:
			old,new = zip(*index_replacements)
			for indices in (self.bottom, self.top):
				for n in range(len(indices)):
					if indices[n] in old:
						indices[n] = new[old.index(indices[n])]		# this uses list.index ... unfortunate collision of naming conventions
		return self
	def _copy(self, index_replacements):
		if (self is zero) or (self is minus) or (self is one):  return self
		else:  return Coeff(self.name, self.bottom, self.top)._replace_indices(index_replacements)
	def __str__(self):
		string = self.name
		if len(self.bottom)>0:
			string += "_{"
			for index in self.bottom:
				string += index.letter		# Sort of assumes this is being called from inside print(<instance of Sum>)
			string += "}"
		if len(self.top)>0:
			string += "^{"
			for index in self.top:
				string += index.letter		# Sort of assumes this is being called from inside print(<instance of Sum>)
			string += "}"
		return string+" "
	def code(self):
		name = self.name
		indices = []
		if len(self.bottom)>0:
			for i in range(len(self.bottom)):
				indices += [ self.bottom[i].letter, self.top[i].letter ]
		return name, indices

		
zero  = Coeff("(0)")
one   = Coeff("(1)")
minus = Coeff("(-1)")



class delta(Coeff):
	def __init__(self, p, q):
		Coeff.__init__(self, "\\delta", (p,q))
		if p.index_set is not q.index_set:  raise Exception("illegal Kronecker delta")		
	def _copy(self, index_replacements):
		return delta(*self.bottom)._replace_indices(index_replacements)
	def code(self):
		print("There should be NO COMMUTATOR WHEN DONE.")
		raise AssertionError


def condense_deltas(delta_list):		# does not need the delta class, but used in conjunction with it
	delta_list_cp = list(delta_list)	# shallow copy
	union = None
	break_out = False
	for p_idx in range(len(delta_list)):
		for q_idx in range(p_idx+1,len(delta_list)):
			p = delta_list[p_idx]
			q = delta_list[q_idx]
			intersection = list(set(p) & set(q))
			union        = list(set(p) | set(q))
			if intersection:
				del delta_list_cp[q_idx]	# must delete higher index first ...
				del delta_list_cp[p_idx]	# ... so that this one still means what we think it means
				break_out = True
				break
		if break_out:  break
	if break_out:  return condense_deltas(delta_list_cp + [union])
	else:          return delta_list_cp
