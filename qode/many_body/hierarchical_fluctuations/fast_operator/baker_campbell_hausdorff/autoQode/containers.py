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
from base import indexed_base, copy
from scalars import Coeff, zero, one, minus, delta



def commute(A,B):  return A.commute(B)		# by calling a member function, this file does not have to import operators (which needs to import this file to return result ... avoid cyclic import)

def reorder(op_product):
	if op_product.factors[0].is_field_op:
		return op_product.factors[0].reorder_primitives(op_product.factors[1:]), None		# hack ... by calling a member function ...
	if   len(op_product.factors)==1:
		return op_product, False
	elif len(op_product.factors)==2:
		op1,op2 = op_product.factors
		return op1.reorder(op2)		# by calling a member function, this file does not have to import operators (which needs to import this file to return result ... avoid cyclic import)
	else:
		n = 0
		N = len(op_product.factors)
		while n<N-1:
			new,reordered = reorder(mult(op_product.factors[n], op_product.factors[n+1]))
			if reordered:  break
			n += 1
		if reordered:
			factors = []
			if n>0:   factors += op_product.factors[:n]
			factors += [new]
			if n<N-2: factors += op_product.factors[n+2:]
			reordered_expansion = expand(mult(*factors))
			if isinstance(reordered_expansion,container):  return reorder_ops(reordered_expansion), True
			else:                                          return reordered_expansion, True		# might be zero or single operator
		else:
			return op_product, False


# advantage of external function is unbound name to give to _pass_through methods

def do_commutators(container_object,index_replacements=[]):
	return container_object._do_commutators(index_replacements)

def reorder_ops(container_object,index_replacements=[]):
	return container_object._reorder_ops(index_replacements)

def reorder_pdts(container_object,index_replacements=[]):
	return container_object._reorder_pdts(index_replacements)

def nest_summations(container_object,index_replacements=[]):
	return container_object._nest_summations(index_replacements)

def resolve_deltas(container_object,index_replacements=[]):
	return container_object._resolve_deltas(index_replacements)

def expand(container_object,index_replacements=[]):
	return container_object._expand(index_replacements)

def sum_case(container_object,index_replacements=[]):
	return container_object._sum_case(index_replacements)

def trim_high_order(container_object,index_replacements=[]):
	return container_object._trim_high_order(index_replacements)



class container(indexed_base):
	def __init__(self):  indexed_base.__init__(self)
	@staticmethod
	def _pass_through_helper(function, items, index_replacements):
		new_items = []
		for item in items:
			if isinstance(item,container):  new_items += [function(item, index_replacements)]		# this should always return an appropriately modified copy
			else:                           new_items += [copy(item, index_replacements)]
		return new_items
	@staticmethod
	def _string_helper(items, begin="\\left( ", sep=None, end=" \\right)"):
		string = ""
		string += begin + str(items[0])
		for item in items[1:]:  string += sep + str(item)
		string += end
		return string
	def _do_commutators(self, index_replacements):   raise NotImplementedError()
	def _reorder_ops(self, index_replacements):      raise NotImplementedError()
	def _reorder_pdts(self, index_replacements):     raise NotImplementedError()
	def _nest_summations(self, index_replacements):  raise NotImplementedError()
	def _resolve_deltas(self, index_replacements):   raise NotImplementedError()
	def _expand(self, index_replacements):           raise NotImplementedError()
	def _sum_case(self, index_replacements):         raise NotImplementedError()
	def multiline(self,_):  return self	# overridden only by add



preamble = """\
\\documentclass[12pt,letterpaper]{article}
\\setlength{\\pdfpagewidth}{30in}
\\setlength{\\pdfpageheight}{100in}
\\setlength{\\textwidth}{13in}
\\setlength{\\textheight}{98in}
\\usepackage{amssymb}
\\usepackage{amsmath}
\\usepackage[utf8]{inputenc}
\\begin{document}
\\begin{eqnarray}
"""

class add(container):
	def __init__(self, *terms):
		container.__init__(self)
		if not terms:  raise Exception("Empty addition is ambiguous")
		self.terms = list(terms)
		self.mline    = False
		self.selfcont = False
	def _pass_through(self, function, index_replacements):
		terms = container._pass_through_helper(function, self.terms, index_replacements)
		return add(*terms)
	def _do_commutators(self, index_replacements):   return self._pass_through(do_commutators, index_replacements)
	def _reorder_ops(self, index_replacements):      return self._pass_through(reorder_ops, index_replacements)
	def _reorder_pdts(self, index_replacements):     return self._pass_through(reorder_pdts, index_replacements)
	def _nest_summations(self, index_replacements):  return self._pass_through(nest_summations, index_replacements)
	def _resolve_deltas(self, index_replacements):   return self._pass_through(resolve_deltas, index_replacements)
	def _sum_case(self, index_replacements):         return self._pass_through(sum_case, index_replacements)
	def _trim_high_order(self, index_replacements):  return self._pass_through(trim_high_order, index_replacements)
	def _copy(self, index_replacements):             return self._pass_through(copy, index_replacements)			# calls copy on container and primitive terms
	def _expand(self, index_replacements):
		terms = container._pass_through_helper(expand, self.terms, index_replacements)
		if len(terms)==1:
			return terms[0]
		else:
			result = []
			for term in terms:
				if isinstance(term, add):
					for subterm in term.terms:  result += [subterm]
				elif term is not zero:
					result += [term]
			if result:  return add(*result)
			else:       return zero			# can happen that zeros leave result empty
	def __str__(self):
		string = ""
		if self.selfcont:
			begin = preamble
			sep = "\n\\end{eqnarray}\n\\begin{eqnarray}\n + "
			end = "\n\\end{eqnarray}\n\\end{document}\n"
		elif self.mline:  begin, sep, end = "", "  +  \\\\ \n", ""
		else:             begin, sep, end = "\\left( ", "  +  ", " \\right)"
		return container._string_helper(self.terms, begin, sep, end)
	def multiline(self, self_contained=False):
		new = add(*self.terms)
		new.mline = True
		if self_contained:  new.selfcont = True
		return new
	def code(self):
		add_returns = []
		if len(self.terms) > 0:
			for each_term in self.terms:
				print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
				add_returns += [each_term.code()]
		return add_returns


class mult(container):
	def __init__(self, *factors):
		container.__init__(self)
		if not factors:  raise Exception("Empty multiplication is ambiguous")
		self.factors = list(factors)
	def _pass_through(self, function, index_replacements):
		factors = container._pass_through_helper(function, self.factors, index_replacements)
		#if len(factors)==1:  return factors[0]
		#else:                return mult(*factors)
		return mult(*factors)
	def _do_commutators(self, index_replacements):   return self._pass_through(do_commutators, index_replacements)
	def _nest_summations(self, index_replacements):  return self._pass_through(nest_summations, index_replacements)
	def _resolve_deltas(self, index_replacements):   return self._pass_through(resolve_deltas, index_replacements)
	def _sum_case(self, index_replacements):         return self._pass_through(sum_case, index_replacements)
	def _copy(self, index_replacements):             return self._pass_through(copy, index_replacements)			# calls copy on container and primitive terms
	@staticmethod
	def _pdt_organizer(factors):
		phase = 1
		kroneckers   = []
		tensors      = []
		fluctuations = []
		for factor in factors:
			if factor is minus:		# intentionally ignores "one"
				phase *= -1
			elif isinstance(factor,delta):
				kroneckers   += [factor]
			elif isinstance(factor,Coeff):
				tensors      += [factor]
			else:					# must be an operator, but can't import that here (circular)
				fluctuations += [factor]
		if phase==-1:  tensors = [minus] + tensors
		products = []
		if kroneckers:    products += [mult(*kroneckers)]
		if tensors:       products += [mult(*tensors)]
		if fluctuations:  products += [mult(*fluctuations)]
		if not products:  products = [one]
		return kroneckers, tensors, fluctuations, products
	def _reorder_ops(self, index_replacements):
		kroneckers, tensors, fluctuations, products = mult._pdt_organizer(container._pass_through_helper(reorder_ops, self.factors, index_replacements))
		products = []
		if kroneckers:    products += [mult(*kroneckers)]
		if tensors:       products += [mult(*tensors)]
		if fluctuations:
			operators,_ = reorder(mult(*fluctuations))
			products += [operators]
		if not products:  products = [one]
		return mult(*products)
	def _reorder_pdts(self, index_replacements):
		kroneckers, tensors, fluctuations, products = mult._pdt_organizer(container._pass_through_helper(reorder_pdts, self.factors, index_replacements))
		return mult(*products)
	def _trim_high_order(self, index_replacements):
		kroneckers, tensors, fluctuations, products = mult._pdt_organizer(container._pass_through_helper(trim_high_order, self.factors, index_replacements))
		if len(fluctuations)>2:
			if fluctuations[0].is_Ex and fluctuations[1].is_Ex and fluctuations[2].is_not_Dx:  return zero		# very CCSD specific (must use member to avoid importing operators.py)
		return mult(*products)
	def _expand(self, index_replacements):
		factors = container._pass_through_helper(expand, self.factors, index_replacements)
		if   len(factors)==1:
			return factors[0]
		elif len(factors)==2:
			A, B = factors
			if (A is zero) or (B is zero):
				return zero
			elif (A is one):
				return B
			elif (B is one):
				return A
			elif isinstance(A,add):
				terms = []
				for a in A.terms:
					terms += [ expand(mult(a,B)) ]
				return add(*terms)
			elif isinstance(B,add):
				terms = []
				for b in B.terms:
					terms += [ expand(mult(A,b)) ]
				return add(*terms)
			else:
				result = []
				for factor in (A,B):		# only 2 things in factors, but loop still convenient
					if isinstance(factor, mult):
						for subfactor in factor.factors:  result += [subfactor]
					else:  result += [factor]
				return mult(*result)
		elif len(factors)>2:
			result = factors[0]
			for factor in factors[1:]:  result = expand(mult(result, factor))
			return result
	def __str__(self):
		return container._string_helper(self.factors, sep=" ")
	def code(self):
		mult_returns = []
		# print("IN MULT OBJ self.factors =", self.factors)
		if len(self.factors) > 0:
			for each_factor in self.factors:
				mult_returns += [each_factor.code()]
		return mult_returns




class comm(container):
	def __init__(self, A, B):
		container.__init__(self)
		self.components = A,B
	def _pass_through(self, function, index_replacements):
		A,B = container._pass_through_helper(function, self.components, index_replacements)
		return comm(A,B)
	def _reorder_ops(self, index_replacements):      return self._pass_through(reorder_ops, index_replacements)
	def _reorder_pdts(self, index_replacements):     return self._pass_through(reorder_pdts, index_replacements)
	def _nest_summations(self, index_replacements):  return self._pass_through(nest_summations, index_replacements)
	def _resolve_deltas(self, index_replacements):   return self._pass_through(resolve_deltas, index_replacements)
	def _sum_case(self, index_replacements):         return self._pass_through(sum_case, index_replacements)
	def _trim_high_order(self, index_replacements):  return self._pass_through(trim_high_order, index_replacements)
	def _copy(self, index_replacements):             return self._pass_through(copy, index_replacements)			# calls copy on container and primitive terms
	def _do_commutators(self, index_replacements):
		A,B = container._pass_through_helper(do_commutators, self.components, index_replacements)
		return commute(A,B)										# a little dirty.  this will fail unless completely expanded
	def _expand(self, index_replacements):
		A,B = container._pass_through_helper(expand, self.components, index_replacements)
		if isinstance(A,Coeff) or isinstance(B,Coeff):
			return zero
		elif isinstance(A,add):
			terms = []
			for a in A.terms:
				terms += [ expand(comm(a,B)) ]
			return add(*terms)
		elif isinstance(B,add):
			terms = []
			for b in B.terms:
				terms += [ expand(comm(A,b)) ]
			return add(*terms)
		elif isinstance(A,mult):
			terms = []
			a     = A.factors[0]
			Atail = A.factors[1:]
			terms += [mult( comm(a,B), mult(*Atail) )]
			for n in range(1, len(A.factors)-1):
				Ahead = A.factors[:n]
				a     = A.factors[n]
				Atail = A.factors[n+1:]
				terms += [mult( mult(*Ahead), comm(a,B), mult(*Atail) )]
			Ahead = A.factors[:-1]
			a     = A.factors[-1]
			terms += [mult( mult(*Ahead), comm(a,B) )]
			return add(*terms)
		elif isinstance(B,mult):
			terms = []
			b     = B.factors[0]
			Btail = B.factors[1:]
			terms += [mult( comm(A,b), mult(*Btail) )]
			for n in range(1, len(B.factors)-1):
				Bhead = B.factors[:n]
				b     = B.factors[n]
				Btail = B.factors[n+1:]
				terms += [mult( mult(*Bhead), comm(A,b), mult(*Btail) )]
			Bhead = B.factors[:-1]
			b     = B.factors[-1]
			terms += [mult( mult(*Bhead), comm(A,b) )]
			return add(*terms)
		else:
			return comm(A,B)
	def __str__(self):
		return container._string_helper(self.components, begin="\\left[ ", sep=" , ", end=" \\right]")
	def code(self):
		print("There should be NO COMMUTATOR WHEN DONE.")
		raise AssertionError
