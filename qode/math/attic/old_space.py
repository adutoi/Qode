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
import  math
import cmath
from . import field_traits



class _operator_base(object):
	"""\
	This empty class exists only so that we can test if something is an operator (e.g., from inside the vector class below),
	regardless of whether it is a primitive operator or a composite operator (defined later).
	"""
	pass



class vector(object):
	"""\
	+ Universal base class for user-implemented vector classes.
	+ Deriving from this class maps the actions v|w and v+w, etc, to static virtual functions that the user implements in the child class.
	+ A child class should implement the functions:
	    - vecs_are_compatible(v,w):  Returns a boolean; True if v and w belong to the same vector space (user creativity needed for 
	                                   implementation, but, if type saftey is not desired, you can just return True and be careful).
	    - dot(v,w):   . . . . . . .  Returns the scalar product of v with w, after checking for compatibility.
	    - sum(v,w):   . . . . . . .  Returns a *new vector* of the same type as v and w (compatibility is checked), which is the sum
	                                   of v and w.
	    - scaled(n,v):  . . . . . .  Returns a *new vector* of the same type as v, which has the value of v scaled by n.
	+ The details of the constructor are completely up to the user.  In general, it should probably use composition to wrap a vector of
            the desired type in an object derived from this base class, with the functions above implemented
	+ The focus of this class is on implementational simplicity, requiring the minimum amount of work from the user.  However, the supported
	    operators below functions also as a list of those things that need to be implemented for any other vector class to be compliant with
	    the expectations of the codes that use this class.  One is free to skip this base class and implement those things directly and 
	    efficiently.
	+ A non-obvious special feature is that the integer 0 is interpreted by this base class as the null vector, which is safe, since the
	    arithmetical identity must exist in all vector spaces.  However, if 0 is thought of as being the null operator (why would one do
	    this, anyway?), and 0|v is written (indistinct from projection onto 0 vector), then the return of 0 should itself be interpreted as
	    vector valued.
	"""
	def __or__(self,other):
		""" Handles (v|w) and v|A, in the case where w or A is either another vector (perhaps 0) or an operator, respectivley. """
		if other is 0:  return 0
		else:
			v = self
			if isinstance(other,_operator_base):
				A = other
				if not other.op_and_vec_are_compatible(A,v):  raise error('Illegal operator on vector from different space')
				return other.back_act_on_vec(v,A)
			else:
				w = other
				if not self.vecs_are_compatible(v,w):  raise error('Illegal contraction of vectors from different spaces')
				value = self.dot(v,w)
				if   v is w and abs(value.imag/value.real)>1e14:  raise error('Something seriously wrong.  Norm^2 of vector is not real.')	# Danger, hardcoded threshold!
				elif v is w:  return value.real
				else:         return value
	def __ror__(self,other):
		""" Handles (0|v).  If the left-hand vector would be of vector class, then __or__ above would take precedence. """
		if other is 0:  return 0
		else:  raise error('Illegal contraction of vector with unknown object')
	def __add__(self,w):
		""" Handles v+w, where w is another vector (perhaps 0). """
		if w is 0:  return self.scaled(1,self)
		else:
			v = self
			if not self.vecs_are_compatible(v,w):  raise error('Illegal linear combination of vectors from different spaces')
			return self.sum(v,w)
	def __radd__(self,w):
		""" Handles 0+v.  If the left-hand vector would be of vector class, then __add__ above would take precedence. """
		if w is 0:  return self.scaled(1,self)
		else:  raise error('Illegal linear combination of vector with unknown object')
	def __sub__(self,w):
		""" Handles v-w, where w is another vector (perhaps 0). """
		if w is 0:  return self.scaled(1,self)
		else:
			v = self
			if not self.vecs_are_compatible(v,w):  raise error('Illegal linear combination of vectors from different spaces')
			return self.sum(v,self.scaled(-1,w))
	def __rsub__(self,w):
		""" Handles 0-v.  If the left-hand vector would be of vector class, then __sub__ above would take precedence. """
		if w is 0:  return self.scaled(-1,self)
		else:  raise error('Illegal linear combination of vector with unknown object')
	def __mul__(self,n):
		""" Handles v*n, where v is a vector and n is a scalar. """
		v = self
		return self.scaled(n,v)
	def __rmul__(self,n):
		""" Handles n*v, where v is a vector and n is a scalar. """
		v = self
		return self.scaled(n,v)
	def __truediv__(self,n):
		""" Handles v/n, where v is a vector and n is a scalar. """
		v = self
		return self.scaled(1./n,v)



class _composite_op(_operator_base):
	"""\
	+ Class for holding added or contracted (multiplied) operators X (e.g., X=A|B or X=A+B), so that they have X|v
	    and v|X on a vector v defined.
	+ Maps X|v and v|X of composite X in terms of + and | for the stored primitive operators (see operator class for details).
	+ Adding and contracting of *other composites* is defined here.
	"""
	def __init__(self,A,B,algebra):
		self.A       = A
		self.B       = B
		self.algebra = algebra
		# compatibility is assumed already checked, so in principle that includes the field too, though this is getting
                # increasingly hairy (see note in vector_set.py)
		self.field = A.field
	def __or__(self,other):
		if isinstance(other,_operator_base):
			A = self
			B = other
			if not self.ops_are_compatible(A,B):  raise error('Illegal contraction of operators on different spaces')
			return _composite_op(A,B,'|')
		else:
			A = self.A
			B = self.B
			v = other
			if not self.op_and_vec_are_compatible(self,v):  raise error('Illegal operator on vector from different space')
			if self.algebra=='+':  return (A|v) + (B|v)
			if self.algebra=='-':  return (A|v) - (B|v)
			if self.algebra=='|':  return (A|(B|v))
			if self.algebra=='*':  return (A|v) * B		# B is a scalar
			if self.algebra=='/':  return (A|v) / B		# B is a scalar
	def __add__(self,B):
		A = self
		if not self.ops_are_compatible(A,B):  raise error('Illegal linear combination of operators on different spaces')
		return _composite_op(A,B,'+')
	def __sub__(self,B):
		A = self
		if not self.ops_are_compatible(A,B):  raise error('Illegal linear combination of operators on different spaces')
		return _composite_op(A,B,'-')
	def __mul__(self,number):
		A = self
		if type(number) not in (float,int,complex):  raise error('Cannot multiply operator by non-scalar')
		return _composite_op(A,number,'*')
	def __rmul__(self,number):
		A = self
		if type(number) not in (float,int,complex):  raise error('Cannot multiply operator by non-scalar')
		return _composite_op(A,number,'*')
	def __neg__(self):
		A = self
		return _composite_op(A,-1,'*')
	def __truediv__(self,number):
		A = self
		if type(number) not in (float,int,complex):  raise error('Cannot divide operator by non-scalar')
		return _composite_op(A,number,'/')
	@staticmethod
	def back_act_on_vec(v,C):
		A = C.A
		B = C.B
		if C.algebra=='+':  return (v|A) + (v|B)
		if C.algebra=='-':  return (v|A) - (v|B)
		if C.algebra=='|':  return ((v|A)|B)
		if C.algebra=='*':  return (v|A) * B	# B is a scalar
		if C.algebra=='/':  return (v|A) / B	# B is a scalar
	@staticmethod
	def op_and_vec_are_compatible(C,v):
		B = C.B
		return B.op_and_vec_are_compatible(B,v)
	@staticmethod
	def ops_are_compatible(C,D):
		B = C.B
		if isinstance(D,_composite_op):  return D.ops_are_compatible(B,D)
		else:                           return B.ops_are_compatible(B,D)



class operator(_operator_base):
	"""\
	+ Universal base class for user-implemented operator classes.
	+ Deriving from this class maps the actions A|v and v|A to static virtual functions that the user implements in the child class.
	+ A child class should implement the functions:
	    - ops_are_compatible(A,B):  . . .  Returns a boolean; True if A and B act on the same vector space (user creativity needed
	                                         for implementation, but, if type saftey is not desired, you can just return True and be
	                                         careful).  The implementor may assume that both A and be are non-composite operators;
                                                 the recursion needed to deal with testing composites is taken care of in this file
	    - op_and_vec_are_compatible(A,v):  Returns a boolean; True if A acts on the vector space to which v belongs (user creativity
                                                 needed for implementation, but, if type saftey is not desired, you can just return True
                                                 and be careful).
	    - act_on_vec(A,v):  . . . . . . .  Returns a *new vector* of the same type as v, which is the action of A onto v from the left
	                                         (compatibility is checked).
	    - back_act_on_vec(v,A):   . . . .  Returns a *new vector* of the same type as v, which is the action of A onto v from the right
                                                 (compatibility is checked).  Note that this function is not called from the operator
	                                         class but rather it is called from the vector class when something like v|A is parsed.
	### Replace the above with something more like add_to so that v+=A|w is all one operation? ... of course, the functions already defined below retain their behavior.
	+ Adding and contracting of primitive operators (returns a composite) is defined here.
	"""
	def __or__(self,other):
		A = self
		if isinstance(other,_operator_base):
			B = other
			if not self.ops_are_compatible(A,B):  raise error('Illegal contraction of operators on different spaces')
			return _composite_op(A,B,'|')
		else:
			v = other
			if not self.op_and_vec_are_compatible(A,v):  raise error('Illegal operator on vector from different space')
			return self.act_on_vec(A,v)
	def __add__(self,B):
		A = self
		if not self.ops_are_compatible(A,B):  raise error('Illegal linear combination of operators on different spaces')
		return _composite_op(A,B,'+')
	def __sub__(self,B):
		A = self
		if not self.ops_are_compatible(A,B):  raise error('Illegal linear combination of operators on different spaces')
		return _composite_op(A,B,'-')
	def __mul__(self,number):
		A = self
		if type(number) not in (float,int,complex):  raise error('Cannot multiply operator by non-scalar')
		return _composite_op(A,number,'*')
	def __rmul__(self,number):
		A = self
		if type(number) not in (float,int,complex):  raise error('Cannot multiply operator by non-scalar')
		return _composite_op(A,number,'*')
	def __neg__(self):
		A = self
		return _composite_op(A,-1,'*')
	def __truediv__(self,number):
		A = self
		if type(number) not in (float,int,complex):  raise error('Cannot divide operator by non-scalar')
		return _composite_op(A,number,'/')



class functionable_op(operator):
	"""\
	+ Allows one to write things like 2*A and exp(A), assuming the implementor can implement member functions
	    - copy
	    - function_on_diags
	+ These sort of assume the user is storing the operator in some diagonalized form, but it is not strictly necessary
	"""
	def __rmul__(self,number):
		func = lambda x: number*x 
		cp = self.copy()
		cp.function_on_diags(func)
		return cp
	def __mul__(self,number):
		func = lambda x: number*x 
		cp = self.copy()
		cp.function_on_diags(func)
		return cp
	def __truediv__(self,number):
		func = lambda x: x/number 
		cp = self.copy()
		cp.function_on_diags(func)
		return cp
	def __neg__(self):
		func = lambda x: -x 
		cp = self.copy()
		cp.function_on_diags(func)
		return cp
	def __pow__(self,number):
		func = lambda x: x**number 
		cp = self.copy()
		cp.function_on_diags(func)
		return cp



def abs(obj):
	""" Overrides the abs() function to seamlessly handle all numbers and allowed (functionable) operators. """
	if isinstance(obj,functionable_op):
		cp = obj.copy()
		cp.function_on_diags(field_traits.abs)
		return cp
	else:   return field_traits.abs(obj)

def conjugate(obj):
	""" Overrides the conjugate() function to seamlessly handle all numbers and allowed (functionable) operators. """
	if isinstance(obj,functionable_op):
		cp = obj.copy()
		cp.function_on_diags(field_traits.conjugate)
		return cp
	else:   return field_traits.conjugate(obj)

def sqrt(obj):
	""" Overrides the sqrt() function to seamlessly handle all numbers and allowed (functionable) operators. """
	if isinstance(obj,functionable_op):
		cp = obj.copy()
		cp.function_on_diags(field_traits.sqrt)
		return cp
	else:   return field_traits.sqrt(obj)

def exp(obj):
	""" Overrides the exp() function to seamlessly handle all numbers and allowed (functionable) operators. """
	if isinstance(obj,functionable_op):
		cp = obj.copy()
		cp.function_on_diags(field_traits.exp)
		return cp
	else:   return field_traits.exp(obj)
