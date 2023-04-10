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



# efficiency layer for scaling of vectors
# another thin wrapper around a vector class which is already suitable for wrapping by the below can be used to keep track of scaling
# factors "externally".  This means that the implementation of add_to as seen by this function can make use of that (perhaps factoring out a global factor)
# to speed up v+=c*w by doing this all in one step internally (externally, it still looks like two steps, but one is cheap).



class vector(object):
	"""\
	+ Universal base class for user-implemented vector classes.
	+ Deriving from this class maps the actions v|w and v+w, etc, to static virtual functions that the user implements in the child class.
	+ A child class should implement the functions:
	    - vecs_are_compatible(v,w):  Returns a boolean; True if v and w belong to the same vector space (user creativity needed for 
	                                   implementation, but, if type saftey is not desired, you can just return True and be careful).
	    - dot(v,w):   . . . . . . .  Returns the scalar product of v with w, after checking for compatibility.
	    - add_to(v,w,c):  . . . . .  In-place addition of c*w to v (c is optional! default is 1), leaving w unchanged (compatibility is
	                                   checked).  No return value needed.
	    - scale(n,v):   . . . . . .  In-place scaling of v by n.  No return value needed.
	    - copy(v,n):  . . . . . . .  Instantiates and returns an independent (deep) copy of vector v, scaled by n (n is optional! default
	                                   is 1).
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
		""" Handles (v|w) and v|A, where w or A is either another vector (perhaps 0) or an operator (returns a new vector), respectivley. """
		if other is 0:  return 0
		else:
			v = self
			if isinstance(other,operator_base):
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
		""" Handles v+w, where w is another vector (perhaps 0).  Returns a new vector. """
		v = self
		if w is 0:  return self.copy(v)
		else:
			if not self.vecs_are_compatible(v,w):  raise error('Illegal linear combination of vectors from different spaces')
			a = self.copy(v)
			self.add_to(a,w)
			return a
	def __radd__(self,w):
		""" Handles 0+v.  If the left-hand vector would be of vector class, then __add__ above would take precedence. """
		v = self
		if w is 0:  return self.copy(v)
		else:  raise error('Illegal linear combination of vector with unknown object')
	def __iadd__(self,w):
		""" Handles v+=w, where w is another vector (perhaps 0).  Returns v. """
		v = self
		if w is not 0:
			if not self.vecs_are_compatible(v,w):  raise error('Illegal linear combination of vectors from different spaces')
			self.add_to(v,w)
		return v
	def __sub__(self,w):
		""" Handles v-w, where w is another vector (perhaps 0). """
		v = self
		if w is 0:  return self.copy(v)
		else:
			if not self.vecs_are_compatible(v,w):  raise error('Illegal linear combination of vectors from different spaces')
			a = self.copy(v)
			self.add_to(a,w,-1)
			return a
	def __rsub__(self,w):
		""" Handles 0-v.  If the left-hand vector would be of vector class, then __sub__ above would take precedence. """
		v = self
		if w is 0:  return self.copy(v,-1)
		else:  raise error('Illegal linear combination of vector with unknown object')
	def __isub__(self,w):
		""" Handles v-=w, where w is another vector (perhaps 0).  Returns v. """
		v = self
		if w is not 0:
			if not self.vecs_are_compatible(v,w):  raise error('Illegal linear combination of vectors from different spaces')
			self.add_to(v,w,-1)
		return v
	def __neg__(self):
		""" Handles -v. """
		v = self
		return self.copy(v,-1)
	def __mul__(self,n):
		""" Handles v*n, where v is a vector and n is a scalar. """
		v = self
		return self.copy(v,n)
	def __rmul__(self,n):
		""" Handles n*v, where v is a vector and n is a scalar. """
		v = self
		return self.copy(v,n)
	def __truediv__(self,n):
		""" Handles v/n, where v is a vector and n is a scalar. """
		v = self
		return self.copy(v,1./n)
