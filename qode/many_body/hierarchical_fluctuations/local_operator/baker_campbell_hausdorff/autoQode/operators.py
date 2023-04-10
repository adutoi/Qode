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
from scalars import delta, zero, one, minus
from containers import add, mult
from indices import comparable



class Operator(indexed_base):
	def __init__(self, name, bottom, top, side, unwrap):
		indexed_base.__init__(self)
		self.name = name
		if unwrap and bottom is not None:  bottom = bottom.index
		if unwrap and top    is not None:  top    = top.index
		if unwrap and side   is not None:  side   = side.index
		self.bottom = bottom	# an Index object (or None)
		self.top    = top	# an Index object (or None)
		self.side   = side	# an Index object (or None) ... specificially interpreted to be an index on which other indices depend
		if (self.bottom is not None) and (self.side is not None):
			if self.bottom.block is None:
				self.bottom.block = self.side
				self.side.blocked_indices += [self.bottom]
		if (self.top    is not None) and (self.side is not None):
			if self.top.block is None:
				self.top.block    = self.side
				self.side.blocked_indices += [self.top]
	def _replace_indices(self, index_replacements):		# private helper function.  in place modification.
		if index_replacements:
			old,new = zip(*index_replacements)
			if self.bottom in old:  self.bottom = new[old.index(self.bottom)]	# this uses list.index ... unfortunate collision of naming conventions
			if self.top    in old:  self.top    = new[old.index(self.top)]		# this uses list.index ... unfortunate collision of naming conventions
			if self.side   in old:  self.side   = new[old.index(self.side)]		# this uses list.index ... unfortunate collision of naming conventions
		return self
	def __str__(self):
		# Sort of assumes this is being called from inside print(<instance of Sum>)
		string = self.name
		if self.bottom is not None:  string +=  "_{" + str(self.bottom) +   "}"
		if self.top    is not None:  string +=  "^{" + str(self.top)    +   "}"
		return string

	def commute(self,other):
		# only makes sense for number conserving operators now
		return commute(self,other)	# is it bad to do this from the base class when this function differentiates between derived classes?
	def reorder(self,other):
		# only makes sense for number conserving operators now
		return reorder_t_t(self,other)	# is it bad to do this from the base class when this function differentiates between derived classes?
	def reorder_primitives(self,others):
		# only makes sense for strings of primitives now
		# is it bad to do this from the base class when this function differentiates between derived classes?
		arguments = [self] + others
		if   len(others)==1:  return reorder1(*arguments)
		elif len(others)==3:  return reorder2(*arguments)
		else:  raise Exception("wrong number of field operators in reordering")
	@staticmethod
	def _wrap1(p, wrap):
		if wrap:  return comparable(p)
		else:     return p
	@staticmethod
	def _wrap2(p, q, wrap):
		if wrap:  return (comparable(p),comparable(q))
		else:     return (p,q)
	@staticmethod
	def _wrap3(p, q, M, wrap):
		if wrap:  return (comparable(p),comparable(q),comparable(M))
		else:     return (p,q,M)



class Ex(Operator):
	def __init__(self, a, i, M, unwrap=False):				# a* i   in order written
		Operator.__init__(self, "{\\hat{E}}", i, a, M, unwrap)		# i on bottom, a on top
		self.is_Ex     = True
		self.is_not_Dx = True
		self.is_field_op = False
	def transition(self, wrap=False):
		a,i,M = self.top, self.bottom, self.side			# a* i   in order written
		return Operator._wrap3(a,i,M,wrap)
	def _copy(self, index_replacements):
		return Ex(*self.transition())._replace_indices(index_replacements)
	def code(self):
		return [ 'Ex',  [ self.bottom.code(), self.top.code() ] ]


class Fo(Operator):
	def __init__(self, i, j, M, unwrap=False):				# i j*   in order written
		Operator.__init__(self, "{\\hat{F}}", i, j, M, unwrap)		# i on bottom, j on top
		self.is_Ex     = False
		self.is_not_Dx = True
		self.is_field_op = False
	def transition(self, wrap=False):
		i,j,M = self.bottom, self.top, self.side			# i j*   in order written
		return Operator._wrap3(i,j,M,wrap)
	def _copy(self, index_replacements):
		return Fo(*self.transition())._replace_indices(index_replacements)
	def code(self):
		return [ 'Fo', [ self.bottom.code(), self.top.code() ] ]

class Fv(Operator):
	def __init__(self, a, b, M, unwrap=False):				# a* b   in order written
		Operator.__init__(self, "{\\hat{F}}", b, a, M, unwrap)		# b on bottom, a on top
		self.is_Ex     = False
		self.is_not_Dx = True
		self.is_field_op = False
	def transition(self, wrap=False):
		a,b,M = self.top, self.bottom, self.side			# a* b   in order written
		return Operator._wrap3(a,b,M,wrap)
	def _copy(self, index_replacements):
		return Fv(*self.transition())._replace_indices(index_replacements)
	def code(self):
		return [ 'Fv',  [ self.bottom.code(), self.top.code() ] ]
		
class Dx(Operator):
	def __init__(self, a, i, M, unwrap=False):				# a i*   in order written
		Operator.__init__(self, "{\\hat{D}}", a, i, M, unwrap)		# a on bottom, i on top
		self.is_Ex     = False
		self.is_not_Dx = False
		self.is_field_op = False
	def transition(self, wrap=False):
		a,i,M = self.bottom, self.top, self.side			# a i*   in order written
		return Operator._wrap3(a,i,M,wrap)
	def _copy(self, index_replacements):
		return Dx(*self.transition())._replace_indices(index_replacements)
	def code(self):
		return [ 'Dx',  [ self.bottom.code(), self.top.code() ] ]



# I like f_p for destruction and f^p for creation.  "f" stands for field operator

class DstOcc(Operator):
	def __init__(self, i, unwrap=False):
		Operator.__init__(self, "{\\hat{f}}", i, None, None, unwrap)
		self.is_field_op = True
	def orbital(self, wrap=False):
		i = self.bottom
		return Operator._wrap1(i,wrap)
	def _copy(self, index_replacements):
		return DstOcc(self.orbital())._replace_indices(index_replacements)


class DstVrt(Operator):
	def __init__(self, a, unwrap=False):
		Operator.__init__(self, "{\\hat{f}}", a, None, None, unwrap)
		self.is_field_op = True
	def orbital(self, wrap=False):
		a = self.bottom
		return Operator._wrap1(a,wrap)
	def _copy(self, index_replacements):
		return DstVrt(self.orbital())._replace_indices(index_replacements)


class CrtOcc(Operator):
	def __init__(self, i, unwrap=False):
		Operator.__init__(self, "{\\hat{f}}", None, None, i, unwrap)
		self.is_field_op = True
	def orbital(self, wrap=False):
		i = self.top
		return Operator._wrap1(i,wrap)
	def _copy(self, index_replacements):
		return CrtOcc(self.orbital())._replace_indices(index_replacements)

class CrtVrt(Operator):
	def __init__(self, a, unwrap=False):
		Operator.__init__(self, "{\\hat{f}}", None, None, a, unwrap)
		self.is_field_op = True
	def orbital(self, wrap=False):
		a = self.top
		return Operator._wrap1(a,wrap)
	def _copy(self, index_replacements):
		return CrtVrt(self.orbital())._replace_indices(index_replacements)




def commute_Ex_Ex(Ex1, Ex2):
	return zero

def commute_Fo_Ex(Fo1, Ex2):
	i,j,M = Fo1.transition()
	c,k,N = Ex2.transition()
	#   F{^j}{_i}   E{^c}{_k}
	return mult(delta(j,k), Ex(c,i,M) )

def commute_Fv_Ex(Fv1, Ex2):
	a,b,M = Fv1.transition()
	c,k,N = Ex2.transition()
	#   F{^a}{_b}   E{^c}{_k}
	return mult(delta(b,c), Ex(a,k,M) )

def commute_Dx_Ex(Dx1, Ex2):
	a,i,M = Dx1.transition()
	b,j,N = Ex2.transition()
	#   D{^i}{_a}   E{^b}{_j}
	return add( mult(delta(a,b), Fo(j,i,M)), mult(delta(i,j), Fv(b,a,M)), mult(minus, delta(i,j), delta(a,b)) )

def commute(T1,T2):
	if isinstance(T2,Ex):
		if   isinstance(T1,Ex):  return commute_Ex_Ex(T1,T2)
		elif isinstance(T1,Fo):  return commute_Fo_Ex(T1,T2)
		elif isinstance(T1,Fv):  return commute_Fv_Ex(T1,T2)
		elif isinstance(T1,Dx):  return commute_Dx_Ex(T1,T2)
		else:  raise Exception("commutation of non-primitives ... did you forget to .expand()")
	else:  raise Exception("commutation of non-primitives ... did you forget to .expand()")








def reorder_e_e(Ebj,Eai):
	# b+ j  a+ i
	b,j,M = Ebj.transition(wrap=True)
	a,i,N = Eai.transition(wrap=True)
	if M==N:
		if i==j or a==b:  return zero, True
		phase = 1
		reordered = False
		if j<i:
			i,j = j,i
			phase *= -1
			reordered = True
		if b<a:
			a,b = b,a
			phase *= -1
			reordered = True
		if phase==-1:  return mult(minus, Ex(b,j,M,unwrap=True), Ex(a,i,M,unwrap=True)), reordered
		else:          return mult(       Ex(b,j,M,unwrap=True), Ex(a,i,M,unwrap=True)), reordered
	else:
		if M>N:  return mult(Ebj,Eai), False
		else:    return mult(Eai,Ebj), True

def reorder_e_fo(Eck,Fji):
	# c+ k  j  i+
	c,k,M = Eck.transition(wrap=True)
	j,i,N = Fji.transition(wrap=True)
	if M==N:
		if j==k:  return zero, True
		phase = 1
		if k<j:
			j,k = k,j
			phase *= -1
		if phase==-1:  return mult(minus, Ex(c,k,M,unwrap=True), Fo(j,i,M,unwrap=True)), True
		else:          return mult(       Ex(c,k,M,unwrap=True), Fo(j,i,M,unwrap=True)), False
	else:
		return mult(Eck,Fji), False

def reorder_e_fv(Eck,Fba):
	# c+ k  b+ a
	c,k,M = Eck.transition(wrap=True)
	b,a,N = Fba.transition(wrap=True)
	if M==N:
		if b==c:  return zero, True
		phase = 1
		if c<b:
			b,c = c,b
			phase *= -1
		if phase==-1:  return mult(minus, Ex(c,k,M,unwrap=True), Fv(b,a,M,unwrap=True)), True
		else:          return mult(       Ex(c,k,M,unwrap=True), Fv(b,a,M,unwrap=True)), False
	else:
		return mult(Eck,Fba), False

def reorder_e_d(Ebj,Dai):
	# b+ j  a  i+
	return mult(Ebj, Dai), False

def reorder_fo_e(Fkj,Eai):
	# k  j+ a+ i
	k,j,M = Fkj.transition(wrap=True)
	a,i,N = Eai.transition(wrap=True)
	if M==N:
		if i==k:
			if i==j:  return Ex(a,k,M,unwrap=True), True
			else:     return zero, True
		else:
			phase = 1
			ordered_ops = []
			if i==j:  ordered_ops += [Ex(a,k,M,unwrap=True)]
			phase *= -1
			if k<i:
				i,k = k,i
				phase *= -1
			if phase==-1:  ordered_ops += [mult(minus, Ex(a,k,M,unwrap=True), Fo(i,j,M,unwrap=True))]
			else:          ordered_ops += [mult(       Ex(a,k,M,unwrap=True), Fo(i,j,M,unwrap=True))]
			return add(*ordered_ops), True
	else:
		return mult(Eai,Fkj), True

def reorder_fo_fo(Flk,Fji):
	# l  k+ j  i+
	l,k,M = Flk.transition(wrap=True)
	j,i,N = Fji.transition(wrap=True)
	if M==N:
		if j==l or i==k:
			if j==k:  return Fo(l,i,M,unwrap=True), True
			else:     return zero, True
		else:
			phase = 1
			reordered = False
			ordered_ops = []
			if l<j or k<i:		# could skip this test but would generate two cancelling terms if already ordered correctly
				reordered = True
				if j==k:  ordered_ops += [Fo(l,i,M,unwrap=True)]
				if l<j:
					j,l = l,j
					phase *= -1
				if k<i:
					i,k = k,i
					phase *= -1
				if j==k:
					if phase==-1:  ordered_ops += [            Fo(l,i,M,unwrap=True) ]		# the new, swapped j,k,l,i  (must be false if j==k test above was true)
					else:          ordered_ops += [mult(minus, Fo(l,i,M,unwrap=True))]		# the new, swapped j,k,l,i  (must be false if j==k test above was true)
			if phase==-1:  ordered_ops += [mult(minus, Fo(l,k,M,unwrap=True), Fo(j,i,M,unwrap=True))]
			else:          ordered_ops += [mult(       Fo(l,k,M,unwrap=True), Fo(j,i,M,unwrap=True))]
			return add(*ordered_ops), reordered
	else:
		if M>N:  return mult(Flk,Fji), False
		else:    return mult(Fji,Flk), True

def reorder_fo_fv(Fji,Fba):
	# j  i+ b+ a
	j,i,M = Fji.transition(wrap=True)
	b,a,N = Fba.transition(wrap=True)
	if M==N:
		return mult(minus, Ex(b,j,M,unwrap=True), Dx(a,i,M,unwrap=True)), True
	else:
		return mult(Fji,Fba), False

def reorder_fo_d(Fkj,Dai):
	# k  j+ a  i+
	k,j,M = Fkj.transition(wrap=True)
	a,i,N = Dai.transition(wrap=True)
	if M==N:
		if i==j:  return zero, True
		phase = 1
		if j<i:
			i,j = j,i
			phase *= -1
		if phase==-1:  return mult(minus, Fo(k,j,M,unwrap=True), Dx(a,i,M,unwrap=True)), True
		else:          return mult(       Fo(k,j,M,unwrap=True), Dx(a,i,M,unwrap=True)), False
	else:
		return mult(Fkj,Dai), False

def reorder_fv_e(Fcb,Eai):
	# c+ b  a+ i
	c,b,M = Fcb.transition(wrap=True)
	a,i,N = Eai.transition(wrap=True)
	if M==N:
		if a==c:
			if a==b:  return Ex(c,i,M,unwrap=True), True
			else:     return zero, True
		else:
			phase = 1
			ordered_ops = []
			if a==b:  ordered_ops += [Ex(c,i,M,unwrap=True)]
			phase *= -1
			if c<a:
				a,c = c,a
				phase *= -1
			if phase==-1:  ordered_ops += [mult(minus, Ex(c,i,M,unwrap=True), Fv(a,b,M,unwrap=True))]
			else:          ordered_ops += [mult(       Ex(c,i,M,unwrap=True), Fv(a,b,M,unwrap=True))]
			return add(*ordered_ops), True
	else:
		return mult(Eai,Fcb), True

def reorder_fv_fo(Fba,Fji):
	# b+ a  j  i+
	b,a,M = Fba.transition(wrap=True)
	j,i,N = Fji.transition(wrap=True)
	if M==N:
		return mult(minus, Ex(b,j,M,unwrap=True), Dx(a,i,M,unwrap=True)), True
	else:
		return mult(Fji,Fba), True

def reorder_fv_fv(Fdc,Fba):
	# d+ c  b+ a
	d,c,M = Fdc.transition(wrap=True)
	b,a,N = Fba.transition(wrap=True)
	if M==N:
		if b==d or a==c:
			if b==c:  return Fv(d,a,M,unwrap=True), True
			else:     return zero, True
		else:
			phase = 1
			reordered = False
			ordered_ops = []
			if d<b or c<a:		# could skip this test but would generate two cancelling terms if already ordered correctly
				reordered = True
				if b==c:  ordered_ops += [Fv(d,a,M,unwrap=True)]
				if d<b:
					b,d = d,b
					phase *= -1
				if c<a:
					a,c = c,a
					phase *= -1
				if b==c:
					if phase==-1:  ordered_ops += [            Fv(d,a,M,unwrap=True) ]		# the new, swapped b,c,d,a  (must be false if b==c test above was true)
					else:          ordered_ops += [mult(minus, Fv(d,a,M,unwrap=True))]		# the new, swapped b,c,d,a  (must be false if b==c test above was true)
			if phase==-1:  ordered_ops += [mult(minus, Fv(d,c,M,unwrap=True), Fv(b,a,M,unwrap=True))]
			else:          ordered_ops += [mult(       Fv(d,c,M,unwrap=True), Fv(b,a,M,unwrap=True))]
			return add(*ordered_ops), reordered
	else:
		if M>N:  return mult(Fdc,Fba), False
		else:    return mult(Fba,Fdc), True
		
def reorder_fv_d(Fcb,Dai):
	# c+ b  a  i+
	c,b,M = Fcb.transition(wrap=True)
	a,i,N = Dai.transition(wrap=True)
	if M==N:
		if a==b:  return zero, True
		phase = 1
		if b<a:
			a,b = b,a
			phase *= -1
		if phase==-1:  return mult(minus, Fv(c,b,M,unwrap=True), Dx(a,i,M,unwrap=True)), True
		else:          return mult(       Fv(c,b,M,unwrap=True), Dx(a,i,M,unwrap=True)), False
	else:
		return mult(Fcb,Dai), False

def reorder_d_e(Dbj,Eai):
	# b  j+ a+ i
	b,j,M = Dbj.transition(wrap=True)
	a,i,N = Eai.transition(wrap=True)
	if M==N:
		ordered_ops = [mult(Ex(a,i,M,unwrap=True), Dx(b,j,M,unwrap=True))]
		if i==j:  ordered_ops += [Fv(a,b,M,unwrap=True)]
		if a==b:  ordered_ops += [Fo(i,j,M,unwrap=True)]
		if i==j and a==b:  ordered_ops += [minus]
		return add(*ordered_ops), True
	else:
		return mult(Eai,Dbj), True

def reorder_d_fo(Dck,Fji):
	# c  k+ j  i+
	c,k,M = Dck.transition(wrap=True)
	j,i,N = Fji.transition(wrap=True)
	if M==N:
		if i==k:
			if j==k:  return Dx(c,i,M,unwrap=True), True
			else:     return zero, True
		else:
			phase = 1
			ordered_ops = []
			if j==k:  ordered_ops += [Dx(c,i,M,unwrap=True)]
			phase *= -1
			if k<i:
				i,k = k,i
				phase *= -1
			if phase==-1:  ordered_ops += [mult(minus, Fo(j,k,M,unwrap=True), Dx(c,i,M,unwrap=True))]
			else:          ordered_ops += [mult(       Fo(j,k,M,unwrap=True), Dx(c,i,M,unwrap=True))]
			return add(*ordered_ops), True
	else:
		return mult(Fji,Dck), True

def reorder_d_fv(Dck,Fba):
	# c  k+ b+ a
	c,k,M = Dck.transition(wrap=True)
	b,a,N = Fba.transition(wrap=True)
	if M==N:
		if a==c:
			if b==c:  return Dx(a,k,M,unwrap=True), True
			else:     return zero, True
		else:
			phase = 1
			ordered_ops = []
			if b==c:  ordered_ops += [Dx(a,k,M,unwrap=True)]
			phase *= -1
			if c<a:
				a,c = c,a
				phase *= -1
			if phase==-1:  ordered_ops += [mult(minus, Fv(b,c,M,unwrap=True), Dx(a,k,M,unwrap=True))]
			else:          ordered_ops += [mult(       Fv(b,c,M,unwrap=True), Dx(a,k,M,unwrap=True))]
			return add(*ordered_ops), True
	else:
		return mult(Fba,Dck), True

def reorder_d_d(Dbj,Dai):
	# b  j+ a  i+
	b,j,M = Dbj.transition(wrap=True)
	a,i,N = Dai.transition(wrap=True)
	if M==N:
		if i==j or a==b:  return zero, True
		phase = 1
		reordered = False
		if j<i:
			i,j = j,i
			phase *= -1
			reordered = True
		if b<a:
			a,b = b,a
			phase *= -1
			reordered = True
		if phase==-1:  return mult(minus, Dx(b,j,M,unwrap=True), Dx(a,i,M,unwrap=True)), reordered
		else:          return mult(       Dx(b,j,M,unwrap=True), Dx(a,i,M,unwrap=True)), reordered
	else:
		if M>N:  return mult(Dbj,Dai), False
		else:    return mult(Dai,Dbj), True



def reorder_t_t(Tpq, Trs):
	if   isinstance(Tpq,Ex):
		if   isinstance(Trs,Ex):
			return reorder_e_e(Tpq,Trs)
		elif isinstance(Trs,Dx):
			return reorder_e_d(Tpq,Trs)
		elif isinstance(Trs,Fo):
			return reorder_e_fo(Tpq,Trs)
		else:
			return reorder_e_fv(Tpq,Trs)
	elif isinstance(Tpq,Dx):
		if   isinstance(Trs,Ex):
			return reorder_d_e(Tpq,Trs)
		elif isinstance(Trs,Dx):
			return reorder_d_d(Tpq,Trs)
		elif isinstance(Trs,Fo):
			return reorder_d_fo(Tpq,Trs)
		else:
			return reorder_d_fv(Tpq,Trs)
	elif isinstance(Tpq,Fo):
		if   isinstance(Trs,Ex):
			return reorder_fo_e(Tpq,Trs)
		elif isinstance(Trs,Dx):
			return reorder_fo_d(Tpq,Trs)
		elif isinstance(Trs,Fo):
			return reorder_fo_fo(Tpq,Trs)
		else:
			return reorder_fo_fv(Tpq,Trs)
	else:
		if   isinstance(Trs,Ex):
			return reorder_fv_e(Tpq,Trs)
		elif isinstance(Trs,Dx):
			return reorder_fv_d(Tpq,Trs)
		elif isinstance(Trs,Fo):
			return reorder_fv_fo(Tpq,Trs)
		else:
			return reorder_fv_fv(Tpq,Trs)













def reorder1(Cp,Dq):
	p = Cp.orbital(wrap=True)
	q = Dq.orbital(wrap=True)
	if   isinstance(Dq,DstVrt):
		a = q
		if   isinstance(Cp,CrtVrt):
			b = p
			# b+ a
			return Fv(b,a,None,unwrap=True)
		elif isinstance(Cp,CrtOcc):
			i = p
			# i+ a
			return mult( Dx(a,i,None,unwrap=True), minus )
		else:
			raise Exception("illegal arguments to reorder")
	elif isinstance(Dq,DstOcc):
		i = q
		if   isinstance(Cp,CrtVrt):
			a = p
			# a+ i
			return Ex(a,i,None,unwrap=True)
		elif isinstance(Cp,CrtOcc):
			j = p
			# j+ i
			ordered_ops = mult( Fo(i,j,None,unwrap=True), minus )
			if i==j:  ordered_ops = add(one, ordered_ops)
			return ordered_ops
		else:
			raise Exception("illegal arguments to reorder")
	else:
		raise Exception("illegal arguments to reorder")



def reorder2(Cp,Cq,Dr,Ds):
	phase = 1
	if (isinstance(Cp,CrtVrt) and isinstance(Cq,CrtVrt)) or (isinstance(Cp,CrtOcc) and isinstance(Cq,CrtOcc)):
		p = Cp.orbital(wrap=True)
		q = Cq.orbital(wrap=True)
		if p==q:  return zero
		if p<q:
			Cp,Cq = Cq,Cp
			phase *= -1
	if (isinstance(Dr,DstVrt) and isinstance(Ds,DstVrt)) or (isinstance(Dr,DstOcc) and isinstance(Ds,DstOcc)):
		r = Dr.orbital(wrap=True)
		s = Ds.orbital(wrap=True)
		if r==s:  return zero
		if r<s:
			Dr,Ds = Ds,Dr
			phase *= -1
	if phase==-1:  return mult(minus, reorder2_dressed(Cp,Cq,Dr,Ds))
	else:          return             reorder2_dressed(Cp,Cq,Dr,Ds)


def reorder2_dressed(Cp,Cq,Dr,Ds):
	p = Cp.orbital(wrap=True)
	q = Cq.orbital(wrap=True)
	r = Dr.orbital(wrap=True)
	s = Ds.orbital(wrap=True)
	if   isinstance(Ds,DstVrt):
		a = s
		if   isinstance(Dr,DstVrt):
			b = r
			if   isinstance(Cq,CrtVrt):
				c = q
				if   isinstance(Cp,CrtVrt):
					d = p
					# d+ c+ b  a
					ordered_ops = mult( Fv(d,b,None,unwrap=True), Fv(c,a,None,unwrap=True), minus )
					if b==c:  ordered_ops = add(Fv(d,a,None,unwrap=True), ordered_ops)
					return ordered_ops
				elif isinstance(Cp,CrtOcc):
					i = p
					# i+ c+ b  a
					return mult( Fv(c,b,None,unwrap=True), Dx(a,i,None,unwrap=True), minus )
				else:
					raise Exception("illegal arguments to reorder")
			elif isinstance(Cq,CrtOcc):
				i = q
				if   isinstance(Cp,CrtVrt):
					c = p
					# c+ i+ b  a
					return mult( Fv(c,b,None,unwrap=True), Dx(a,i,None,unwrap=True) )
				elif isinstance(Cp,CrtOcc):
					j = p
					# j+ i+ b  a
					return mult( Dx(b,j,None,unwrap=True), Dx(a,i,None,unwrap=True), minus )
				else:
					raise Exception("illegal arguments to reorder")
			else:
				raise Exception("illegal arguments to reorder")
		elif isinstance(Dr,DstOcc):
			i = r
			if   isinstance(Cq,CrtVrt):
				b = q
				if   isinstance(Cp,CrtVrt):
					c = p
					# c+ b+ i  a
					return mult( Ex(c,i,None,unwrap=True), Fv(b,a,None,unwrap=True), minus )
				elif isinstance(Cp,CrtOcc):
					j = p
					# j+ b+ i  a
					ordered_ops = mult( Ex(b,i,None,unwrap=True), Dx(a,j,None,unwrap=True), minus )
					if i==j:  ordered_ops = add(mult(Fv(b,a,None,unwrap=True),minus), ordered_ops)
					return ordered_ops
				else:
					raise Exception("illegal arguments to reorder")
			elif isinstance(Cq,CrtOcc):
				j = q
				if   isinstance(Cp,CrtVrt):
					b = p
					# b+ j+ i  a
					ordered_ops = mult( Ex(b,i,None,unwrap=True), Dx(a,j,None,unwrap=True) )
					if i==j:  ordered_ops = add(Fv(b,a,None,unwrap=True), ordered_ops)
					return ordered_ops
				elif isinstance(Cp,CrtOcc):
					k = p
					# k+ j+ i  a
					ordered_ops = mult( Fo(i,k,None,unwrap=True), Dx(a,j,None,unwrap=True), minus )
					if i==j:  ordered_ops = add(mult(Dx(a,k,None,unwrap=True),minus), ordered_ops)
					if i==k:  ordered_ops = add(Dx(a,j,None,unwrap=True), ordered_ops)
					return ordered_ops
				else:
					raise Exception("illegal arguments to reorder")
			else:
				raise Exception("illegal arguments to reorder")
		else:
			raise Exception("illegal arguments to reorder")
	elif isinstance(Ds,DstOcc):
		i = s
		if   isinstance(Dr,DstVrt):
			a = r
			if   isinstance(Cq,CrtVrt):
				b = q
				if   isinstance(Cp,CrtVrt):
					c = p
					# c+ b+ a  i
					return mult( Ex(c,i,None,unwrap=True), Fv(b,a,None,unwrap=True) )
				elif isinstance(Cp,CrtOcc):
					j = p
					# j+ b+ a  i
					ordered_ops = mult( Ex(b,i,None,unwrap=True), Dx(a,j,None,unwrap=True) )
					if i==j:  ordered_ops = add(Fv(b,a,None,unwrap=True), ordered_ops)
					return ordered_ops
				else:
					raise Exception("illegal arguments to reorder")
			elif isinstance(Cq,CrtOcc):
				j = q
				if   isinstance(Cp,CrtVrt):
					b = p
					# b+ j+ a  i
					ordered_ops = mult( Ex(b,i,None,unwrap=True), Dx(a,j,None,unwrap=True), minus )
					if i==j:  ordered_ops = add(mult(Fv(b,a,None,unwrap=True),minus), ordered_ops)
					return ordered_ops
				elif isinstance(Cp,CrtOcc):
					k = p
					# k+ j+ a  i
					ordered_ops = mult( Fo(i,k,None,unwrap=True), Dx(a,j,None,unwrap=True) )
					if i==j:  ordered_ops = add(Dx(a,k,None,unwrap=True), ordered_ops)
					if i==k:  ordered_ops = add(mult(Dx(a,j,None,unwrap=True),minus), ordered_ops)
					return ordered_ops
				else:
					raise Exception("illegal arguments to reorder")
			else:
				raise Exception("illegal arguments to reorder")
		elif isinstance(Dr,DstOcc):
			j = r
			if   isinstance(Cq,CrtVrt):
				a = q
				if   isinstance(Cp,CrtVrt):
					b = p
					# b+ a+ j  i
					return mult( Ex(b,j,None,unwrap=True), Ex(a,i,None,unwrap=True), minus )
				elif isinstance(Cp,CrtOcc):
					k = p
					# k+ a+ j  i
					ordered_ops = mult( Ex(a,j,None,unwrap=True), Fo(i,k,None,unwrap=True), minus )
					if i==k:  ordered_ops = add(Ex(a,j,None,unwrap=True), ordered_ops)
					if j==k:  ordered_ops = add(mult(Ex(a,i,None,unwrap=True),minus), ordered_ops)
					return ordered_ops
				else:
					raise Exception("illegal arguments to reorder")
			elif isinstance(Cq,CrtOcc):
				k = q
				if   isinstance(Cp,CrtVrt):
					a = p
					# a+ k+ j  i
					ordered_ops = mult( Ex(a,j,None,unwrap=True), Fo(i,k,None,unwrap=True) )
					if i==k:  ordered_ops = add(mult(Ex(a,j,None,unwrap=True),minus), ordered_ops)
					if j==k:  ordered_ops = add(Ex(a,i,None,unwrap=True), ordered_ops)
					return ordered_ops
				elif isinstance(Cp,CrtOcc):
					l = p
					# l+ k+ j  i
					ordered_ops = mult( Fo(j,l,None,unwrap=True), Fo(i,k,None,unwrap=True), minus )
					if i==k:  ordered_ops = add(Fo(j,l,None,unwrap=True), ordered_ops)
					if j==k:  ordered_ops = add(mult(Fo(i,l,None,unwrap=True),minus), ordered_ops)
					if j==l:  ordered_ops = add(Fo(i,k,None,unwrap=True), ordered_ops)
					if i==k and j==l:  ordered_ops = add(minus, ordered_ops)
					return ordered_ops
				else:
					raise Exception("illegal arguments to reorder")
			else:
				raise Exception("illegal arguments to reorder")
		else:
			raise Exception("illegal arguments to reorder")
	else:
		raise Exception("illegal arguments to reorder")
