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
from ..many_body.hierarchical_fluctuations.nested_operator import ZeroCommutator
from .state import create, annihilate, op_string

class excitation(object):
	def __init__(self,a,i):
		self.E_transition = (a,i)
		self.index = (a,i)
		self.type = "Ex"
		self.is_excitation    = True
		self.is_deexcitation  = False
		self.is_occ_rearrange = False
		self.is_vrt_rearrange = False
		self.is_flat          = False
		self.CrtOcc = []
		self.DstOcc = [i]
		self.CrtVrt = [a]
		self.DstVrt = []
	def __call__(self,the_state):
		return op_string(*self._ops())(the_state)
	def __str__(self):
		a,i = self.E_transition
		return "{}* {}  ".format(a,i)
	def _ops(self):
		a,i = self.E_transition
		return [create(a),annihilate(i)]

class deexcitation(object):
	def __init__(self,a,i):
		self.D_transition = (a,i)
		self.index = (a,i)
		self.type = "Dx"
		self.is_excitation    = False
		self.is_deexcitation  = True
		self.is_occ_rearrange = False
		self.is_vrt_rearrange = False
		self.is_flat          = False
		self.CrtOcc = [i]
		self.DstOcc = []
		self.CrtVrt = []
		self.DstVrt = [a]
	def __call__(self,the_state):
		return op_string(*self._ops())(the_state)
	def __str__(self):
		a,i = self.D_transition
		return "{}  {}* ".format(a,i)
	def _ops(self):
		a,i = self.D_transition
		return [annihilate(a),create(i)]

class occ_rearrange(object):
	def __init__(self,i,j):
		self.Fo_transition = (i,j)
		self.index = (i,j)
		self.type = "Fo"
		self.is_excitation    = False
		self.is_deexcitation  = False
		self.is_occ_rearrange = True
		self.is_vrt_rearrange = False
		self.is_flat          = True
		self.CrtOcc = [j]
		self.DstOcc = [i]
		self.CrtVrt = []
		self.DstVrt = []
	def __call__(self,the_state):
		return op_string(*self._ops())(the_state)
	def __str__(self):
		i,j = self.Fo_transition
		return "{}  {}* ".format(i,j)
	def _ops(self):
		i,j = self.Fo_transition
		return [annihilate(i),create(j)]

class vrt_rearrange(object):
	def __init__(self,a,b):
		self.Fv_transition = (a,b)
		self.index = (a,b)
		self.type = "Fv"
		self.is_excitation    = False
		self.is_deexcitation  = False
		self.is_occ_rearrange = False
		self.is_vrt_rearrange = True
		self.is_flat          = True
		self.CrtOcc = []
		self.DstOcc = []
		self.CrtVrt = [a]
		self.DstVrt = [b]
	def __call__(self,the_state):
		return op_string(*self._ops())(the_state)
	def __str__(self):
		a,b = self.Fv_transition
		return "{}* {}  ".format(a,b)
	def _ops(self):
		a,b = self.Fv_transition
		return [create(a),annihilate(b)]



def reorder_CpAq(virtuals):
	def reorder(p,q,phase=1):
		if q in virtuals:
			a = q
			if p in virtuals:
				b = p
				# b+ a
				return [[vrt_rearrange(b,a), +phase]]
			else:
				i = p
				# i+ a
				return [[deexcitation(a,i), -phase]]
		else:
			i = q
			if p in virtuals:
				a = p
				# a+ i
				return [[excitation(a,i), +phase]]
			else:
				j = p
				# j+ i
				ordered_ops = [[occ_rearrange(i,j), -phase]]
				if i==j:  ordered_ops += [[+phase]]
				return ordered_ops
	return reorder

def reorder_CpCqArAs(virtuals):
	def reorder(p,q,r,s,phase=1):
		if p==q or r==s:  return [[0]]
		if p<q:
			p,q = q,p
			phase *= -1
		if r<s:
			r,s = s,r
			phase *= -1
		if s in virtuals:
			a = s
			if r in virtuals:
				b = r
				if q in virtuals:
					c = q
					if p in virtuals:
						d = p
						# d+ c+ b  a
						ordered_ops = [[vrt_rearrange(d,b), vrt_rearrange(c,a), -phase]]
						if b==c:  ordered_ops += [[vrt_rearrange(d,a), +phase]]
						return ordered_ops
					else:
						i = p
						# i+ c+ b  a
						return [[vrt_rearrange(c,b), deexcitation(a,i), -phase]]
				else:
					i = q
					if p in virtuals:
						c = p
						# c+ i+ b  a
						return [[vrt_rearrange(c,b), deexcitation(a,i), +phase]]
					else:
						j = p
						# j+ i+ b  a
						return [[deexcitation(b,j), deexcitation(a,i), -phase]]
			else:
				i = r
				if q in virtuals:
					b = q
					if p in virtuals:
						c = p
						# c+ b+ i  a
						return [[excitation(c,i), vrt_rearrange(b,a), -phase]]
					else:
						j = p
						# j+ b+ i  a
						ordered_ops = [[excitation(b,i), deexcitation(a,j), -phase]]
						if i==j:  ordered_ops += [[vrt_rearrange(b,a), -phase]]
						return ordered_ops
				else:
					j = q
					if p in virtuals:
						b = p
						# b+ j+ i  a
						ordered_ops = [[excitation(b,i), deexcitation(a,j), +phase]]
						if i==j:  ordered_ops += [[vrt_rearrange(b,a), +phase]]
						return ordered_ops
					else:
						k = p
						# k+ j+ i  a
						ordered_ops = [[occ_rearrange(i,k), deexcitation(a,j), -phase]]
						if i==j:  ordered_ops += [[deexcitation(a,k), -phase]]
						if i==k:  ordered_ops += [[deexcitation(a,j), +phase]]
						return ordered_ops
		else:
			i = s
			if r in virtuals:
				a = r
				if q in virtuals:
					b = q
					if p in virtuals:
						c = p
						# c+ b+ a  i
						return [[excitation(c,i), vrt_rearrange(b,a), +phase]]
					else:
						j = p
						# j+ b+ a  i
						ordered_ops = [[excitation(b,i), deexcitation(a,j), +phase]]
						if i==j:  ordered_ops += [[vrt_rearrange(b,a), +phase]]
						return ordered_ops
				else:
					j = q
					if p in virtuals:
						b = p
						# b+ j+ a  i
						ordered_ops = [[excitation(b,i), deexcitation(a,j), -phase]]
						if i==j:  ordered_ops += [[vrt_rearrange(b,a), -phase]]
						return ordered_ops
					else:
						k = p
						# k+ j+ a  i
						ordered_ops = [[occ_rearrange(i,k), deexcitation(a,j), +phase]]
						if i==j:  ordered_ops += [[deexcitation(a,k), +phase]]
						if i==k:  ordered_ops += [[deexcitation(a,j), -phase]]
						return ordered_ops
			else:
				j = r
				if q in virtuals:
					a = q
					if p in virtuals:
						b = p
						# b+ a+ j  i
						return [[excitation(b,j), excitation(a,i), -phase]]
					else:
						k = p
						# k+ a+ j  i
						ordered_ops = [[excitation(a,j), occ_rearrange(i,k), -phase]]
						if i==k:  ordered_ops += [[excitation(a,j), +phase]]
						if j==k:  ordered_ops += [[excitation(a,i), -phase]]
						return ordered_ops
				else:
					k = q
					if p in virtuals:
						a = p
						# a+ k+ j  i
						ordered_ops = [[excitation(a,j), occ_rearrange(i,k), +phase]]
						if i==k:  ordered_ops += [[excitation(a,j), -phase]]
						if j==k:  ordered_ops += [[excitation(a,i), +phase]]
						return ordered_ops
					else:
						l = p
						# l+ k+ j  i
						ordered_ops = [[occ_rearrange(j,l), occ_rearrange(i,k), -phase]]
						if i==k:  ordered_ops += [[occ_rearrange(j,l), +phase]]
						if j==k:  ordered_ops += [[occ_rearrange(i,l), -phase]]
						if j==l:  ordered_ops += [[occ_rearrange(i,k), +phase]]
						if i==k and j==l:  ordered_ops += [[-phase]]
						return ordered_ops
	return reorder







def reorder_e_e(Ebj,Eai,phase=1):
	# b+ j  a+ i
	b,j = Ebj.E_transition
	a,i = Eai.E_transition
	if i==j or a==b:  return [[0]]
	if j<i:
		i,j = j,i
		phase *= -1
	if b<a:
		a,b = b,a
		phase *= -1
	return [[excitation(b,j), excitation(a,i), phase]]

def reorder_e_fo(Eck,Fji,phase=1):
	# c+ k  j  i+
	c,k = Eck.E_transition
	j,i = Fji.Fo_transition
	if j==k:  return [[0]]
	if k<j:
		j,k = k,j
		phase *= -1
	return [[excitation(c,k), occ_rearrange(j,i), phase]]

def reorder_e_fv(Eck,Fba,phase=1):
	# c+ k  b+ a
	c,k = Eck.E_transition
	b,a = Fba.Fv_transition
	if b==c:  return [[0]]
	if c<b:
		b,c = c,b
		phase *= -1
	return [[excitation(c,k), vrt_rearrange(b,a), phase]]

def reorder_e_d(Ebj,Dai,phase=1):
	# b+ j  a  i+
	b,j = Ebj.E_transition		# act as quick&dirty type checks
	a,i = Dai.D_transition
	return [[Ebj, Dai, phase]]

def reorder_fo_e(Fkj,Eai,phase=1):
	# k  j+ a+ i
	k,j = Fkj.Fo_transition
	a,i = Eai.E_transition
	if i==k:
		if i==j:  return [[excitation(a,k), +phase]]
		else:     return [[0]]
	else:
		ordered_ops = []
		if i==j:  ordered_ops += [[excitation(a,k), +phase]]
		phase *= -1
		if k<i:
			i,k = k,i
			phase *= -1
		ordered_ops += [[excitation(a,k), occ_rearrange(i,j), phase]]
		return ordered_ops

def reorder_fo_fo(Flk,Fji,phase=1):
	# l  k+ j  i+
	l,k = Flk.Fo_transition
	j,i = Fji.Fo_transition
	if j==l or i==k:
		if j==k:  return [[occ_rearrange(l,i), +phase]]
		else:     return [[0]]
	else:
		ordered_ops = []
		if l<j or k<i:		# could skip this test but would generate two cancelling terms if already ordered correctly
			if j==k:  ordered_ops += [[occ_rearrange(l,i), +phase]]
			if l<j:
				j,l = l,j
				phase *= -1
			if k<i:
				i,k = k,i
				phase *= -1
			if j==k:  ordered_ops += [[occ_rearrange(l,i), -phase]]		# the new, swapped j,k,l,i  (must be false if j==k test above was true)
		ordered_ops += [[occ_rearrange(l,k), occ_rearrange(j,i), phase]]
		return ordered_ops

def reorder_fo_fv(Fji,Fba,phase=1):
	# j  i+ b+ a
	j,i = Fji.Fo_transition
	b,a = Fba.Fv_transition
	return [[excitation(b,j), deexcitation(a,i), -phase]]

def reorder_fo_d(Fkj,Dai,phase=1):
	# k  j+ a  i+
	k,j = Fkj.Fo_transition
	a,i = Dai.D_transition
	if i==j:  return [[0]]
	if j<i:
		i,j = j,i
		phase *= -1
	return [[occ_rearrange(k,j), deexcitation(a,i), phase]]

def reorder_fv_e(Fcb,Eai,phase=1):
	# c+ b  a+ i
	c,b = Fcb.Fv_transition
	a,i = Eai.E_transition
	if a==c:
		if a==b:  return [[excitation(c,i), +phase]]
		else:     return [[0]]
	else:
		ordered_ops = []
		if a==b:  ordered_ops += [[excitation(c,i), +phase]]
		phase *= -1
		if c<a:
			a,c = c,a
			phase *= -1
		ordered_ops += [[excitation(c,i), vrt_rearrange(a,b), phase]]
		return ordered_ops

def reorder_fv_fo(Fba,Fji,phase=1):
	# b+ a  j  i+
	b,a = Fba.Fv_transition
	j,i = Fji.Fo_transition
	return [[excitation(b,j), deexcitation(a,i), -phase]]

def reorder_fv_fv(Fdc,Fba,phase=1):
	# d+ c  b+ a
	d,c = Fdc.Fv_transition
	b,a = Fba.Fv_transition
	if b==d or a==c:
		if b==c:  return [[vrt_rearrange(d,a), +phase]]
		else:     return [[0]]
	else:
		ordered_ops = []
		if d<b or c<a:		# could skip this test but would generate two cancelling terms if already ordered correctly
			if b==c:  ordered_ops += [[vrt_rearrange(d,a), +phase]]
			if d<b:
				b,d = d,b
				phase *= -1
			if c<a:
				a,c = c,a
				phase *= -1
			if b==c:  ordered_ops += [[vrt_rearrange(d,a), -phase]]		# the new, swapped b,c,d,a  (must be false if b==c test above was true)
		ordered_ops += [[vrt_rearrange(d,c), vrt_rearrange(b,a), phase]]
		return ordered_ops
		
def reorder_fv_d(Fcb,Fai,phase=1):
	# c+ b  a  i+
	c,b = Fcb.Fv_transition
	a,i = Fai.D_transition
	if a==b:  return [[0]]
	if b<a:
		a,b = b,a
		phase *= -1
	return [[vrt_rearrange(c,b), deexcitation(a,i), phase]]

def reorder_d_e(Dbj,Eai,phase=1):
	# b  j+ a+ i
	b,j = Dbj.D_transition
	a,i = Eai.E_transition
	ordered_ops = [[excitation(a,i), deexcitation(b,j), +phase]]
	if i==j:  ordered_ops += [[vrt_rearrange(a,b), +phase]]
	if a==b:  ordered_ops += [[occ_rearrange(i,j), +phase]]
	if i==j and a==b:  ordered_ops += [[-phase]]
	return ordered_ops

def reorder_d_fo(Dck,Fji,phase=1):
	# c  k+ j  i+
	c,k = Dck.D_transition
	j,i = Fji.Fo_transition
	if i==k:
		if j==k:  return [[deexcitation(c,i), +phase]]
		else:     return [[0]]
	else:
		ordered_ops = []
		if j==k:  ordered_ops += [[deexcitation(c,i), +phase]]
		phase *= -1
		if k<i:
			i,k = k,i
			phase *= -1
		ordered_ops += [[occ_rearrange(j,k), deexcitation(c,i), phase]]
		return ordered_ops

def reorder_d_fv(Dck,Fba,phase=1):
	# c  k+ b+ a
	c,k = Dck.D_transition
	b,a = Fba.Fv_transition
	if a==c:
		if b==c:  return [[deexcitation(a,k), +phase]]
		else:     return [[0]]
	else:
		ordered_ops = []
		if b==c:  ordered_ops += [[deexcitation(a,k), +phase]]
		phase *= -1
		if c<a:
			a,c = c,a
			phase *= -1
		ordered_ops += [[vrt_rearrange(b,c), deexcitation(a,k), phase]]
		return ordered_ops

def reorder_d_d(Dbj,Dai,phase=1):
	# b  j+ a  i+
	b,j = Dbj.D_transition
	a,i = Dai.D_transition
	if i==j or a==b:  return [[0]]
	if j<i:
		i,j = j,i
		phase *= -1
	if b<a:
		a,b = b,a
		phase *= -1
	return [[deexcitation(b,j), deexcitation(a,i), phase]]


def reorder_t_t(Tpq, Trs, phase=1):
	if   Tpq.is_excitation:
		if   Trs.is_excitation:
			return reorder_e_e(Tpq,Trs,phase)
		elif Trs.is_deexcitation:
			return reorder_e_d(Tpq,Trs,phase)
		elif Trs.is_occ_rearrange:
			return reorder_e_fo(Tpq,Trs,phase)
		else:
			return reorder_e_fv(Tpq,Trs,phase)
	elif Tpq.is_deexcitation:
		if   Trs.is_excitation:
			return reorder_d_e(Tpq,Trs,phase)
		elif Trs.is_deexcitation:
			return reorder_d_d(Tpq,Trs,phase)
		elif Trs.is_occ_rearrange:
			return reorder_d_fo(Tpq,Trs,phase)
		else:
			return reorder_d_fv(Tpq,Trs,phase)
	elif Tpq.is_occ_rearrange:
		if   Trs.is_excitation:
			return reorder_fo_e(Tpq,Trs,phase)
		elif Trs.is_deexcitation:
			return reorder_fo_d(Tpq,Trs,phase)
		elif Trs.is_occ_rearrange:
			return reorder_fo_fo(Tpq,Trs,phase)
		else:
			return reorder_fo_fv(Tpq,Trs,phase)
	else:
		if   Trs.is_excitation:
			return reorder_fv_e(Tpq,Trs,phase)
		elif Trs.is_deexcitation:
			return reorder_fv_d(Tpq,Trs,phase)
		elif Trs.is_occ_rearrange:
			return reorder_fv_fo(Tpq,Trs,phase)
		else:
			return reorder_fv_fv(Tpq,Trs,phase)



def commute_t_e(Tpq, Trs, phase=1):
	if   Trs.is_excitation:
		a,i = Trs.E_transition
		if   Tpq.is_excitation:
			# [ b+ j , a+ i ]			
			raise ZeroCommutator()
		elif Tpq.is_occ_rearrange:
			j,k = Tpq.Fo_transition
			# [ j k+ , a+ i ]
			if i==k:  return [excitation(a,j)]
			else:     raise ZeroCommutator()
		elif Tpq.is_vrt_rearrange:
			b,c = Tpq.Fv_transition
			# [ b+ c , a+ i ]
			if a==c:  return [excitation(b,i)]
			else:     raise ZeroCommutator()
		else:
			b,j = Tpq.D_transition
			# [ b j+ , a+ i ]
			if   a==b and i==j:  return [occ_rearrange(i,j), vrt_rearrange(a,b), -1]
			elif a==b:           return [occ_rearrange(i,j)]
			elif i==j:           return [vrt_rearrange(a,b)]
			else:                raise ZeroCommutator()
	else:
		raise Exception("commutator case not yet implemented")
