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
from .reorder import excitation, deexcitation, occ_rearrange, vrt_rearrange, reorder_t_t, commute_t_e



class transition_indices(object):
	"""\
	Just a containter for the fluctuation indices.  In principle, only the members Excitations, Flat_transitions, Deexcitations need to be defined and accessible, and these are generally
	expected from the outside to just be storage for static lists that were put in, however, manipulations are easier if we store raw orbital indices and then generate what we need.
	"""
	def __init__(self, CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range):
		self._CrtOcc_range = sorted(CrtOcc_range, reverse=True)
		self._CrtVrt_range = sorted(CrtVrt_range, reverse=True)
		self._DstOcc_range = sorted(DstOcc_range, reverse=True)
		self._DstVrt_range = sorted(DstVrt_range, reverse=True)
		self._build()
	def _build(self):
		self.Excitations      = []
		self.Flat_transitions = []
		self.Deexcitations    = []
		for a in self._CrtVrt_range:
			for i in self._DstOcc_range:
				self.Excitations      += [(a,i)]
		for a in self._CrtVrt_range:
			for b in self._DstVrt_range:
				self.Flat_transitions += [(a,b)]
		for i in self._DstOcc_range:
			for j in self._CrtOcc_range:
				self.Flat_transitions += [(i,j)]
		for a in self._DstVrt_range:
			for i in self._CrtOcc_range:
				self.Deexcitations    += [(a,i)]



class OV_partitioning(object):
	"""\
	The job of this class is to story system-specific manipulations of the indexing of fluctuations classes
	"""
	def __init__(self, occ_orbs, vrt_orbs):
		self.occ_orbs = sorted(occ_orbs, reverse=True)
		self.vrt_orbs = sorted(vrt_orbs, reverse=True)
		self.id = id(self)			# use this in conjunction with == instead of 'is' to avoid trouble with multiprocessing
	def __eq__(self, other):
		if self.id==other.id:  return True	# if running under a single process, this will be identical to (self is other), since return of id() is constant over object life
		else:                  return False
	def __neq__(self, other):
		if self==other:  return False		# uses __eq__ above
		else:            return True
	def all_transition_indices(self):
		return transition_indices(self.occ_orbs, self.occ_orbs, self.vrt_orbs, self.vrt_orbs)
	def order_indices(self, Tx_indices):
		if   isinstance(Tx_indices, transition_indices):  return Tx_indices				# this class never has the chance to become disordered
		elif isinstance(Tx_indices, list):                return sorted(Tx_indices, reverse=True)	# tuples have inherent python ordering
		else:                                             raise ValueError()
	def append_indices(self, Tx_indices, transition):
		CrtOcc_range = list( set(Tx_indices._CrtOcc_range) | set(transition.CrtOcc) )
		DstOcc_range = list( set(Tx_indices._DstOcc_range) | set(transition.DstOcc) )
		CrtVrt_range = list( set(Tx_indices._CrtVrt_range) | set(transition.CrtVrt) )
		DstVrt_range = list( set(Tx_indices._DstVrt_range) | set(transition.DstVrt) )
		return transition_indices(CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
	def combine_indices(self, Tx_indices_1, Tx_indices_2):
		CrtOcc_range = list( set(Tx_indices_1._CrtOcc_range) | set(Tx_indices_2._CrtOcc_range) )
		DstOcc_range = list( set(Tx_indices_1._DstOcc_range) | set(Tx_indices_2._DstOcc_range) )
		CrtVrt_range = list( set(Tx_indices_1._CrtVrt_range) | set(Tx_indices_2._CrtVrt_range) )
		DstVrt_range = list( set(Tx_indices_1._DstVrt_range) | set(Tx_indices_2._DstVrt_range) )
		return transition_indices(CrtOcc_range, DstOcc_range, CrtVrt_range, DstVrt_range)
	def slice_indices(self, Tx_indices, transition):
		if   transition.is_excitation:
			a,i = transition.E_transition
			na = Tx_indices._CrtVrt_range.index(a)
			ni = Tx_indices._DstOcc_range.index(i)
			return transition_indices(Tx_indices._CrtOcc_range,        Tx_indices._DstOcc_range[ni+1:], Tx_indices._CrtVrt_range[na+1:], Tx_indices._DstVrt_range)
		elif transition.is_vrt_rearrange:
			a,b = transition.Fv_transition
			na = Tx_indices._CrtVrt_range.index(a)
			nb = Tx_indices._DstVrt_range.index(b)
			return transition_indices(Tx_indices._CrtOcc_range,        [],                              Tx_indices._CrtVrt_range[na+1:], Tx_indices._DstVrt_range[nb+1:])
		elif transition.is_occ_rearrange:
			i,j = transition.Fo_transition
			ni = Tx_indices._DstOcc_range.index(i)
			nj = Tx_indices._CrtOcc_range.index(j)
			return transition_indices(Tx_indices._CrtOcc_range[nj+1:], Tx_indices._DstOcc_range[ni+1:], [],                              Tx_indices._DstVrt_range)
		elif transition.is_deexcitation:
			a,i = transition.D_transition
			na = Tx_indices._DstVrt_range.index(a)
			ni = Tx_indices._CrtOcc_range.index(i)
			return transition_indices(Tx_indices._CrtOcc_range[ni+1:], [],                              [],                              Tx_indices._DstVrt_range[na+1:])
	def excitation(self,index):
		a,i = index
		return excitation(a,i)
	def flat_transition(self,index):
		p,q = index
		if   (p in self.vrt_orbs) and (q in self.vrt_orbs):  return vrt_rearrange(p,q)
		elif (p in self.occ_orbs) and (q in self.occ_orbs):  return occ_rearrange(p,q)
		else:  raise Exception("Indices not appropriate for flat transition")
	def deexcitation(self,index):
		a,i = index
		return deexcitation(a,i)
	def make_transition(self, t_type, index):
		if   t_type=="Ex":  return self.excitation(index)
		elif t_type=="Fo":  return self.flat_transition(index)	# could just have one "Fx" case here and make .type of both above equal to "Fx"
		elif t_type=="Fv":  return self.flat_transition(index)	# but this is more readable in a dump of the operator where it is actually used.
		elif t_type=="Dx":  return self.deexcitation(index)
		else:  raise Exception("uh oh")
	@staticmethod
	def reorder(t1, t2):
		return reorder_t_t(t1, t2)
	@staticmethod
	def commute(t, e):
		if not e.is_excitation:  raise Exception("sorry, only commutation with excitations on the right is currently supported")
		return commute_t_e(t, e)
