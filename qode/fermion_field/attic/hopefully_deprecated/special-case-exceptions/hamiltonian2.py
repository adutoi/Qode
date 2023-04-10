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
import numpy as np
from qode.fermion_field.state import state
from qode.fermion_field.nested_operator2 import nested_operator as operator
from qode.fermion_field.reorder2 import *



def _digest_V_mat(V_mat_raw, num_spin_orb):
	V_mat = np.zeros( shape=(num_spin_orb, num_spin_orb, num_spin_orb, num_spin_orb) )
	for p in range(num_spin_orb):
		for q in range(num_spin_orb):
			for r in range(num_spin_orb):
				for s in range(num_spin_orb):				# Yikes!  V_mat_raw should come in as 4-index tensor!
					V_mat[p,q,r,s]  = V_mat_raw[p * num_spin_orb + q, s * num_spin_orb + r]
					V_mat[p,q,r,s] -= V_mat_raw[p * num_spin_orb + q, r * num_spin_orb + s]
					V_mat[p,q,r,s] /= 4
	return V_mat



def hamiltonian(E_nuc, h_mat, V_mat_raw, occ_strings, occ_orbs, vrt_orbs):
	num_orb = len(occ_orbs) + len(vrt_orbs)
	V_mat = _digest_V_mat(V_mat_raw,num_orb)
	H = operator(occ_strings,occ_orbs,occ_orbs,vrt_orbs,vrt_orbs)
	reorder_H1 = reorder_CpAq(vrt_orbs)
	reorder_H2 = reorder_CpCqArAs(vrt_orbs)
	H.add_term([E_nuc])
	for p in range(num_orb):
		for q in range(num_orb):
			terms = reorder_H1(p,q,h_mat[p,q])
			for term in terms:  H.add_term(term)
	for p in range(num_orb):
		for q in range(num_orb):
			for r in range(num_orb):
				for s in range(num_orb):
					terms = reorder_H2(p,q,r,s,V_mat[p,q,r,s])
					for term in terms:  H.add_term(term)
	return H



class inv_fock(object):
	def __init__(self, occ_orbs, vrt_orbs, fock_mat, occ_strings):
		self.occ_orbs    = occ_orbs
		self.vrt_orbs    = vrt_orbs
		self.occ_strings = occ_strings
		self.fock_diag   = [ fock_mat[i,i] for i in range(fock_mat.shape[0]) ]	# could just get these from Hamiltonian
	def __call__(self, HBar):
		new_operator = operator(self.occ_strings, [], self.occ_orbs, self.vrt_orbs, [])
		self._apply_denominators(new_operator, HBar, 0)
		return new_operator
	def _apply_denominators(self, new_operator, old_operator, shift):	# Very dirty, looking at another objects privates
		if old_operator._Ex is not None:
			target = new_operator._get_Ex_block()
			for a in old_operator._CrtVrt_range:
				for i in old_operator._DstOcc_range:
					denom = self.fock_diag[a] - self.fock_diag[i] + shift
					target[a][i]._C = old_operator._Ex[a][i]._C / denom			# Dividing ... 
					self._apply_denominators(target[a][i], old_operator._Ex[a][i], denom)	# ... and then recurring misses reference, which is desired.
