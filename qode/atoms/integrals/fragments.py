#    (C) Copyright 2019 Anthony D. Dutoi
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

# WARNING:  right now, this code expects the fragments to know both the basis and the number of orbitals for that fragment
# in that basis, which is a potential place to not have synchronization.

# 8.Jun.2023:  changed name semiMO to fragMO.  Might have broken some calling code (can be fixed with 4-char substitution semi->frag)



import numpy
from ...util.dynamic_array import dynamic_array, cached, mask
from .. import nuc_repulsion
from .external_engines import psi4_ints
from . import spatial_to_spin
import timeit



def key_sorter(fragments):
    if   isinstance(fragments,list):
        sort_keys = sorted
    elif isinstance(fragments,dict):
        ordering = fragments.keys()
        def sort_keys(keys):  return sorted(keys, key=ordering.index)
    elif isinstance(fragments,dynamic_array):
        sort_keys = sorted				# this is a nasty temp hack until a better range-like object is implemented, which can also order indices (see intro in ...util.dynamic_array)
    return sort_keys

def keys(fragments):
    if   isinstance(fragments,list):
        all_keys = list(range(len(fragments)))
    elif isinstance(fragments,dict):
        all_keys = list(fragments.keys())
    elif isinstance(fragments,dynamic_array):
        all_keys = fragments.ranges[0]
    return all_keys





def _blocking(fragments, rng, spin_orbs):
    if rng is None:  rng = keys(fragments)
    spin_orbs = 2 if spin_orbs else 1
    blocks = []
    begin = 0
    for mm in rng:
        end = begin + spin_orbs*fragments[mm].basis.n_spatial_orb
        blocks += [(mm,begin,end)]
        begin = end
    return blocks
def unblock_1(psiB, fragments, ranges=None, spin_orbs=False):
    if ranges is None:  ranges = [None]
    blocks = _blocking(fragments, ranges[0], spin_orbs)
    _, _, dim = blocks[-1]
    psi = numpy.zeros(dim)
    for m,b,e in blocks:  psi[b:e] = psiB[m]
    return psi
def block_1(psi, fragments, ranges=None, spin_orbs=False):
    if ranges is None:  ranges = [None]
    blocks = _blocking(fragments, ranges[0], spin_orbs)
    psiB = {}
    for m,b,e in blocks:  psiB[m] = psi[b:e]
    return psiB
def unblock_2(Hb, fragments, ranges=None, spin_orbs=False):     # Hb could be S or T (blocked)
    if ranges is None:  ranges = [None,None]
    blocks_list = tuple(_blocking(fragments, rng, spin_orbs) for rng in ranges)
    dims = tuple(blocks[-1][-1] for blocks in blocks_list)
    H = numpy.zeros(dims)
    for m1,b1,e1 in blocks_list[0]:
        for m2,b2,e2 in blocks_list[1]:
            H[b1:e1, b2:e2] = Hb[m1,m2]
    return H
def block_2(H, fragments, ranges=None, spin_orbs=False):	# H could be S or T
    if ranges is None:  ranges = [None,None]
    blocks_list = tuple(_blocking(fragments, rng, spin_orbs) for rng in ranges)
    Hb = {}
    for m1,b1,e1 in blocks_list[0]:
        for m2,b2,e2 in blocks_list[1]:
            Hb[m1,m2] = H[b1:e1, b2:e2]
    return Hb
def unblock_last2(Ub, fragments, ranges=None, spin_orbs=False):
    if ranges is None:  ranges = [None,None,None]
    blocks_list = tuple(_blocking(fragments, rng, spin_orbs) for rng in ranges)
    dims = tuple(blocks[-1][-1] for blocks in blocks_list[1:])
    U = {}
    for m1,b1,e1 in blocks_list[0]:
        U[m1] = numpy.zeros(dims)
        for m2,b2,e2 in blocks_list[1]:
            for m3,b3,e3 in blocks_list[2]:
                U[m1][b2:e2, b3:e3] = Ub[m1,m2,m3]
    return U
def block_last2(U, fragments, ranges=None, spin_orbs=False):
    if ranges is None:  ranges = [None,None,None]
    blocks_list = tuple(_blocking(fragments, rng, spin_orbs) for rng in ranges)
    Ub = {}
    for m1,b1,e1 in blocks_list[0]:
        for m2,b2,e2 in blocks_list[1]:
            for m3,b3,e3 in blocks_list[2]:
                Ub[m1,m2,m3] = U[m1][b2:e2, b3:e3]
    return Ub
def unblock_4(Vb, fragments, ranges=None, spin_orbs=False):
    if ranges is None:  ranges = [None,None,None,None]
    blocks_list = tuple(_blocking(fragments, rng, spin_orbs) for rng in ranges)
    dims = tuple(blocks[-1][-1] for blocks in blocks_list)
    V = numpy.zeros(dims)
    for m1,b1,e1 in blocks_list[0]:
        for m2,b2,e2 in blocks_list[1]:
            for m3,b3,e3 in blocks_list[2]:
                for m4,b4,e4 in blocks_list[3]:
                    V[b1:e1, b2:e2, b3:e3, b4:e4] = Vb[m1,m2,m3,m4]
    return V
def block_4(V, fragments, ranges=None, spin_orbs=False):
    if ranges is None:  ranges = [None,None,None,None]
    blocks_list = tuple(_blocking(fragments, rng, spin_orbs) for rng in ranges)
    Vb = {}
    for m1,b1,e1 in blocks_list[0]:
        for m2,b2,e2 in blocks_list[1]:
            for m3,b3,e3 in blocks_list[2]:
                for m4,b4,e4 in blocks_list[3]:
                    Vb[m1,m2,m3,m4] = V[b1:e1, b2:e2, b3:e3, b4:e4]
    return Vb







class Nuc_repulsion(object):
    def __init__(self, fragments, rule_wrappers=None):
        self.fragments = fragments
        frag_ids = keys(self.fragments)
        sort_keys = key_sorter(fragments)
        if rule_wrappers is None:  rule_wrappers = []
        frag_calcs = dynamic_array([cached, Nuc_repulsion._compute_rule(fragments, sort_keys)], ranges=[None])	# doesn't hurt to cache (small data, each inefficiently computed)
        self.matrix = dynamic_array(rule_wrappers + [Nuc_repulsion._resolve_rule(frag_calcs,sort_keys)], ranges=[frag_ids,frag_ids])	# diagonal elements give internal repulsion energy of fragment, off-diagonals give only inter-fragment repulsion
    def total_repulsion_energy(self):
        frag_ids = keys(self.fragments)
        total = 0
        for i,m1 in enumerate(frag_ids):
            for m2 in frag_ids[i:]:
                total += self.matrix[m1,m2]
        return total
    @staticmethod
    def _compute_rule(fragments, sort_keys):
        def nuc_energy(frag_indices):
            if sort_keys(frag_indices)!=sort_keys(set(frag_indices)):  raise AssertionError    # The indices should all be unique (but in any order the user likes)
            atoms = []
            for i in frag_indices:  atoms += fragments[i].atoms
            return nuc_repulsion.energy(atoms)
        return nuc_energy
    @staticmethod
    def _resolve_rule(all_repulsions, sort_keys):
        def get_repulsion(m1,m2):
            if m1==m2:
                return all_repulsions[((m1,),)]
            else:
                m1,m2 = tuple(sort_keys((m1,m2)))	# avoid redundant calculations; unpacking ensures only called for dimers
                return all_repulsions[((m1,m2),)] - all_repulsions[((m1,),)] - all_repulsions[((m2,),)]
        return get_repulsion



class AO_integrals(object):
    def __init__(self, fragments, integrals_engine=psi4_ints, rule_wrappers=None):
        self.fragments = fragments
        frag_ids = keys(self.fragments)
        sort_keys = key_sorter(fragments)
        if rule_wrappers is None:  rule_wrappers = []
        frag_calcs = dynamic_array([cached, AO_integrals._compute_rule(fragments, sort_keys, integrals_engine)], ranges=[None])	# definitely cache because so expensive and called repeatedly by arrays below ... no rule_wrappers because elements different in structure than those below
        self.S   = dynamic_array(rule_wrappers + [AO_integrals._block_rule("S", frag_calcs, sort_keys)], ranges=[frag_ids,frag_ids])		# caching not important because just block finder (and each grabbed only once to do fragMO transform)
        self.T   = dynamic_array(rule_wrappers + [AO_integrals._block_rule("T", frag_calcs, sort_keys)], ranges=[frag_ids,frag_ids])
        self.U   = dynamic_array(rule_wrappers + [AO_integrals._block_rule("U", frag_calcs, sort_keys)], ranges=[frag_ids,frag_ids,frag_ids])
        self.V   = dynamic_array(rule_wrappers + [AO_integrals._block_rule("V", frag_calcs, sort_keys)], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
    @staticmethod
    def _compute_rule(fragments, sort_keys, integrals_engine):
        def the_integrals(frag_indices):

            start_time = timeit.default_timer()

            if sort_keys(frag_indices)!=sort_keys(set(frag_indices)):  raise AssertionError    # The indices should all be unique (but in any order the user likes)
            subsystem = [fragments[i] for i in frag_indices]
            U = []    # multiple distinct matrices, one for each nucleus
            for i in range(len(subsystem)):    # runs over the attraction to nuclei of individual fragments
                basis = None
                geometry = ""
                for j,frag in enumerate(subsystem):
                    if basis is None:
                        basis = frag.basis.AOcode
                    else:
                        if frag.basis.AOcode!=basis:  raise AssertionError("right now, all fragments must have the same basis")	# see also, warning at top of file
                    geometry += AO_integrals._add_fragment(frag.atoms, ghost=(i!=j))
                S, T, Ui, V, _ = integrals_engine.AO_ints(geometry, basis, NucPotentialOnly=True)
                U += [Ui]

            elapsed = timeit.default_timer() - start_time
            print('elapsed time (after integrals loops) =', elapsed, flush=True)

            S, T, _, V, _ = integrals_engine.AO_ints(geometry, basis)	# relies on most recent geometry from loop above
            return {"S":block_2(S,subsystem), "T":block_2(T,subsystem), "U":block_last2(U,subsystem), "V":block_4(V,subsystem)}

            elapsed = timeit.default_timer() - start_time
            print('elapsed time (after second int call) =', elapsed, flush=True)

        return the_integrals
    @staticmethod
    def _add_fragment(atoms, ghost=False):
        ghost = "@" if ghost else ""
        fragment = ""
        for atom in atoms:
            element, (x,y,z) = atom  #("element position")    # atom is a struct
            fragment += "{:s} {:.16f} {:.16f} {:.16f}\n".format(ghost+element, x, y, z)
        return fragment
    @staticmethod
    def _block_rule(ints_type, all_integrals, sort_keys):
        def get_block(*frag_indices):
            ordered = sort_keys(set(frag_indices))           # "ascending" order with duplicates removed
            mapping = {m:i for i,m in enumerate(ordered)}    # where a given value from indices occurs in the ordered list
            mapped  = [mapping[m] for m in frag_indices]     # absolute fragment indices as given, mapped to relative indices of subsystem in standard order
            ordered,mapped =  tuple(ordered),tuple(mapped)   # ordered is suitable for passing to the integrals engine, and mapped helps find the desired block
            return all_integrals[(ordered,)][ints_type][mapped]   # all_integrals is a "1-D" array whose index is a variable-length tuple
        return get_block



class fragMO_integrals(object):		# actually more general than MOs ... any (fragment-local) transformations you like, including right not the same as left ... would like to let either right or left or individual transforms be 1, but then checking is more involved, so let sleeping dogs lie for now
    def __init__(self, AOints, right=True, left=True, cache=False, rule_wrappers=None):
        self.fragments = AOints.fragments
        frag_ids = keys(self.fragments)
        if   (left is True or left is False) and (right is True or right is False):  raise AssertionError
        elif left  is True:  left  = right
        elif right is True:  right = left
        cache = [] if cache else [cached]
        if rule_wrappers is None:  rule_wrappers = []
        self.S = dynamic_array(cache + rule_wrappers + [fragMO_integrals._transform2_rule(    AOints.S, left, right)], ranges=[frag_ids,frag_ids])
        self.T = dynamic_array(cache + rule_wrappers + [fragMO_integrals._transform2_rule(    AOints.T, left, right)], ranges=[frag_ids,frag_ids])
        self.U = dynamic_array(cache + rule_wrappers + [fragMO_integrals._transformLast2_rule(AOints.U, left, right)], ranges=[frag_ids,frag_ids,frag_ids])
        self.V = dynamic_array(cache + rule_wrappers + [fragMO_integrals._transform4_rule(    AOints.V, left, right)], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
    @staticmethod
    def _transform2_rule(H, left, right):    # H is overlap or kinetic
        def trans(m1,m2):
            H_ = H[m1,m2]
            if left  is not False:  H_ = left[m1].T @ H_
            if right is not False:  H_ =              H_ @ right[m2]
            return H_
        return trans
    @staticmethod
    def _transformLast2_rule(U, left, right):
        def trans(m1,m2,m3):
            U_ = U[m1,m2,m3]
            if left  is not False:  U_ = left[m2].T @ U_
            if right is not False:  U_ =              U_ @ right[m3]
            return U_
        return trans
    @staticmethod
    def _transform4_rule(V, left, right):
        def trans(m1,m2,m3,m4):
            V_ = V[m1,m2,m3,m4]
            if left  is not False:
                for C in left[m2],left[m1]:    V_ = numpy.tensordot(C, V_, axes=([0],[1]))	# cycle through the former tensor axes (this assumes everything is real)
            if right is not False:
                for C in right[m3],right[m4]:  V_ = numpy.tensordot(V_, C, axes=([2],[0]))	# cycle through the latter tensor axes (this assumes everything is real)
            return V_
        return trans



class spin_orb_integrals(object):
    def __init__(self, spatial_ints, blocking, antisymmetrize=True, cache=False, rule_wrappers=None):
        self.fragments = spatial_ints.fragments
        frag_ids = keys(self.fragments)
        cache = [] if cache else [cached]
        if rule_wrappers is None:  rule_wrappers = []
        self.S = dynamic_array(cache + rule_wrappers + [spin_orb_integrals._spin2_rule(    spatial_ints.S, blocking)], ranges=[frag_ids,frag_ids])
        self.T = dynamic_array(cache + rule_wrappers + [spin_orb_integrals._spin2_rule(    spatial_ints.T, blocking)], ranges=[frag_ids,frag_ids])
        self.U = dynamic_array(cache + rule_wrappers + [spin_orb_integrals._spin2_rule(    spatial_ints.U, blocking)], ranges=[frag_ids,frag_ids,frag_ids])
        if antisymmetrize:  self.V = dynamic_array(cache + rule_wrappers + [spin_orb_integrals._spinAnti4_rule(spatial_ints.V, blocking)], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
        else:               self.V = dynamic_array(cache + rule_wrappers + [spin_orb_integrals._spin4_rule(    spatial_ints.V, blocking)], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
        if hasattr(spatial_ints, "V_half"):
            if antisymmetrize:  self.V_half = dynamic_array(cache + rule_wrappers + [spin_orb_integrals._spinAnti4half_rule(spatial_ints.V_half, blocking)], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
            else:               self.V_half = dynamic_array(cache + rule_wrappers + [spin_orb_integrals._spin4_rule(    spatial_ints.V_half, blocking)], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
    @staticmethod
    def _spin2_rule(H, blocking):    # H is S, T, or U
        def spin2(*m1m2):	# could be m1,m2,m3 in the case of U, but antisymmetrization of block is the same
            return spatial_to_spin.one_electron(H[m1m2], blocking)
        return spin2
    @staticmethod
    def _spin4_rule(V, blocking):
        def spin4(m1,m2,m3,m4):
            return spatial_to_spin.two_electron(V[m1,m2,m3,m4], blocking)
        return spin4
    @staticmethod
    def _spinAnti4_rule(V, blocking):
        def spinAnti4(m1,m2,m3,m4):
            tmp1 = spatial_to_spin.two_electron(V[m1,m2,m3,m4], blocking)
            tmp2 = spatial_to_spin.two_electron(V[m1,m2,m4,m3], blocking)
            return (1/4.) * (tmp1 - numpy.swapaxes(tmp2,2,3))			# premature optimization is the root of all evil ... this will possibly call non-cached result from next layer twice ... must be possible to fix at this layer
        return spinAnti4
    @staticmethod
    def _spinAnti4half_rule(V, blocking):
        def spinAnti4half(m1,m2,m3,m4):
            tmp1 = spatial_to_spin.two_electron(V[m1,m2,m3,m4], blocking)
            tmp2 = spatial_to_spin.two_electron(V[m1,m2,m4,m3], blocking)
            tmp3 = spatial_to_spin.two_electron(V[m2,m1,m3,m4], blocking)
            tmp4 = spatial_to_spin.two_electron(V[m2,m1,m4,m3], blocking)
            return (1/4.) * ((tmp1-numpy.swapaxes(tmp2,2,3)) - numpy.swapaxes(tmp3-numpy.swapaxes(tmp4,2,3),0,1))			# premature optimization is the root of all evil ... this will possibly call non-cached result from next layer twice ... must be possible to fix at this layer
        return spinAnti4half



class bra_transformed(object):
    def __init__(self, transform, frag_integrals, cache=False, rule_wrappers=None, halftrans_rule_wrappers=None):
        self.fragments = frag_integrals.fragments
        frag_ids = keys(self.fragments)
        cache = [] if cache else [cached]
        if rule_wrappers is None:  rule_wrappers = []
        if halftrans_rule_wrappers is None:  halftrans_rule_wrappers = [cached]
        else:                                halftrans_rule_wrappers = [cached] + halftrans_rule_wrappers
        self.S  = dynamic_array(  cache + rule_wrappers + [bra_transformed._transformH_rule( transform, frag_integrals.S)], ranges=[frag_ids,frag_ids])
        self.T  = dynamic_array(  cache + rule_wrappers + [bra_transformed._transformH_rule( transform, frag_integrals.T)], ranges=[frag_ids,frag_ids])
        self.U  = dynamic_array(  cache + rule_wrappers + [bra_transformed._transformU_rule( transform, frag_integrals.U)], ranges=[frag_ids,frag_ids,frag_ids])
        V_h     = dynamic_array(halftrans_rule_wrappers + [bra_transformed._transformV2_rule(transform, frag_integrals.V)], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])    # half transformed ... must cache or else very inefficient?
        self.V  = dynamic_array(  cache + rule_wrappers + [bra_transformed._transformV1_rule(transform, V_h             )], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
        self.V_half = dynamic_array(    cache + rule_wrappers + [bra_transformed._halfV_rule(transform, V_h             )], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
    @staticmethod
    def _transformH_rule(transform, H):		# H can be S or T
        frag_ids = transform.ranges[1]	# The transform might not touch every block of H
        def transform_act(m1,m2):
            H_ = numpy.zeros(H[m1,m2].shape)
            for mm in frag_ids:  H_ += numpy.dot(transform[m1,mm], H[mm,m2])
            return H_
        return transform_act
    @staticmethod
    def _transformU_rule(transform, U):
        frag_ids = transform.ranges[1]	# The transform might not touch every block of U
        def transform_act(m1,m2,m3):
            U_ = numpy.zeros(U[m1,m2,m3].shape)
            for mm in frag_ids:  U_ += numpy.dot(transform[m2,mm], U[m1,mm,m3])
            return U_
        return transform_act
    @staticmethod
    def _transformV2_rule(transform, V):
        frag_ids = transform.ranges[1]	# The transform might not touch every block of V
        def transform_act(m1,m2,m3,m4):
            d1, d2, d3, d4 = V[m1,m2,m3,m4].shape
            V_h = numpy.zeros((d2,d1,d3,d4))
            for mm in frag_ids:  V_h += numpy.tensordot(transform[m2,mm], V[m1,mm,m3,m4], axes=([1],[1]))	# permutes first two indices (only inside the fragment block)
            return V_h
        return transform_act
    @staticmethod
    def _transformV1_rule(transform, V_h):
        frag_ids = transform.ranges[1]	# The transform might not touch every block of V_h
        def transform_act(m1,m2,m3,m4):
            d2, d1, d3, d4 = V_h[m1,m2,m3,m4].shape
            V_ = numpy.zeros((d1,d2,d3,d4))
            for mm in frag_ids:  V_ += numpy.tensordot(transform[m1,mm], V_h[mm,m2,m3,m4], axes=([1],[1]))	# undoes above permutation
            return V_
        return transform_act
    @staticmethod
    def _halfV_rule(transform, V_h):
        def transform_act(m1,m2,m3,m4):
            return numpy.swapaxes(V_h[m1,m2,m3,m4], 0, 1)
        return transform_act



class ket_transformed(object):
    def __init__(self, transform, frag_integrals, cache=False, rule_wrappers=None, halftrans_rule_wrappers=None):
        self.fragments = frag_integrals.fragments
        frag_ids = keys(self.fragments)
        cache = [] if cache else [cached]
        if rule_wrappers is None:  rule_wrappers = []
        if halftrans_rule_wrappers is None:  halftrans_rule_wrappers = [cached]
        else:                                halftrans_rule_wrappers = [cached] + halftrans_rule_wrappers
        self.S  = dynamic_array(  cache + rule_wrappers + [ket_transformed._transformH_rule( transform, frag_integrals.S)], ranges=[frag_ids,frag_ids])
        self.T  = dynamic_array(  cache + rule_wrappers + [ket_transformed._transformH_rule( transform, frag_integrals.T)], ranges=[frag_ids,frag_ids])
        self.U  = dynamic_array(  cache + rule_wrappers + [ket_transformed._transformU_rule( transform, frag_integrals.U)], ranges=[frag_ids,frag_ids,frag_ids])
        V_h     = dynamic_array(halftrans_rule_wrappers + [ket_transformed._transformV2_rule(transform, frag_integrals.V)], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])    # half transformed ... must cache or else very inefficient?
        self.V  = dynamic_array(  cache + rule_wrappers + [ket_transformed._transformV1_rule(transform, V_h             )], ranges=[frag_ids,frag_ids,frag_ids,frag_ids])
    @staticmethod
    def _transformH_rule(transform, H):		# H can be S or T
        frag_ids = transform.ranges[0]	# The transform might not touch every block of H
        def transform_act(m1,m2):
            H_ = numpy.zeros(H[m1,m2].shape)
            for mm in frag_ids:  H_ += numpy.dot(H[m1,mm], transform[mm,m2])
            return H_
        return transform_act
    @staticmethod
    def _transformU_rule(transform, U):
        frag_ids = transform.ranges[0]	# The transform might not touch every block of U
        def transform_act(m1,m2,m3):
            U_ = numpy.zeros(U[m1,m2,m3].shape)
            for mm in frag_ids:  U_ += numpy.dot(U[m1,m2,mm], transform[mm,m3])
            return U_
        return transform_act
    @staticmethod
    def _transformV2_rule(transform, V):
        frag_ids = transform.ranges[0]	# The transform might not touch every block of V
        def transform_act(m1,m2,m3,m4):
            d1, d2, d3, d4 = V[m1,m2,m3,m4].shape
            V_h = numpy.zeros((d2,d1,d3,d4))
            for mm in frag_ids:  V_h += numpy.tensordot(V[m1,m2,mm,m4], transform[mm,m3], axes=([2],[0]))	# permutes first two indices (only inside the fragment block)
            return V_h
        return transform_act
    @staticmethod
    def _transformV1_rule(transform, V_h):
        frag_ids = transform.ranges[0]	# The transform might not touch every block of V_h
        def transform_act(m1,m2,m3,m4):
            d2, d1, d3, d4 = V_h[m1,m2,m3,m4].shape
            V_ = numpy.zeros((d1,d2,d3,d4))
            for mm in frag_ids:  V_ += numpy.tensordot(V_h[m1,m2,m3,mm], transform[mm,m4], axes=([2],[0]))	# undoes above permutation
            return V_
        return transform_act






def _add_rule(A, B, Aranges, Branges):
    def add(m1,m2):
        if (m1 in Aranges[0]) and (m2 in Aranges[1]) and (m1 in Branges[0]) and (m2 in Branges[1]):
            return A[m1,m2] + B[m1,m2]
        elif (m1 in Aranges[0]) and (m2 in Aranges[1]):
            return A[m1,m2]
        elif (m1 in Branges[0]) and (m2 in Branges[1]):
            return B[m1,m2]
        else:
            return None		# should I get the fragments structure so I can allocate a block of zeros of the correct size?
    return add

def _subtract_rule(A, B, Aranges, Branges):
    def add(m1,m2):
        if (m1 in Aranges[0]) and (m2 in Aranges[1]) and (m1 in Branges[0]) and (m2 in Branges[1]):
            return A[m1,m2] - B[m1,m2]
        elif (m1 in Aranges[0]) and (m2 in Aranges[1]):
            return  A[m1,m2]
        elif (m1 in Branges[0]) and (m2 in Branges[1]):
            return -B[m1,m2]
        else:
            return None		# should I get the fragments structure so I can allocate a block of zeros of the correct size?
    return add

def _mat_mul_rule(A, B):
    try:     frag_ids_tmp = list(A.ranges[1])   # At least one of these ...
    except:  frag_ids_tmp = list(B.ranges[0])   # ... has to be convertable to a list.
    frag_ids = []
    for frag_id in frag_ids_tmp:
        if (frag_id in A.ranges[1]) and (frag_id in B.ranges[0]):  frag_ids += [frag_id]	# assume that undefined columns-of-A or rows-of-B are zero
    def multiply(m1,m2):
        d1 = A[m1,frag_ids[0]].shape[0]    # if not same for all frag_ids ...
        d2 = B[frag_ids[0],m2].shape[1]    # ... lines below will fail (as is the only sensible option).
        C = numpy.zeros((d1,d2))
        for mm in frag_ids:  C += numpy.dot(A[m1,mm], B[mm,m2])
        return C
    return multiply


class range_union(object):
    def __init__(self, rangeA, rangeB):
        self.rangeA = rangeA
        self.rangeB = rangeB
    def __contains__(self,key):
        if (key in self.rangeA) or (key in self.rangeB):  return True
        else:  return False

def add(A,B, rules=None):
    if rules is None:  rules = []
    ranges = range_union(A.ranges[0],B.ranges[0]), range_union(A.ranges[1],B.ranges[1])		# assume that undefined blocks of A or B are zero
    return dynamic_array(rules + [_add_rule(A,B,A.ranges,B.ranges)], ranges=ranges)

def subtract(A,B, rules=None):
    if rules is None:  rules = []
    ranges = range_union(A.ranges[0],B.ranges[0]), range_union(A.ranges[1],B.ranges[1])		# assume that undefined blocks of A or B are zero
    return dynamic_array(rules + [_subtract_rule(A,B,A.ranges,B.ranges)], ranges=ranges)

def mat_mul(A,B, rules=None):
    if rules is None:  rules = []
    return dynamic_array(rules + [_mat_mul_rule(A,B)], ranges=[A.ranges[0],B.ranges[1]])



def _zeros1_rule_raw(fragments, spin_orbs):
    spin_orbs = 2 if spin_orbs else 1
    def allocate(mm):
        dim = spin_orbs * fragments[mm].basis.n_spatial_orb
        return numpy.zeros(dim)
    return allocate

def zeros1(fragments, rules=None, ranges=None, spin_orbs=False):
    frag_ids = keys(fragments)
    if rules is None:  rules = []
    if ranges is None:  ranges = [frag_ids]
    return dynamic_array(rules + [cached, _zeros1_rule_raw(fragments, spin_orbs)], ranges=ranges)	# cached does not mean it cannot be modified, it just means it is not reallocated (reset when a new block is requested)

def _zeros2_rule_raw(fragments, spin_orbs):
    spin_orbs = 2 if spin_orbs else 1
    def allocate(m1,m2):
        d1 = spin_orbs * fragments[m1].basis.n_spatial_orb
        d2 = spin_orbs * fragments[m2].basis.n_spatial_orb
        return numpy.zeros((d1,d2))
    return allocate

def zeros2(fragments, rules=None, ranges=None, spin_orbs=False):
    frag_ids = keys(fragments)
    if rules is None:  rules = []
    if ranges is None:  ranges = [frag_ids,frag_ids]
    return dynamic_array(rules + [cached, _zeros2_rule_raw(fragments, spin_orbs)], ranges=ranges)	# cached does not mean it cannot be modified, it just means it is not reallocated (reset when a new block is requested)



def _Id_rule_raw(fragments, spin_orbs):
    spin_orbs = 2 if spin_orbs else 1
    def allocate(m1,m2):
        d1 = spin_orbs * fragments[m1].basis.n_spatial_orb
        d2 = spin_orbs * fragments[m2].basis.n_spatial_orb
        if m1==m2:  return numpy.identity(d1)
        else:       return numpy.zeros((d1,d2))
    return allocate

def Id(fragments, rules=None, ranges=None, spin_orbs=False):
    frag_ids = keys(fragments)
    if rules is None:  rules = []
    if ranges is None:  ranges = [frag_ids,frag_ids]
    return dynamic_array(rules + [cached, _Id_rule_raw(fragments, spin_orbs)], ranges=ranges)	# cached does not mean it cannot be modified, it just means it is not reallocated (reset when a new block is requested)


def column_rule(matrix, m2, i2):
    def get_block(m1):
        return matrix[m1,m2][:,i2]
    return get_block

class mat_as_columns(object):
    def __init__(self, matrix):
        self.matrix = matrix
    def __call__(self, m2, i2):
        return dynamic_array(column_rule(self.matrix, m2, i2), ranges=[self.matrix.ranges[0]])


def row_rule(matrix, m1, i1):
    def get_block(m2):
        return matrix[m1,m2][i1,:]
    return get_block

class mat_as_rows(object):
    def __init__(self, matrix):
        self.matrix = matrix
    def __call__(self, m1, i1):
        return dynamic_array(row_rule(self.matrix, m1, i1), ranges=[self.matrix.ranges[1]])



def as_raw_mat(matrix, fragments, ranges=None, spin_orbs=False):
    if ranges is None:  ranges = matrix.ranges
    return unblock_2(matrix, fragments, ranges, spin_orbs)

def as_frag_blocked_mat(matrix, fragments, ranges=None, spin_orbs=False):
    frag_ids = keys(fragments)
    if ranges is None:  ranges = [frag_ids,frag_ids]
    M = block_2(matrix, fragments, ranges, spin_orbs)
    return mask(M, ranges=ranges)


def copy_rule(to_copy):
    def copy_block(mm):
        return numpy.array(to_copy[mm])
    return copy_block


class space_traits(object):
    def __init__(self, field):
        self.field = field
    def dot(self,v,w):
        result = self.field(0)
        vlen = len(v.ranges[0])
        wlen = len(w.ranges[0])
        if vlen<wlen:  frag_ids = v.ranges[0]   # if the longer one is not a superset
        else:          frag_ids = w.ranges[0]   # the line below will fail
        for mm in frag_ids:  result += numpy.dot(v[mm], w[mm])
        return result
    def add_to(self,v,w,c=1):
        frag_ids = w.ranges[0]   # if v does not have slots to resolve all of w, the below will fail
        for mm in frag_ids:
            tmp = v[mm]		# not sure why
            tmp += c*w[mm]	# v[mm] += c*w[mm] sometimes fails ... because dynamic_array do not have __setelem__ defined(?) ...  and they should not.
    def scale(self,n,v):
        frag_ids = v.ranges[0]
        for mm in frag_ids:
            tmp = v[mm]		# not sure why
            tmp *= n		# v[mm] *= n sometimes fails
    def copy(self,v):
        frag_ids = v.ranges[0]
        u = dynamic_array([cached, copy_rule(v)], ranges=[frag_ids])
        for mm in frag_ids:  u[mm]	# forces the copy to take place right now
        return u
    def check_member(self,v):  pass
