#    (C) Copyright 2023 Anthony D. Dutoi
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

from .base      import evaluate, raw, scalar_value, resolve_ellipsis, timings_start, timings_record
from .tensors   import summable_tensor, tensor_sum, primitive_tensor
from .heuristic import heuristic    # how to order contraction executions in a network

_backend_contract_path = False    # if True, let backend handle finding the optimal contraction path upon evaluate() call
_warned = False                   # have we warned the user yet against asking for individual tensor elements?

def backend_contract_path(TrueFalse):
    global _backend_contract_path
    _backend_contract_path = TrueFalse



# Barebones theory (written much later after a forensic debug battle).  A tensor_network object contains
# two fundamental pieces of information (and other incidental info).  One is a list of contractions.
# Each contraction itself is a list of two-tuples; each two-tuple identifies a tensor and the index
# of that tensor involved in the contraction.  So if a contraction contains two two-tuples, then 
# two indices are contracted with one another, but there might be more for unusual contractions.
# The ordering of the list of contractions is irrelevant, as is the ordering of the list of two-tuples
# that defines a given contraction. The other piece of information is a list of free indices of the result,
# which correspond to uncontracted indices of the tensors in the network.  This list is ordered, and each
# free index is itself a list of two-tuples with the same tensor-index structure as the two-tuples that
# define a contraction.  In the most usual case, each such list corresponding to a free index will be of
# length one, but if there is more than one tensor index that corresponds to a single free index it is because
# those indices are set equal to each other and reduced to a single free index.

class tensor_network(summable_tensor):
    def __init__(self, scalar, contractions, free_indices, backend, contract):
        shape = []
        for free_index in free_indices:
            tens0, pos0 = free_index[0]      # always instantiated by function that checks congruence
            shape += [tens0.shape[pos0]]
        summable_tensor.__init__(self, tuple(shape), backend, contract)
        self._scalar = scalar
        self._contractions = contractions
        self._free_indices = free_indices
        #
        by_id = {}    # populated by calls to _hashable() below
        def _hashable(prim_lists):
            nonlocal by_id
            new_prim_lists = []
            for prim_list in prim_lists:
                new_prim_list = []
                for tens,pos in prim_list:
                    tens_id = id(tens)
                    by_id[tens_id] = tens
                    new_prim_list += [(tens_id, pos)]
                new_prim_lists += [new_prim_list]
            return new_prim_lists
        contractions = _hashable(self._contractions)
        free_indices = _hashable(self._free_indices)
        self._hashable = by_id, contractions, free_indices    # hashable (redundant) catalogs for internal use only
        # Nested tuples of hashable values like that below can automatically be sorted lexicographically, and they can be tested 
        # for equality.  Although false *negatives* can happen, if they are equal, the result is guaranteed to be the same.  This
        # can be useful in sums where permutationally (anti)symmetric tensors are built from asymmetric primitives and then contracted.
        # Redundant terms can be recognized in n*log(n) time by sorting and grouping.
        self._result_hash = ( tuple(sorted(tuple(sorted((id(tens._raw_tensor), pos) for tens,pos in prim_list)) for prim_list in self._contractions)),
                                     tuple(tuple(sorted((id(tens._raw_tensor), pos) for tens,pos in prim_list)) for prim_list in self._free_indices) )
    def _increment(self, result):
        result += raw(self)
        return
    def _evaluate(self):
        if self._scalar==0:    # usually a bad test, but in this case we really mean it.  If it is not exactly zero, there is something to do, and zero can happen (ie, a = 0 * b)
            return primitive_tensor(self._backend.zeros(self.shape), self._backend, self._contract)
        timings_start()
        # It is assumed that all of the tensors in the network are represented by distinct
        # objects, even if they point to the same underlying data.  This is enforced by
        # the contract function, which is the only way for a user to make a tensor_network.
        # Also, the individual primitives have all been forced to have unit scalar by adjusting the overall scalar.
        by_id, contractions, free_indices = self._hashable
        #
        if _backend_contract_path:    # All the tensors in one big group
            do_scalar_mult = False if self._scalar==1 else True
            do_reduction   = True
            tens_group = tuple(sorted({tens for contraction in contractions for tens,_ in contraction}
                                    | {tens for free_index  in free_indices for tens,_ in free_index }))
            contraction_groups  = {tens_group: contractions}
            index_reduct_groups = {tens_group: free_indices}
            target = tens_group
        else:                         # Use tensornet contraction path
            def _group_by_tensors(prim_lists, allow_singles=False):
                prim_list_groups = {}
                for prim_list in prim_lists:
                    if len(prim_list)>1 or allow_singles:
                        tens_group = tuple(sorted({tens for tens,_ in prim_list}))    # tuples of sorted ids can be used as dict keys
                        if tens_group not in prim_list_groups:
                            prim_list_groups[tens_group] = []
                        prim_list_groups[tens_group] += [prim_list]
                return prim_list_groups
            contraction_groups  = _group_by_tensors(contractions, allow_singles=True)
            index_reduct_groups = _group_by_tensors(free_indices)
            shapes = {tens:by_id[tens].shape for tens in by_id}
            timings_record("tensor_network._evaluate")
            timings_start()
            do_scalar_mult, do_reduction, target = heuristic(self._scalar, contraction_groups, index_reduct_groups, shapes)
            timings_record("heuristic")
        #
        if do_scalar_mult or do_reduction:
            timings_start()
            if do_reduction:    # as if do_scalar_mult is False, which it will be if using tensornet contraction path
                # print("do_reduction")
                scalar = self._scalar
                other_contractions = []
                for group,contraction_sublist in contraction_groups.items():
                    if group!=target:
                        other_contractions += contraction_sublist
                mapping = {tens:[None]*len(by_id[tens].shape) for tens in target}
                def _letter(idx):
                    if idx<26:  return "abcdefghijklmnopqrstuvwxyz"[idx]
                    else:       return str(idx)    # obfuscated if ever printed for >26 indices (?!)
                if target in contraction_groups:
                    for i,contraction in enumerate(contraction_groups[target]):
                        for tens,pos in contraction:
                            mapping[tens][pos] = _letter(i)
                i = 0
                if target in index_reduct_groups:
                    for free_index in index_reduct_groups[target]:
                        for tens,pos in free_index:
                            mapping[tens][pos] = i
                        i += 1
                for tens in mapping:
                    for j in range(len(mapping[tens])):
                        if mapping[tens][j] is None:
                            mapping[tens][j] = i
                            i += 1
            else:
                # print("do_scalar_mult")
                other_contractions = []
                for group,contraction_sublist in contraction_groups.items():
                    other_contractions += contraction_sublist
                mapping = {target:list(range(len(by_id[target].shape)))}
            args = [(by_id[tens]._raw_tensor, *indices) for tens,indices in mapping.items()]
            if do_scalar_mult:
                args += [self._scalar]
                scalar = 1
            timings_record("tensor_network._evaluate")
            #
            timings_start()
            new_tens = primitive_tensor(self._backend.contract(*args), self._backend, self._contract)
            timings_record("backend.contract")
            if _backend_contract_path:
                return new_tens    # bottom out immediately ... there must be a more elegant way of switching between all these options!
            #
            timings_start()
            def _map_indices(prim_lists):
                new_prim_lists = []
                for prim_list in prim_lists:
                    new_prim_list = []
                    done = []
                    for tens,pos in prim_list:
                        if (tens,pos) not in done:
                            done += [(tens,pos)]
                            if tens in mapping:
                                new_prim_list += [(new_tens, mapping[tens][pos])]
                            else:
                                new_prim_list += [(by_id[tens], pos)]
                    new_prim_lists += [new_prim_list]
                return new_prim_lists
            new_contractions = _map_indices(other_contractions)    # guaranteed safe, ...
            new_free_indices = _map_indices(free_indices)          # ... even if new_tens is 0-dim
            if len(new_tens.shape)==0:
                scalar *= scalar_value(new_tens)    # in no way not a scalar (unlike a 1x1x1x... tensor).  note that new_tens itself is now forgotten
            timings_record("tensor_network._evaluate")
            return evaluate(tensor_network(scalar, new_contractions, new_free_indices, self._backend, self._contract))    # recur
        else:    # must be a single tensor (?) or an outer product
            if len(free_indices)>0:
                timings_start()
                mapping = {}
                for i,free_index in enumerate(free_indices):
                    tens, pos = free_index[0]    # guaranteed to be only one entry per index now
                    if tens not in mapping:
                        mapping[tens] = [None]*len(by_id[tens].shape)
                    mapping[tens][pos] = i
                if self._scalar!=1:                                     # eventually deprecate.
                    raise RuntimeError("scalar should be 1, right?")    # pretty sure this has to be true
                args = [(by_id[tens]._raw_tensor, *indices) for tens,indices in mapping.items()]    # self._scalar is 1 by now?
                timings_record("tensor_network._evaluate")
                timings_start()
                Z = primitive_tensor(self._backend.contract(*args), self._backend, self._contract)                     # bottom out (might give a 0-dim tensor; this is intended)
                timings_record("backend.contract")
                return Z
            else:    # there is nothing left but the scalar
                return primitive_tensor(self._backend.scalar_tensor(self._scalar), self._backend, self._contract)
    def __getitem__(self, indices):
        indices = resolve_ellipsis(indices, self.shape)
        full = slice(None)    # the slice produced by [:] with no limits
        scalar = self._scalar
        by_id, contractions, free_indices = self._hashable
        slice_indices = {}
        for index,free_index in zip(indices,free_indices):
            if index!=full:
                for tens,pos in free_index:
                    if tens not in slice_indices:
                        slice_indices[tens] = [full]*len(by_id[tens].shape)
                    slice_indices[tens][pos] = index
        mappings = {}
        for tens,tens_slice in slice_indices.items():
            new_tens = by_id[tens][tuple(tens_slice)]
            j, mapping= 0, {}
            for i,index in enumerate(tens_slice):
                if isinstance(index,slice):
                    mapping[i] = j
                    j += 1
            if j==0:                  # this tensor has no indices left in contractions or free_indices
                scalar *= new_tens    # fully indexed should be a scalar.  will be forgotten when removed from contractions and free_indices
            else:
                mappings[tens] = (new_tens, mapping)
        def _map_indices(prim_list):
            new_prim_list = []
            for tens,pos in prim_list:
                if tens in mappings:
                    new_tens, mapping = mappings[tens]
                    new_prim_list += [(new_tens, mapping[pos])]    # should logically never happen that a missing entry is requested 
                else:
                    new_prim_list += [(by_id[tens],pos)]
            return new_prim_list
        new_free_indices, new_contractions = [], []
        for index,free_index in zip(indices,free_indices):
            if isinstance(index,slice):
                new_free_indices += [_map_indices(free_index)]
        for contraction in contractions:
            new_contractions += [_map_indices(contraction)]
        new = tensor_network(scalar, new_contractions, new_free_indices, self._backend, self._contract)
        if len(new_free_indices)==0:
            global _warned
            if not _warned:
                print("Accessing individual elements of a tensor network is really inefficient.  Consider alternatives.")    # How user can suppress this altogether?
                _warned = True
            return scalar_value(new)
        else:
            return new
    def __imul__(self, x):
        self._scalar *= x
        return self
