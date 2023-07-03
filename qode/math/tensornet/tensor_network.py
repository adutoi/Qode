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

from .base      import tensor_base, evaluate, extract, scalar_value
from .tensors   import summable_tensor, tensor_sum, primitive_tensor
from .heuristic import heuristic    # how to order contraction executions in a network

warned = False    # have we warned the user yet against asking for individual tensor elements?



class tensor_network(summable_tensor):
    def __init__(self, scalar, contractions, free_indices, backend):
        self._backend = backend
        shape = []
        for free_index in free_indices:
            tens0, pos0 = free_index[0]      # always instantiated by function that checks congruence
            shape += [tens0.shape[pos0]]
        self.shape = tuple(shape)
        self._scalar = scalar
        self._contractions = contractions
        self._free_indices = free_indices
        #
        by_id = {}
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
    def _increment(self, result):
        result += extract(self)
        return
    def _evaluate(self):
        # It is assumed that all of the tensors in the network are represented by distinct
        # objects, even if they point to the same underlying data.  This is enforced by
        # the contract function, which is the only way for a user to make a tensor_network.
        # also, the individual primitives have all been forced to have unit scalar by adjusting the overall scalar
        by_id, contractions, free_indices = self._hashable
        def _group_by_tensors(prim_lists, allow_singles=False):
            prim_list_groups = {}
            for prim_list in prim_lists:
                if len(prim_list)>1 or allow_singles:
                    tens_group = tuple(sorted({tens for tens,_ in prim_list}))    # tuples of sorted ids can be used as dict keys
                    if tens_group not in prim_list_groups:
                        prim_list_groups[tens_group] = []
                    prim_list_groups[tens_group] += [prim_list]
            return prim_list_groups
        # print("contractions", contractions)
        # print("free_indices", free_indices)
        contraction_groups  = _group_by_tensors(contractions, allow_singles=True)
        index_reduct_groups = _group_by_tensors(free_indices)
        shapes = {tens:by_id[tens].shape for tens in by_id}
        #
        do_scalar_mult, do_reduction, target = heuristic(self._scalar, contraction_groups, index_reduct_groups, shapes)
        #
        if do_scalar_mult or do_reduction:
            if do_scalar_mult:
                # print("do_scalar_mult")
                scalar = 1
                other_contractions = []
                for group,contraction_sublist in contraction_groups.items():
                    other_contractions += contraction_sublist
                mapping = {target:list(range(len(by_id[target].shape)))}
            else:
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
            args = [(by_id[tens]._raw_tensor, *indices) for tens,indices in mapping.items()]
            if do_scalar_mult:  args += [self._scalar]
            #
            new_tens = primitive_tensor(self._backend.contract(*args), self._backend)
            #
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
            return evaluate(tensor_network(scalar, new_contractions, new_free_indices, self._backend))    # recur
        else:
            if len(free_indices)>0:
                mapping = {}
                for i,free_index in enumerate(free_indices):
                    tens, pos = free_index[0]    # guaranteed to be only one entry per index now
                    if tens not in mapping:
                        mapping[tens] = [None]*len(by_id[tens].shape)
                    mapping[tens][pos] = i
                args = [(by_id[tens]._raw_tensor, *indices) for tens,indices in mapping.items()] + [self._scalar]
                return primitive_tensor(self._backend.contract(*args), self._backend)                             # bottom out (might give a 0-dim tensor; this is intended)
            else:    # there is nothing left but the scalar
                return primitive_tensor(self._backend.scalar_tensor(self._scalar), self._backend)
    def __getitem__(self, indices):
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
        new = tensor_network(scalar, new_contractions, new_free_indices, self._backend)
        if len(new_free_indices)==0:
            global warned
            if not warned:
                print("Accessing individual elements of a tensor network is really inefficient.  Consider alternatives.")    # How user can suppress this altogether?
                warned = True
            return scalar_value(new)
        else:
            return new
    def __imul__(self, x):
        self._scalar *= x
        return self

