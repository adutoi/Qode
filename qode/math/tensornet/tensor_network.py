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

from primitive_tensor import primitive_tensor, thin_wrapper
from heuristic        import heuristic



# adds __mul__, __rmul__, and __neg__ for child classes that implement           __imul__ and tensornet_copy
# adds __sub__                        for child classes that implement __add__,  __imul__ and tensornet_copy
# adds __isub__                       for child classes that implement __iadd__, __imul__ and tensornet_copy
class tensor_network_base(object):
    def __mul__(self, x):
        new = self.tensornet_copy()    # long name so that intended use case only works with child classes
        new *= x
        return new
    def __rmul__(self, x):
        return self * x
    def __neg__(self):
        return self * -1
    def __sub__(self, other):
        return self + (-other)
    def __isub__(self, other):
        self += -other
        return



class tensor_network_sum(tensor_network_base):
    def __init__(self, tensor_terms):
        self.tensor_terms = []
        for term in tensor_terms:
            try:
                sub_tensor_terms = term.tensor_terms
            except AttributeError:
                new_terms = [term.tensornet_copy()]
            else:
                new_terms = [sub_term.tensornet_copy() for sub_term in sub_tensor_terms]
            self.tensor_terms += new_terms
    def tensornet_copy(self):
        return tensor_network_sum([self])
    def tensornet_evaluate(self):
        result = self.tensor_terms[0].tensornet_evaluate()
        for term in self.tensor_terms[1:]:
            result += term.tensornet_evaluate()
        return result
    def __iadd__(self, other):
        try:
            other_tensor_terms = other.tensor_terms
        except AttributeError:
            new_terms = [other.tensornet_copy()]
        else:
            new_terms = [other_term.tensornet_copy() for other_term in other_tensor_terms]
        self.tensor_terms += new_terms
        return
    def __add__(self, other):
        new = self.tensornet_copy()
        new += other
        return new_sum
    def __imul__(self, x):
        for term in self.tensor_terms:
            term *= x
        return



class tensor_network(tensor_network_base):
    def __init__(self, scalar, contractions, free_indices, backend):
        self.scalar       = scalar                #
        self.contractions = list(contractions)    # only accessed when we are sure it is a tensor_network
        self.free_indices = list(free_indices)    #
        self.tensornet_backend = backend
        shape = []
        for free_index in free_indices:
            tens0, pos0 = free_index[0]      # always instantiated by function that checks congruence
            shape += [tens0.shape[pos0]]
        self.tensornet_shape = tuple(shape)  # shape needs to act like unwrapped ...
        self.shape = self.tensornet_shape    # ... or wrapped primitive_tensor
    def tensornet_copy(self):
        return tensor_network(self.scalar, self.contractions, self.free_indices, self.tensornet_backend)
    def tensornet_evaluate(self):
        by_id = {}    # to resolve hashable references to the tensors
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
        contractions = _hashable(self.contractions)
        free_indices = _hashable(self.free_indices)
        #
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
        #
        do_scalar_mult, do_reduction, target = heuristic(self.scalar, contraction_groups, index_reduct_groups, by_id)
        #
        if do_scalar_mult or do_reduction:
            if do_scalar_mult:
                new_scalar = 1
                other_contractions = []
                for group,contraction_sublist in contraction_groups.items():
                    other_contractions += contraction_sublist
                mapping = {target:list(range(len(by_id[target].shape)))}
            else:
                new_scalar = self.scalar
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
            args = [(by_id[tens].raw, *indices) for tens,indices in mapping.items()]
            if do_scalar_mult:  args += [self.scalar]
            #
            new_tens = thin_wrapper(primitive_tensor(self.tensornet_backend.contract(*args), self.tensornet_backend))
            #
            def _map_indices(prim_lists):
                new_prim_lists = []
                for prim_list in prim_lists:
                    new_prim_list = []
                    done = []
                    for tens,pos in prim_list:
                        if tens not in done:
                            done += [tens]
                            if tens in mapping:
                                new_prim_list += [(new_tens, mapping[tens][pos])]
                            else:
                                new_prim_list += [(by_id[tens], pos)]
                    new_prim_lists += [new_prim_list]
                return new_prim_lists
            new_contractions = _map_indices(other_contractions)    # guaranteed safe, ...
            new_free_indices = _map_indices(free_indices)          # ... even if new_tens is 0-dim
            if len(new_tens.shape)==0:
                new_scalar *= new_tens.scalar_value()    # in no way not a scalar (unlike a 1x1x1x... tensor).  note that new_tens itself is now forgotten
            return tensor_network(new_scalar, new_contractions, new_free_indices, self.tensornet_backend).tensornet_evaluate()    # recur
        else:
            mapping = {}
            for i,free_index in enumerate(free_indices):
                tens, pos = free_index[0]    # guaranteed to be only one entry per index now
                if tens not in mapping:
                    mapping[tens] = [None]*len(by_id[tens].shape)
                mapping[tens][pos] = i
            args = [(by_id[tens].raw, *indices) for tens,indices in mapping.items()] + [self.scalar]
            return primitive_tensor(self.tensornet_backend.contract(*args), self.tensornet_backend)                               # bottom out
    def __add__(self, other):
        return tensor_network_sum([self, other])
    def __imul__(self, x):
        self.scalar *= x
        return
