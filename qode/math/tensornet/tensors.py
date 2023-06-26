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

from copy      import copy, deepcopy
from textwrap  import indent
from base      import tensor_base, evaluate, increment, raw, scalar, backend, scalar_value, shape
from heuristic import heuristic    # how to order contraction executions in a network



# With a syntax like
#    tensor = primitive_tensor_wrapper(tensornet.numpy_backend)
#    A = tensor(numpy.array(...))
# the user can easily wrap all the raw tensors they like in their chosen backend.
def primitive_tensor_wrapper(backend):
    def wrapper(raw_tensor):
        return primitive_tensor(raw_tensor, backend)
    return wrapper

# in an informational sense, this is the most basic tensor, since all it does is store a list
# of tensor_networks and primitive_tensors.
class tensor_sum(tensor_base):
    def __init__(self, tensor_terms=None):
        if tensor_terms==None:  tensor_terms = []    # for instantiation of empty sum as accumulator
        self.tensornet_backend = None
        self.tensornet_shape   = None
        self.tensor_terms      = []
        for term in tensor_terms:
            if self.tensornet_backend is None:
                try:
                    self.tensornet_backend = backend(term)
                    self.tensornet_shape   = shape(term)
                except:
                    raise TypeError("only tensornet tensors can be summed")
            if backend(term) is not self.tensornet_backend or shape(term)!=self.tensornet_shape:
                raise ValueError("only tensornet tensors with the same backend and shape can be summed")
            try:
                tensor_subterms = term.tensor_terms
            except AttributeError:
                new_terms = [copy(term)]                                        # we want copies, ...
            else:
                new_terms = [copy(sub_term) for sub_term in tensor_subterms]    # ... in case we use *=
            self.tensor_terms += new_terms
    def tensornet_evaluate(self):
        result = raw(evaluate(self.tensor_terms[0]))
        for term in self.tensor_terms[1:]:
            increment(result, term)    # move actual math out of here and let child class decide how
        return primitive_tensor(result, backend(self))
    def __imul__(self, x):    # enables __mul__, __rmul__, __neg__, and therefore also __sub__
        for term in self.tensor_terms:
            term *= x    # will change only scalar prefactors, not original data provided by user
        return self
    def __iadd__(self, other):
        try:
            other_backend = backend(other)
        except:
            raise TypeError("only tensornet tensors can be added to a tensornet tensor_sum")
        if other_backend is not backend(self):
            raise ValueError("only tensornet tensors with the same backend can be added")
        try:
            other_tensor_terms = other.tensor_terms
        except AttributeError:
            new_terms = [copy(other)]                                              # we want copies, ...
        else:
            new_terms = [copy(other_term) for other_term in other_tensor_terms]    # ... in case we use *=
        self.tensor_terms += new_terms
        return self
    def __isub__(self, other):    # enabled by __iadd__ and __imul__ (which enabled __neg__)
        self += -other
        return self



# The tensornet type for the primitive tensors that the user sees and uses and builds networks from.
# It should mostly act like the raw tensor type of the chosen backend (ie, forwarding calls), but with
# some extra tensornet-specific functionality that is uniform across backends.  In order to accomplish
# this, the tensor importantly knows its backend, via a provided module (implemented by the user if not
# already provided for that backend type).  
class primitive_tensor(tensor_base):
    def __init__(self, raw_tensor, backend, _scalar=1):
        # long member names to stay out of way of __getattr__ (and never seen by user)
        self.tensornet_raw_tensor = raw_tensor             # only ever accessed by wrapper below
        self.tensornet_backend    = backend
        self.tensornet_scalar     = _scalar                # here so that we can define *= without changing original data
        self.tensornet_shape      = backend.shape(raw_tensor)
    def tensornet_scalar_value(self):                      # only ever accessed by wrapper below
        return scalar(self) * backend(self).scalar_value(raw(self))
    def tensornet_increment(self, result):
        result += scalar(self) * raw(self)       # do not make a copy just to use as an increment, but ...
        return
    def tensornet_evaluate(self):
        return primitive_tensor(deepcopy(raw(self)), backend(self), _scalar=scalar(self))    # ... this might get modified and so should be independent
    #def __getattr__(self, name):
    #    return getattr(raw(self), name)    # forward all other calls for user convenience
    # calls to __XXX__ functions are not passed through by __getattr__ ... for good reason (would return wrong types)
    def __str__(self):
        return "tensornet.primitive_tensor(\n{}\n)".format(indent(str(raw(self)), "    "))
    # __iadd__ and __isub__ would be confusing since __add__ and __sub__ make a tensor_sum (*)
    # but the increment operators should be of input type.  (* this is for the best because it
    # is more flexible; the user can choose to do a hard data-level add outside of tensornet)
    def __imul__(self, x):    # enables __mul__, __rmul__, __neg__, and therefore also __sub__
        self.tensornet_scalar *= x
        return self



class tensor_network(tensor_base):
    def __init__(self, scalar, contractions, free_indices, backend):
        self.tensornet_scalar = scalar
        self.contractions = contractions
        self.free_indices = free_indices
        self.tensornet_backend = backend
        new_shape = []
        for free_index in free_indices:
            tens0, pos0 = free_index[0]      # always instantiated by function that checks congruence
            new_shape += [shape(tens0)[pos0]]
        self.tensornet_shape = tuple(new_shape)
    def tensornet_increment(self, result):
        result += raw(evaluate(self))
        return
    def tensornet_evaluate(self):
        # It is assumed that all of the tensors in the network are represented by distinct
        # objects, even if they point to the same underlying data.  This is enforced by
        # the contract function, which is the only way for a user to make a tensor_network.
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
        shapes = {tens:shape(by_id[tens]) for tens in by_id}
        #
        do_scalar_mult, do_reduction, target = heuristic(scalar(self), contraction_groups, index_reduct_groups, shapes)
        #
        if do_scalar_mult or do_reduction:
            if do_scalar_mult:
                new_scalar = 1
                other_contractions = []
                for group,contraction_sublist in contraction_groups.items():
                    other_contractions += contraction_sublist
                mapping = {target:list(range(len(shape(by_id[target]))))}
            else:
                new_scalar = scalar(self)
                other_contractions = []
                for group,contraction_sublist in contraction_groups.items():
                    if group!=target:
                        other_contractions += contraction_sublist
                mapping = {tens:[None]*len(shape(by_id[tens])) for tens in target}
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
            args = [(raw(by_id[tens]), *indices) for tens,indices in mapping.items()]
            if do_scalar_mult:  args += [scalar(self)]
            #
            new_tens = primitive_tensor(backend(self).contract(*args), backend(self))
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
            if len(shape(new_tens))==0:
                new_scalar *= scalar_value(new_tens)    # in no way not a scalar (unlike a 1x1x1x... tensor).  note that new_tens itself is now forgotten
            return evaluate(tensor_network(new_scalar, new_contractions, new_free_indices, backend(self)))    # recur
        else:
            mapping = {}
            for i,free_index in enumerate(free_indices):
                tens, pos = free_index[0]    # guaranteed to be only one entry per index now
                if tens not in mapping:
                    mapping[tens] = [None]*len(shape(by_id[tens]))
                mapping[tens][pos] = i
            args = [(raw(by_id[tens]), *indices) for tens,indices in mapping.items()] + [scalar(self)]
            return primitive_tensor(backend(self).contract(*args), backend(self))                             # bottom out
    def __imul__(self, x):
        self.tensornet_scalar *= x
        return self
