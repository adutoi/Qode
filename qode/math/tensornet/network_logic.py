#    (C) Copyright 2023, 2025 Anthony D. Dutoi
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

# This is where the wild things are.  All the actual implementation code is buried in here,
# separated from the interface C++ style.  Is still kind of a heap of chaos, but at least
# the highlevel stuff is there to guide what one expects here.  Clean up someday . . .

import copy
from ...util import struct
from .base import evaluate, scalar_value, timings_start, timings_record, ContractionError
from .heuristic import heuristic    # how to order contraction executions in a network

_warned = False                   # have we warned the user yet against asking for individual tensor elements?



class logic(object):    # only used to be lazy about not changing indentation
    @staticmethod
    def init(backend, scalar, contractions, free_indices):
        values = struct()
        shape = []
        for free_index in free_indices:
            tens0, pos0 = free_index[0]      # always instantiated by function that checks congruence
            shape += [tens0.shape[pos0]]
        values.backend = backend
        values.shape = tuple(shape)
        values.scalar = scalar
        values.contractions = contractions
        values.free_indices = free_indices
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
        contractions = _hashable(values.contractions)
        free_indices = _hashable(values.free_indices)
        values.hashable = by_id, contractions, free_indices    # hashable (redundant) catalogs for internal use only
        # Nested tuples of hashable values like that below can automatically be sorted lexicographically, and they can be tested 
        # for equality.  Although false *negatives* can happen, if they are equal, the result is guaranteed to be the same.  This
        # can be useful in sums where permutationally (anti)symmetric tensors are built from asymmetric primitives and then contracted.
        # Redundant terms can be recognized in n*log(n) time by sorting and grouping.
        values.result_hash = ( tuple(sorted(tuple(sorted((id(tens._raw_tensor), pos) for tens,pos in prim_list)) for prim_list in values.contractions)),
                                     tuple(tuple(sorted((id(tens._raw_tensor), pos) for tens,pos in prim_list)) for prim_list in values.free_indices) )
        return values
    @staticmethod
    def build(*tensors_indices):
        timings_start()
        scalar, contractions, new_contractions, free_indices, free_indices_as_dict, backend = 1, [], {}, [], {}, None
        for i,(tens, indices) in enumerate(tensors_indices):
            try:
                if backend is None:
                    backend = tens._backend
                if tens._backend is not backend:
                    raise ValueError(f"argument {i} to tensornet.contract._contract has differing backend than those prior")
            except AttributeError:
                raise  TypeError(f"argument {i} to tensornet.contract._contract does not reference a tensornet tensor")
            if len(indices)!=len(tens.shape):
                raise ValueError(f"argument {i} to tensornet.contract._contract has wrong number of indices specified")
            try:
                c = tens._contractions
            except AttributeError:
                tens = copy.copy(tens)    # must be a primitive tensor, so, make distinct object (for hashing), ...
                scalar *= tens._scalar    # ... accumulate scalar factors, ...
                tens._scalar = 1          # ... and renormalize copy (for aesthetics; expected and never referenced in a network_tensor)
            else:
                contractions += c         # inherit contractions from input tensors ...
                scalar *= tens._scalar    # ... and accumulate scalar factors
            for pos,val in enumerate(indices):
                collector = free_indices_as_dict if isinstance(val,int) else new_contractions
                try:
                    if val not in collector:
                        collector[val] = []    # open a list to collect all same-labeled indices
                except:
                    raise TypeError("index label {} (starting from 0) in argument {} to tensornet.contract._contract is not hashable".format(pos,i))
                collector[val] += [(tens, pos)]
        # a helper function to merge networks
        def _resolve_primitive_indices(index_list):    # input indices to be set equal (reduced free indices or contracted together)
            axis_length = None
            prim_list = []
            for tens,pos in index_list:
                if axis_length is None:
                    axis_length = tens.shape[pos]
                elif axis_length!=tens.shape[pos]:
                    raise ValueError("incompatible axis lengths")
                try:
                    prim_list_other = tens._free_indices[pos]
                except AttributeError:
                    prim_list += [(tens, pos)]      # must be a primitive tensor, so just copy over, or ...
                else:
                    prim_list += prim_list_other    # ... else use free index specifications from input tensors
            return prim_list    # contains only references to primitive tensors
        # resolve free indices in terms of primitive tensors
        for i in range(len(free_indices_as_dict)):
            try:
                free_index = free_indices_as_dict[i]
            except:
                raise ValueError("specification of free indices in arguments to tensornet.contract._contract has a gap")
            try:
                free_indices += [_resolve_primitive_indices(free_index)]
            except ValueError:
                raise ValueError("incompatible lengths for reduction to free axis {} (starting from 0) in tensornet.contract._contract".format(i))
        # resolve contracted indices in terms of primitive tensors
        for dummy,contraction in new_contractions.items():
            try:
                contractions += [_resolve_primitive_indices(contraction)]
            except ValueError:
                raise ValueError("incompatible lengths for summation over \"{}\" in tensornet.contract._contract".format(dummy))
        timings_record("contract")
        return struct(backend=backend, scalar=scalar, contractions=contractions, free_indices=free_indices)

    @staticmethod
    def evaluate(the_tensor, tensor_network, primitive_tensor, _backend_contract_path):
        if the_tensor._scalar==0:    # usually a bad test, but in this case we really mean it.  If it is not exactly zero, there is something to do, and zero can happen (ie, a = 0 * b)
            return primitive_tensor(the_tensor._backend, the_tensor._backend.zeros(the_tensor.shape))
        timings_start()
        # It is assumed that all of the tensors in the network are represented by distinct
        # objects, even if they point to the same underlying data.  This is enforced by
        # the contract function, which is the only way for a user to make a tensor_network.
        # Also, the individual primitives have all been forced to have unit scalar by adjusting the overall scalar.
        by_id, contractions, free_indices = the_tensor._hashable
        #
        if _backend_contract_path:    # All the tensors in one big group
            do_scalar_mult = False if the_tensor._scalar==1 else True
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
            do_scalar_mult, do_reduction, target = heuristic(the_tensor._scalar, contraction_groups, index_reduct_groups, shapes)
            timings_record("heuristic")
        #
        if do_scalar_mult or do_reduction:
            timings_start()
            if do_reduction:    # as if do_scalar_mult is False, which it will be if using tensornet contraction path
                # print("do_reduction")
                scalar = the_tensor._scalar
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
                args += [the_tensor._scalar]
                scalar = 1
            timings_record("tensor_network._evaluate")
            #
            timings_start()
            new_tens = primitive_tensor(the_tensor._backend, the_tensor._backend.contract(*args))
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
            return evaluate(tensor_network(the_tensor._backend, scalar, new_contractions, new_free_indices))    # recur
        else:    # must be a single tensor (?) or an outer product
            if len(free_indices)>0:
                timings_start()
                mapping = {}
                for i,free_index in enumerate(free_indices):
                    tens, pos = free_index[0]    # guaranteed to be only one entry per index now
                    if tens not in mapping:
                        mapping[tens] = [None]*len(by_id[tens].shape)
                    mapping[tens][pos] = i
                if the_tensor._scalar!=1:                                     # eventually deprecate.
                    raise RuntimeError("scalar should be 1, right?")    # pretty sure this has to be true
                args = [(by_id[tens]._raw_tensor, *indices) for tens,indices in mapping.items()]    # the_tensor._scalar is 1 by now?
                timings_record("tensor_network._evaluate")
                timings_start()
                Z = primitive_tensor(the_tensor._backend, the_tensor._backend.contract(*args))                     # bottom out (might give a 0-dim tensor; this is intended)
                timings_record("backend.contract")
                return Z
            else:    # there is nothing left but the scalar
                return primitive_tensor(the_tensor._backend, the_tensor._backend.scalar_tensor(the_tensor._scalar))

    @staticmethod
    def subscript(the_tensor, indices, tensor_network):
        full = slice(None)    # the slice produced by [:] with no limits
        scalar = the_tensor._scalar
        by_id, contractions, free_indices = the_tensor._hashable
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
        new = tensor_network(the_tensor._backend, scalar, new_contractions, new_free_indices)
        if len(new_free_indices)==0:
            global _warned
            if not _warned:
                print("Accessing individual elements of a tensor network is really inefficient.  Consider alternatives.")    # How user can suppress this altogether?
                _warned = True
            return scalar_value(new)
        else:
            return new

    @staticmethod
    def resolve(scalar, tensors_indices, tensor_network):
        # This function turns a contractions of sums into sums of contractions and passes the terms off to tensor_network
        # Implicit type checking is essentially delegated to tensor_network.
        outer_terms = [[]]
        for tens, indices in tensors_indices:
            try:
                inner_factors = tens._terms
            except:
                for outer_term in outer_terms:
                    outer_term += [(tens, indices)]
            else:
                new_outer_terms = []
                for outer_term in outer_terms:
                    for inner_factor in inner_factors:
                        new_outer_terms += [outer_term + [(inner_factor, indices)]]
                outer_terms = new_outer_terms
        terms = []
        for outer_term in outer_terms:
            terms += [scalar * tensor_network(**logic.build(*outer_term))]
        #
        #result_hashes = sorted((term._result_hash, i) for i,term in enumerate(terms))
        #terms, terms_ = [], terms
        #previous = None
        #for result_hash,i in result_hashes:
        #    if result_hash==previous:
        #        terms[-1]._scalar += terms_[i]._scalar
        #    else:
        #        terms += [terms_[i]]
        #    previous = result_hash
        #
        if len(terms)==1:
            return terms[0]
        else:
            the_sum = terms[0] + terms[1]
            for term in terms[2:]:
                the_sum += term
            return the_sum

    @staticmethod
    def get_shape(tensors_indices, check_backend):
        # dims_catalog has the following structure
        #   {index: struct(length, priors=[(arg,dim), ...]), ...}
        # where index is the int or "letter" identifier of the eventual free or contraction index, repsectively,
        # length is the previously recorded length of the associated dimension, and priors records where such prior
        # lengths were found, in which arg is the positional identifier of the tensor in the tensors_indices input,
        # and dim is the position of the respective index on that tensor.  A given tensor (arg) may show up
        # multiple times in a given priors array.
        dims_catalog = {}
        for arg,(tensor,indices) in enumerate(tensors_indices):
            check_backend(tensor)    # otherwise raises exception
            if len(indices)!=len(tensor.shape):
                raise ValueError(f"argument {arg} to tensornet.contraction_expression has wrong number of indices specified ({len(indices)} given, {len(tensor.shape)} expected).")
            for dim,index in enumerate(indices):
                if index not in dims_catalog:
                    dims_catalog[index] = struct(length=tensor.shape[dim], priors=[(arg,dim)])
                else:
                    if tensor.shape[dim]==dims_catalog[index].length:
                        dims_catalog[index].priors += [(arg,dim)]
                    else:
                        priors = "\n".join([f"dimension {dim_} of argument {arg_}" for arg_,dim_ in dims_catalog[index].priors])
                        raise ContractionError(f"\ncontraction_expression:  error for identifier \"{index}\" (all enumerations start at 0)\ndimension {dim} of argument {arg} has incompatible length ({tensor.shape[dim]}) for contraction or reduction with dimension(s) of length {dims_catalog[index].length} in:\n{priors}")
        free_indices = sorted([index for index in dims_catalog if isinstance(index,int)])
        shape = [] if (len(free_indices)==0) else [None]*(free_indices[-1]+1)
        for index in free_indices:
            shape[index] = dims_catalog[index].length
        return shape
