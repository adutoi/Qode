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

from primitive_tensor import thin_wrapper
from tensor_network   import tensor_network, tensor_network_sum



# This handles the case where the tensor_factors are single terms (not sums) for contract() below
def _contract(*tensor_factors):
    # collect arguments
    scalar, contractions, new_contractions, free_indices, free_indices_as_dict, backend = 1, [], {}, [], {}, None
    for i,factor in enumerate(tensor_factors):
        try:
            tens, *indices = factor
        except:
            try:
                scalar *= factor
            except:
                raise TypeError("argument {} to contract._contract is malformed".format(i))
        else:
            try:
                shape = tens.tensornet_shape
            except:
                raise  TypeError("argument {} to contract._contract does not reference a tensornet tensor".format(i))
            if backend is None:
                backend = tens.tensornet_backend
            if tens.tensornet_backend is not backend:
                raise ValueError("argument {} to contract._contract has differing backend than those prior".format(i))
            if len(indices)!=len(shape):
                raise ValueError("argument {} to contract._contract has wrong number of indices specified".format(i))
            try:
                c = tens.contractions
            except AttributeError:
                tens = thin_wrapper(tens)      # must be a primitive tensor, so wrap and skip other stuff
            else:
                contractions += c              # inherit contractions from input tensors ...
                scalar       *= tens.scalar    # ... and accumulate scalar factors
            for pos,val in enumerate(indices):
                collector = free_indices_as_dict if isinstance(val,int) else new_contractions
                try:
                    if val not in collector:
                        collector[val] = []    # open a list to collect all same-labeled indices
                except:
                    raise TypeError("index label {} (starting from 0) in argument {} to contract._contract is not hashable".format(pos,i))
                collector[val] += [(tens, pos)]
    # a helper function to merge networks
    def _resolve_primitive_indices(index_list):    # input indices to be set equal (reduced free indices or contracted together)
        axis_length = None
        prim_list = []
        for tens,pos in index_list:
            if axis_length is None:
                axis_length = tens.shape[pos]     # at this point, it is either a wrapped ...
            elif axis_length!=tens.shape[pos]:    # ... primitive_tensor or a tensor_network
                raise ValueError("incompatible axis lengths")
            try:
                prim_list_other = tens.free_indices[pos]
            except AttributeError:
                prim_list += [(tens, pos)]      # must be a primitive tensor, so just copy over, ...
            else:
                prim_list += prim_list_other    # ... or use free index specifications from input tensors
        return prim_list    # contains only references to primitive tensors
    # resolve free indices in terms of primitive tensors
    for i in range(len(free_indices_as_dict)):
        try:
            free_index = free_indices_as_dict[i]
        except:
            raise ValueError("specification of free indices in arguments to contract._contract has a gap")
        try:
            free_indices += [_resolve_primitive_indices(free_index)]
        except ValueError:
            raise ValueError("incompatible lengths for reduction to free axis {} (starting from 0) in contract._contract".format(i))
    # resolve contracted indices in terms of primitive tensors
    for dummy,contraction in new_contractions.items():
        try:
            contractions += [_resolve_primitive_indices(contraction)]
        except ValueError:
            raise ValueError("incompatible lengths for summation over \"{}\" in contract._contract".format(dummy))
    return tensor_network(scalar, contractions, free_indices, backend)



# Usage: new_tens_network = contract((tens1,"p",0), (tens2,"p",1), scalar, (tens3,"p","p",1), (tens4,2,3))
# This function turns a contractions of sums into sums of contractions and passes the terms off to _contact() above
# Implicit type checking is essentially delegated to _contract().
def contract(*tensor_factors):
    outer_terms = [[]]
    for factor in tensor_factors:
        try:
            tens, *indices = factor
            inner_factors = tens.tensor_terms
        except:
            for outer_term in outer_terms:
                outer_term += [factor]
        else:
            new_outer_terms = []
            for outer_term in outer_terms:
                for inner_factor in inner_factors:
                    new_outer_terms += [outer_term + [(inner_factor, *indices)]]
            outer_terms = new_outer_terms
    tensor_terms = []
    for outer_term in outer_terms:
        tensor_terms += [_contract(*outer_term)]
    if len(tensor_terms)==1:
        return tensor_terms[0]
    else:
        return tensor_network_sum(tensor_terms)
