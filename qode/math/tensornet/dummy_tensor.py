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



class tensor(object):
    def __init__(self, scalar, name, shape):
        self.name    = name
        self.shape   = tuple(shape)
        self.scalar  = scalar
    def item(self):
        if len(self.shape)!=0:
            raise ValueError("can only convert an array with no indices to a scalar")
        return self.scalar
    def __str__(self):
        return str(self.scalar) + " * " + self.name + str(self.shape)



def contract(*tensor_factors):
    name, scalar, contractions, free_indices = "", 1, {}, {}
    local_scalar = 1
    contract_string = ""
    for i,factor in enumerate(tensor_factors):
        try:
            tens, *indices = factor
        except:
            try:
                scalar *= factor
                local_scalar *= factor
            except:
                raise TypeError("argument {} to dummy_tensor.contract is malformed".format(i))
        else:
            try:
                name   += tens.name
                dims    = tens.shape
                scalar *= tens.scalar
            except:
                raise  TypeError("argument {} to dummy_tensor.contract does not reference a dummy tensor".format(i))
            if len(indices)!=len(dims):
                raise ValueError("argument {} to dummy_tensor.contract has wrong number of indices specified".format(i))
            for pos,val in enumerate(indices):
                collector = free_indices if isinstance(val,int) else contractions
                try:
                    if val not in collector:
                        collector[val] = []    # open a list to collect all same-labeled indices
                except:
                    raise TypeError("index label {} (starting from 0) in argument {} to dummy_tensor.contract is not hashable".format(pos,i))
                collector[val] += [(tens, pos)]
            if tens.scalar!=1:  contract_string += str(tens.scalar) + "*"
            contract_string += tens.name + "_"
            for x in indices:  contract_string += str(x)
            contract_string += str(tens.shape).replace(" ", "")
            contract_string += " "
    if local_scalar!=1:
        contract_string = str(local_scalar) + " * " + contract_string
    #
    def _check_indices(prim_list):    # input indices to be set equal (reduced free indices or contracted together)
        axis_length = None
        for tens,pos in prim_list:
            if axis_length is None:
                axis_length = tens.shape[pos]
            elif axis_length!=tens.shape[pos]:
                raise ValueError("incompatible axis lengths")
    #
    for i in range(len(free_indices)):
        try:
            free_index = free_indices[i]
        except:
            raise ValueError("specification of free indices in arguments to dummy_tensor.contract has a gap")
        try:
            _check_indices(free_index)
        except ValueError:
            raise ValueError("incompatible lengths for reduction to free axis {} (starting from 0) in dummy_tensor.contract".format(i))
    # 
    for dummy,contraction in contractions.items():
        try:
            _check_indices(contraction)
        except ValueError:
            raise ValueError("incompatible lengths for summation over \"{}\" in dummy_tensor.contract".format(dummy))
    shape = []
    free = ""
    for i in range(len(free_indices)):
        tens0, pos0 = free_indices[i][0]
        shape  += [tens0.shape[pos0]]
        free += str(i)
    free += str(tuple(shape)).replace(" ", "")
    contract_string = name + "_" + free + " = " + contract_string
    if scalar!=1:
        contract_string = str(scalar) + "*" + contract_string
    print(contract_string)
    return tensor(scalar, name, shape)
