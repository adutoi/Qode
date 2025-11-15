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

# The main purpose of this class is that it announces to the user what contractions are
# being done by the backend.  But for the purposes of having the code run and do something
# sensical when numerical results are required, all non-zero-dimensional tensors are assumed
# to be zero.

class dummy(object):
    def __init__(self, scalar, name, shape):
        self.name    = name
        self.shape   = tuple(shape)
        self.scalar  = scalar
    def data(self):
        if len(self.shape)==0:
            return self.scalar
        tensor = [0] * shape[-1]
        for d in reversed(shape)[1:]:
            tensor = [inner] * d
        return tensor

class functions(object):
    name = "dummy"
    @staticmethod
    def scalar_value(tensor):
        if len(tensor.shape)!=0:
            raise RuntimeError(f"cannot convert {len(tensor.shape)}-dimensional tensor")
        return tensor.scalar
    @staticmethod
    def scalar_tensor(scalar):
        return dummy(1, f"({scalar})", tuple())
    @staticmethod
    def shape(tensor):
        return tensor.shape
    @staticmethod
    def zeros(shape):
        return dummy(1, f"(0)", shape)
    @staticmethod
    def increment(tensor, delta):
        if tensor.shape!=delta.shape:
            raise ValueError("cannot add tensors of different shapes")
        tensor.name = f"({tensor.scalar}*{tensor.name} + {delta.scalar}*{delta.name})"
        tensor.scalar = 1
        return
    @staticmethod
    def mult(scalar, tensor):
        return dummy(scalar*tensor.scalar, tensor.name, tensor.shape)
    @staticmethod
    def slice(tensor, indices):
        if len(indices)!=len(tensor.shape):
            raise ValueError(f"{len(indices)} were given to obtain element of {len(tensor.shape)}-dimensional tensor")
        return 0
    @staticmethod
    def str(tensor):
        shape = "x".join(str(_) for _ in tensor.shape)
        return f"{tensor.scalar} * {tensor.name}[{shape}]"
    @staticmethod
    def contract(*tensor_factors):
        ####
        # args = []
        # for factor in tensor_factors:
        #     try:
        #         tensor,*indices = factor
        #     except:
        #         args += [factor]
        #     else:
        #         args += [(id(tensor),*indices)]
        # print("backend.contract called with", *args)
        ####
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
                    raise TypeError(f"argument {i} to dummy_backend.functions.contract is malformed")
            else:
                try:
                    name   += tens.name
                    scalar *= tens.scalar
                    if len(indices)!=len(tens.shape):
                        raise ValueError(f"argument {i} to dummy_backend.functions.contract has wrong number of indices specified")
                except AttributeError:
                    raise TypeError(f"argument {i} to dummy_backend.functions.contract does not reference a dummy tensor")
                for pos,val in enumerate(indices):
                    collector = free_indices if isinstance(val,int) else contractions
                    try:
                        if val not in collector:
                            collector[val] = []    # open a list to collect all same-labeled indices
                    except:
                        raise TypeError("index label {pos} (starting from 0) in argument {i} to dummy_tensor.contract is not hashable")
                    collector[val] += [(tens, pos)]
                if tens.scalar!=1:  contract_string += str(tens.scalar) + "*"
                contract_string += tens.name + "_" + "".join(str(_) for _ in indices) + "[" + "x".join(str(_) for _ in tens.shape) + "] "
        if local_scalar!=1:
            contract_string = str(local_scalar) + " * " + contract_string
        #
        def _check_indices(prim_list):    # input indices to be set equal (reduced free indices or contracted together)
            axis_length = None
            for tens,pos in prim_list:
                if axis_length is None:
                    axis_length = tens.shape[pos]
                if tens.shape[pos]!=axis_length:
                    raise ValueError("incompatible axis lengths")
        #
        for i in range(len(free_indices)):
            try:
                free_index = free_indices[i]
            except:
                raise ValueError("specification of free indices in arguments to dummy_backend.functions.contract has a gap")
            try:
                _check_indices(free_index)
            except ValueError:
                raise ValueError("incompatible lengths for reduction to free axis {i} (starting from 0) in dummy_backend.functions.contract")
        # 
        for idx,contraction in contractions.items():
            try:
                _check_indices(contraction)
            except ValueError:
                raise ValueError("incompatible lengths for summation over \"{idx}\" in dummy_backend.functions.contract")
        shape = []
        #free = "_"
        for i in range(len(free_indices)):
            tens0, pos0 = free_indices[i][0]
            shape  += [tens0.shape[pos0]]
        #    free += str(i)
        free = "[" + "x".join(str(_) for _ in shape) + "]"
        contract_string = name + free + " = " + contract_string
        if scalar!=1:
            contract_string = str(scalar) + "*" + contract_string
        print(contract_string)
        return dummy(scalar, name, shape)
