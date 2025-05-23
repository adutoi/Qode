#    (C) Copyright 2018, 2023, 2025 Anthony D. Dutoi and Yuhong Liu
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
import numpy
from ...util.PyC import Double
from . import field_op



class _auxilliary(object):
    def __init__(self, parent):
        self.parent = parent
    def basis_vec(self, occupied):
        v = numpy.zeros(len(self.parent.configs), dtype=Double.numpy, order="C")
        i = field_op.find_index_by_occ(occupied, self.parent.configs)
        if i==-1:
            raise RuntimeError("requested basis state is not a member of the set of allowed configurations")
        v[i] = 1
        return v
    def complete_basis(self):
        n_config = len(self.parent.configs)
        basis = []
        for i in range(n_config):
            v = numpy.zeros(n_config, dtype=Double.numpy, order="C")
            v[i] = 1
            basis += [v]
        return basis
        


class CI_space_traits(object):
    def __init__(self, configs):
        self.field = Double.numpy
        self.aux = _auxilliary(self)
        self.configs = field_op.packed_configs(configs)
    def check_member(self,v):
        pass
    def check_lin_op(self,op):
        return False
    @staticmethod
    def copy(v):
        return v.copy()
    @staticmethod
    def scale(c,v):
        v *= c
    @staticmethod
    def add_to(v,w,c=1):
        if c==1:  v += w
        else:     v += c*w
    @staticmethod
    def dot(v,w):
        return v.dot(w)
    def act_on_vec(self, op, v):
        return op(v, self.configs)
    def back_act_on_vec(self, v, op):
        return op(v, self.configs)
    def act_on_vec_block(self, op, v_block):
        return [ op(v, self.configs) for v in v_block ]
    def back_act_on_vec_block(self, v_block, op):
        return [ op(v, self.configs) for v in v_block ]
    @staticmethod
    def dot_vec_blocks(v_block,w_block):
        return numpy.array([[CI_space_traits.dot(v,w) for v in v_block] for w in w_block])
