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
from qode.math.tensornet import dummy_tensor, evaluate
from qode.math.tensornet.backends.dummy_backend import dummy

p, q, r, s = "pqrs"

A  = dummy_tensor.init(dummy(1,"A",(5,4)))
P  = dummy_tensor.init(dummy(1,"P",(5,5)))
Q  = dummy_tensor.init(dummy(1,"Q",(4,4)))
Z0 = dummy_tensor.init(dummy(1,"Z0",(5,)))
Z1 = dummy_tensor.init(dummy(1,"Z1",(4,)))
Z2 = dummy_tensor.init(dummy(1,"Z2",(5,)))
Z3 = dummy_tensor.init(dummy(1,"Z3",(4,)))

AA = A(0,1) @ A(2,3) - A(0,3) @ A(2,1)

PZ0 = P(0,p) @ Z0(p)
QZ1 = Q(0,p) @ Z1(p)
PZ2 = P(0,p) @ Z2(p)
QZ3 = Q(0,p) @ Z3(p)

dot = 2 * AA(p,q,r,s) @ PZ0(p) @ QZ1(q) @ PZ2(r) @ QZ3(s)

X = evaluate(dot)
print(X)
