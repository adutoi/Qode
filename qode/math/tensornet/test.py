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

import dummy_tensor
import dummy_backend
from contract import contract
from primitive_tensor import primitive_tensor_wrapper

tensor = primitive_tensor_wrapper(dummy_backend)



A = tensor(dummy_tensor.tensor(1,"A",(5,4)))
#AA = contract((A,0,1),(A,2,3)) + contract((A,0,3),(A,2,1))
AA = contract((A,0,1),(A,2,3))

P = tensor(dummy_tensor.tensor(1,"P",(5,5)))
Q = tensor(dummy_tensor.tensor(1,"Q",(4,4)))

Z0 = tensor(dummy_tensor.tensor(1,"Z0",(5,)))
Z1 = tensor(dummy_tensor.tensor(1,"Z1",(4,)))
Z2 = tensor(dummy_tensor.tensor(1,"Z2",(5,)))
Z3 = tensor(dummy_tensor.tensor(1,"Z3",(4,)))

PZ0 = contract((P,0,"p"),(Z0,"p"))
QZ1 = contract((Q,0,"p"),(Z1,"p"))
PZ2 = contract((P,0,"p"),(Z2,"p"))
QZ3 = contract((Q,0,"p"),(Z3,"p"))

dot = contract(2, (AA,"p","q","r","s"), (PZ0,"p"), (QZ1,"q"), (PZ2,"r"), (QZ3,"s"))

X = dot.tensornet_evaluate()
print(X)
