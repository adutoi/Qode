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

import numpy
import numpy_backend
from contract import contract
from tensors import primitive_tensor_wrapper
from base import evaluate, raw

tensor = primitive_tensor_wrapper(numpy_backend)

p,q,r,s,t,u = 'pqrstu'


A = numpy.random.random((10, 10))
B = numpy.random.random((10, 10))
C = A @ B

A_ = tensor(A)
B_ = tensor(B)

pdt = contract((A_,0,p), (B_,p,1))
C_ = raw(evaluate(pdt))

norm  = numpy.linalg.norm(C)
error = numpy.linalg.norm(C-C_)
print(norm, error)
