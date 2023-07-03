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
from qode.math.tensornet import np_tensor, contract, evaluate, increment, extract, scalar_value

p,q,r,s = 'pqrs'

# Some original "raw tensors"
A_ = numpy.random.random((10, 10))
B_ = numpy.random.random((10,  8))
C_ = numpy.random.random(( 8, 10))
D_ = numpy.random.random(( 8,  8))

# Prepare "primitive tensors" for use by tensornet
A = np_tensor(A_)
B = np_tensor(B_)
C = np_tensor(C_)
D = np_tensor(D_)

# - The top-level utility is this contract function.  It does no computation at this point,
#   but successive calls build up a "tensor network" of contractions.
# - Each tensor argument comes as a tuple whose 0th element is the tensor itself (primitive
#   or another network), and the remaining members are either (1) a label whose recurrance
#   specifies indices to be contracted, or (2) an integer whose order specifies the
#   permutation of the free indices in the contracted product.
# - Contract will raise an exception if a logical error occurs with the shapes, like
#   contracting indices of incompatible lengths.
# - No limit on the number of arguments and the result will not depend in any way (not even
#   efficiency) on how they are ordered.
# - Scalar arguments can also be given in place of tensors.
# - A contraction "label" can be any hashable object that is not an integer (most probably
#   a letter or a string).
AB     = contract((A,0,p), (B,p,1))
CD     = contract((C,p,1), (D,0,p))
ABCD1  = contract((AB,0,p), (CD,p,1))
ABCD1 *= 2

# Let us test the above.  The evaluate() function will cause the contraction to be performed
# returning a new primitive tensor containing the result (without altering the original).
# The user can use this to force the computations of intermediates that might be used in
# multiple further expressions.  The extract() function calls evaluate() and then returns
# the result as a raw tensor.
test = 2 * numpy.einsum("pr,rt,sq,ts->pq", A_, B_, C_, D_)
print("relative error in first  check:", numpy.linalg.norm(extract(ABCD1) - test) / numpy.linalg.norm(test))

# We could have also contructed ABCD in once line as
ABCD2 = contract(2, (A,0,p), (B,p,q), (C,r,1), (D,q,r))
print("relative error in second check:", numpy.linalg.norm(extract(ABCD2) - test) / numpy.linalg.norm(test))

# And just to illustrate that +, - and * (with scalars) all work.  The result of the multiplications (which are fast)
# another tensor network whose internal scalar factor is different.  Similarly, the result of the subtraction is an
# abstracted structure that keeps track of the operation but performs no computation.
zero = 3*ABCD1 - 3*ABCD2
print("relative error in third  check:", numpy.linalg.norm(extract(zero)) / numpy.linalg.norm(extract(ABCD1)))

# This illustrates that there are not many restrictions on the order of the arguments or the placement of the indices.
# The redundant "1" means that the latter index of the contracted product is resolved by setting these two equal to
# each other.  Also single as well as multiple (not just double) occurances of labels can be contracted.
# Note that the same index can occur multiply on any given tensor, and that outer products (A shares not indices)
# are allowed
ABCD3 = contract((A,q,0), (B,1,p), 3, (C,p,1), (D,p,p))
test = 3 * numpy.einsum("qr,sp,ps,pp->rs", A_, B_, C_, D_)
print("relative error in fourth check:", numpy.linalg.norm(extract(ABCD3) - test) / numpy.linalg.norm(test))

# But the real utility is for something like this.  E represents 4-index tensor, but its contractions with vectors
# can be done faster if done internally as contractions with the factors
E = contract((A,0,1), (D,2,3))
V_ = numpy.random.random((10,))
U_ = numpy.random.random((8,))
V = np_tensor(V_)
U = np_tensor(U_)
test = numpy.einsum("pq,rs,p,q,r,s->", A_, D_, V_, V_, U_, U_)
VVEUU = contract((E,p,q,r,s), (V,p), (V,q), (U,r), (U,s))
print("relative error in fifth  check:", numpy.linalg.norm(extract(VVEUU) - test) / numpy.linalg.norm(test))

# indexing and slicing
#from .tensors  import primitive_tensor_wrapper, tensor_sum    # tensor_sum() can initialize an empty accumulator for += use
#from .backends import dummy_backend, numpy_backend
#dummy_tensor





