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
from qode.math.tensornet import np_tensor, tensor_sum, evaluate, increment, extract, scalar_value

p,q,r,s = 'pqrs'                   # lower the number of quotes we need to type

def check(xpt, ref, raw=False):    # package up the repetitive logic of checking accuracy
    if raw:
        return (xpt - ref) / ref
    else:
        return numpy.linalg.norm(xpt - ref) / numpy.linalg.norm(ref)



# Some original "raw tensors"
A_ = numpy.random.random((10, 10))
B_ = numpy.random.random((10,  8))
C_ = numpy.random.random(( 8, 10))
D_ = numpy.random.random(( 8,  8))
U_ = numpy.random.random((8,))
V_ = numpy.random.random((10,))

# Prepare "primitive tensors" for use by tensornet with wrapping function for numpy (np)
A = np_tensor(A_)
B = np_tensor(B_)
C = np_tensor(C_)
D = np_tensor(D_)
U = np_tensor(U_)
V = np_tensor(V_)

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
AB     = A(0,p) * B(p,1)
CD     = C(p,1) * D(0,p)
ABCD1  = AB(0,p) * CD(p,1)
ABCD1 *= 2

# Let us test the above.  The evaluate() function will cause the contraction to be performed
# returning a new primitive tensor containing the result (without altering the original).
# The user can use this to force the computations of intermediates that might be used in
# multiple further expressions.  The extract() function calls evaluate() [see below]
# and then returns the result as a raw tensor.
test = 2 * numpy.einsum("pr,rt,sq,ts->pq", A_, B_, C_, D_)
print("relative error in check  1:", check(extract(ABCD1), test))

# We could have also contructed ABCD in once line as
ABCD2 = 2 * A(0,p) * B(p,q) * C(r,1) * D(q,r)
print("relative error in check  2:", check(extract(ABCD2), test))

# And just to illustrate that +, - and * (with scalars) all work.  The result of the multiplications (which are fast)
# another tensor network whose internal scalar factor is different.  Similarly, the result of the subtraction is an
# abstracted structure that keeps track of the operation but performs no computation.
zero = 3*ABCD1 - 3*ABCD2
print("relative error in check  3:", numpy.linalg.norm(extract(zero)) / numpy.linalg.norm(extract(ABCD1)))

# This illustrates that there are not many restrictions on the order of the arguments or the placement of the indices.
# The redundant "1" means that the latter index of the contracted product is resolved by setting these two equal to
# each other.  Also single as well as multiple (not just double) occurances of labels can be contracted.
# Note that the same index can occur multiply on any given tensor, and that outer products (A shares not indices)
# are allowed
ABCD3 = A(q,0) * B(1,p) * 3 * C(p,1) * D(p,p)
test = 3 * numpy.einsum("qr,sp,ps,pp->rs", A_, B_, C_, D_)
print("relative error in check  4:", check(extract(ABCD3), test))

# But the real utility is for something like this.  E represents 4-index tensor, but its contractions with vectors
# can be done faster if done internally as contractions with the factors
E = A(0,1) * D(2,3)
test = numpy.einsum("pq,rs,p,q,r,s->", A_, D_, V_, V_, U_, U_).item()
VVEUU = scalar_value(E(p,q,r,s) * V(p) * V(q) * U(r) * U(s))
print("relative error in check  5:", (VVEUU - test) / test)

# You can take elements of a primitive tensor, ...
print("relative error in check  6:", (A[0,0] - A_[0,0]) / A_[0,0])
# ... or of a tensor networ.  This will reevaluate everything for every element, so it warns you the first
# time about this being inefficient.  It does not evaluate the whole tensor however, but first restricts
# the indices in the primitives before contraction.  So, in this particular case, it actually is not inefficent.
print("relative error in check  7:", (E[0,0,0,0] - A_[0,0]*D_[0,0]) / (A_[0,0]*D_[0,0]))
print("relative error in check  8:", (E[1,1,1,1] - A_[1,1]*D_[1,1]) / (A_[1,1]*D_[1,1]))

# Slicing also works, ...
print("relative error in check  9:", check(extract(A[0, 1:5]), A_[0, 1:5]))
# ... within contractions ...
AC = A[:,:8](0,p) * C(p,1)
test = numpy.einsum("pq,qr->pr", A_[:,:8], C_)
print("relative error in check 10:", check(extract(AC), test))
# ... and on tensor networks (where again it first acts to slice the primitives in the network immediately (before evaluation requested)
test = numpy.einsum("p,q->pq", A_[:,0], D_[:,1])
print("relative error in check 11:", check(extract(E[:,0,:,1]), test))

# One thing that is good to be aware of is that, if one contracts a sum (such as F below), the contractions
# of the terms are always done first and then added (most likely desired and far simpler dispatching algorithm), ...
F = A + 5*B(0,p)*C(p,1)
test = numpy.einsum("p,pq,q->", V_, A_, V_) +  5 * numpy.einsum("p,pq,qr,r->", V_, B_, C_, V_)
print("relative error in check 12:", check(extract(V(p) * F(p,q) * V(q)), test))
# ... but the user has control over this by forcing evaluation of the intermediate first.
# The evaluate() function performs all internal contractions and index reductions and returns
# a primitive (not a raw) tensor.
print("relative error in check 13:", check(extract(V(p) * evaluate(F)(p,q) * V(q)), test))

# Should one need to loop over terms to build a sum, and empty accumulator can be started as so:
F = tensor_sum()
for term in [A, 5*B(0,p)*C(p,1)]:    # pretend these would be generated by some algorithm:
    F += term
print("relative error in check 14:", check(extract(V(p) * F(p,q) * V(q)), test))

# Another way to do this, if you know you want a concrete evaluation immediately.
# This does not yet save memory with an in-place build, but the point is to leave the
# door open for that (presently just a thin wrapper around +=extract(...)).
F_ = numpy.zeros(A.shape)
for term in [A, 5*B(0,p)*C(p,1)]:
    increment(F_, term)
print("relative error in check 15:", check(numpy.einsum("pq,p,q->", F_, V_, V_), test)), 



# Finally, the concept of a backend is local to the tensors.  Rather than have a global setting,
# any two tensors with the same back end may enter into tensornet expressions with each other.
# The long way of doing this for a numpy backend would be as follows.  The default is not to make
# a copy of the original tensor data, but simply to wrap it.  As long as the user has no intention
# of chaning the backing data, tensornet will not either, and it is more efficient, but if you 
# are getting weird results and want to be sure, or you just know you can't promise not to
# modify the original, the use the copy_data flag.
# The np_tensor function used above, is equivalent to the first option below and is provided
# for convenience. 
from qode.math.tensornet import primitive_tensor_wrapper, numpy_backend
np_tensor_ = primitive_tensor_wrapper(numpy_backend)                    # this ...
np_tensor_ = primitive_tensor_wrapper(numpy_backend, copy_data=True)    # ... or this
A = np_tensor_(A_)

# One can see now that providing a substitute for numpy_backend above is all one needs for a custom
# backend.  This should be a module that has the following 3 functions defined in it.
#
# def copy_data(tensor):
#    # Returns a tensor referencing independent numerical data that is a copy of that of tensor
#    # for numpy this just returns numpy.array(tensor)
#
# def scalar_value(tensor):
#    # Assuming the tensor has no indices left, return the internal value as a scalar type (like float)
#    # For numpy, this just returns tensor.item()
#
# def scalar_tensor(scalar):
#    # Return the scalar in the data type of a 0-index tensor
#    # For numpy, this just returns numpy.array(scalar)
# def shape(tensor):
#    # Give the shape (dimensions) of tensor as a tuple of integers
#    # For numpy, this just returns tensor.shape
#
# def contract(*tensor_factors):
#    # The arguments of this function are exactly of the same format as the tensornet contract()
#    # function.  The difference is that it is likely to be called with a much smaller number
#    # of arguments, as the tensornet machiney grinds down the contraction network stepwise.
#    # A user is free to make this as general or as limited as necessary.  For example, if only
#    # pairwise contractions of indices on two distinct tensors is every necessary (by-and-large
#    # the dominant use case), then this routine will only every be called as such.  The implementor
#    # can just raise a NotImplementedError if a different use case is ever given.
#    # For numpy, the general use case can be mapped onto the arguments of a valid einsum call.
#
# Additionally, the raw tensor data type should understand the operations +, -, *, [], where
# addition/subtraction is to other tensors of the same shape, and * means with a scalar.  The
# item operator should be able to handle tuples that include slices.

# PS - In this directory, you see a dummy_tensor module, which is supported by the dummy_backend
# module provided by tensornet.  This does not do any math but just prints out the contractions
# that are being performed.  Originally used to debug tensornet, it might come in handy for users
# to verify that their code is being executed as expected (and possibly find problems with the
# contraction-order heuristic).



### Pure outer product
