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

### If everything is working, running this should print out a bunch of random numbers that are
### < ~1e-15, and a single instance of a warning about extracting individual elements.  The
### usage is explained in the comments below.

import sys
import math
import numpy
import torch
import tensorly
from qode.math.tensornet import np_tensor, tl_tensor, tensor_sum, evaluate, increment, raw, scalar_value, contract

p,q,r,s = 'pqrs'        # lower the number of quotes we need to type

#######################################################################################
### In tensornet, backend is a tensor-local concept, but this file does not make use of
### multiple backends, so here is an example of how one makes backend a "global" choice
### (by wrapping the user functions that deal with raw tensors).

#tensorly.set_backend("pytorch")    # uncomment if using tensorly with pytorch backend

#zeros, einsum, prim_tensor =    numpy.zeros,    numpy.einsum, np_tensor    # uncomment if numpy directly as backend
zeros, einsum, prim_tensor = tensorly.zeros, tensorly.einsum, tl_tensor    # uncomment if using tensorly (numpy or pytorch)

def check(xpt, ref):    # package up the repetitive logic of checking accuracy
    #return numpy.linalg.norm(xpt - ref) / numpy.linalg.norm(ref)    # uncomment if numpy directly as backend
    return tensorly.norm(xpt - ref) / tensorly.norm(ref)            # uncomment if using tensorly (numpy or pytorch) ... btw, what tensorly calls l^2 is Frobenius.

def random_tensor(*shape):
    #return numpy.random.random(shape)                     # uncomment if using numpy directly as backend
    return tensorly.tensor(numpy.random.random(shape))    # uncomment if using tensorly with numpy backend
    #return tensorly.tensor(torch.rand(*shape))            # uncomment if using tensorly with pytorch backend

#######################################################################################



    ### The Basics

# Some original "raw" tensors
A_ = random_tensor(10, 10)
B_ = random_tensor(10,  8)
C_ = random_tensor( 8, 10)
D_ = random_tensor( 8,  8)
U_ = random_tensor(8)
V_ = random_tensor(10)

# Prepare "primitive" tensors for tensornet using wrapping function for raw arrays
# The original data is never modified by tensornet calls and functions.
A = prim_tensor(A_)
B = prim_tensor(B_)
C = prim_tensor(C_)
D = prim_tensor(D_)
U = prim_tensor(U_)
V = prim_tensor(V_)

# - The top-level utility is this contract function.  It does no computation at this point,
#   but successive calls build up a tensor "network" of contractions.
# - Each argument to contract is generated by a call on the input tensor (primitive
#   or network), whose own arguments are either (1) a label whose recurrance
#   specifies indices to be contracted, or (2) an integer whose order specifies the
#   permutation of the free indices in the contracted product.
# - Contract will raise an exception if a logical error occurs with the shapes, like
#   contracting indices of incompatible lengths.
# - There is no limit on the number of arguments and (to within machine precision)
#   the result will not depend in any way on how they are ordered (not even efficiency).
# - Scalar arguments can also be given in place of tensors, which multipy the result.
# - A contraction "label" (often a letter) can be any hashable object that is not an integer.
AB1    = contract(A(0,p), B(p,1))
CD1    = contract(C(p,1), D(0,p))
ABCD1  = contract(AB1(0,p), CD1(p,1))
ABCD1 *= 2

# The above can also be done by using the @ operator to separate the input tensors.  Just using
# the call "()" notation on a tensor already eventually implies that the contract function will
# be called, and the @ operator collects arguments for this.  This is encouraged for easy reading,
# but it is the least well tested part of the module (so maybe check against explicit contract()
# if you suspect a problem).  More will be said below, but one side effect of this syntax is
# that the range checking is postponed because the calling of contract() is delayed ("lazy").
AB2    =   A(0,p) @ B(p,1)
CD2    =   C(p,1) @ D(0,p)
ABCD2  = AB2(0,p) @ CD2(p,1)
ABCD2 *= 2

# Let us test the above.  The evaluate() function (more on this further below) will cause the
# actual numerical contractions to be performed (in an optimized order) returning a new primitive
# tensor containing the result.  The user can use evaluate() to force the computations of
# intermediates that might be used in multiple further expressions, for example.  The raw()
# function immediately below calls evaluate() and subsequently returns the result as a raw tensor,
# which we can test against by contracting the orginal data.
test = 2 * einsum("pr,rt,sq,ts->pq", A_, B_, C_, D_)
print("relative error in check  1:", check(raw(ABCD1), test))
print("relative error in check  2:", check(raw(ABCD2), test))

# We could have also contructed ABCD in once line as.  Note that scalar multiplication uses the
# familiar * operator.
ABCD3 = 2 * A(0,p) @ B(p,q) @ C(r,1) @ D(q,r)
print("relative error in check  3:", check(raw(ABCD3), test))

# This illustrates that +, - and * all work on tensors networks, and the call notation "()" can
# also be used just to permute indices.  The result of a scalar multiplication (which is fast)
# is another tensor network whose internal scalar factor is modified.  Similarly, the result of the
# subtraction is an abstracted structure that keeps track of the operation but performs no 
# computation  until evaluate() or raw() is used.
ABCD_p = 3*ABCD1 - 3*ABCD1(1,0)
test = 3*test - 6*einsum("pr,rt,sq,ts->qp", A_, B_, C_, D_)
print("relative error in check  4:", check(raw(ABCD_p), test))

# This illustrates that there are not many restrictions on the order of the arguments or the 
# placement of the indices.  The redundant "1" means that the 1 index of the contracted product
# is resolved by setting these two equal to each other.  Also single as well as multiple (not just
# double) occurances of contraction labels are allowed.  Note that the same index can occur
# multiply on any given tensor, and that outer products (A shares no indices with others) are allowed.
ABCD3 = A(q,0) @ B(1,p) * 3 @ C(p,1) @ D(p,p)
test = 3 * einsum("qr,sp,ps,pp->rs", A_, B_, C_, D_)
print("relative error in check  5:", check(raw(ABCD3), test))

# The real utility of tensornet is for something like this.  E represents 4-index tensor, but
# its contractions with vectors can be done faster if done internally as contractions with the
# factors.  The scalar_value() function extracts the scalar-typed value from a 0-index tensor
# (like ndarray.item() in numpy).
E = A(0,1) @ D(2,3)
VVEUU = scalar_value(E(p,q,r,s) @ V(p) @ V(q) @ U(r) @ U(s))
test = einsum("pq,rs,p,q,r,s->", A_, D_, V_, V_, U_, U_).item()
print("relative error in check  6:", (VVEUU - test) / test)

# You can take elements of a primitive tensor, ...
print("relative error in check  7:", (A[0,0] - A_[0,0]) / A_[0,0])
# ... or of a tensor network.  For a network, this will reevaluate everything for every element
# requested, so it warns the user the first time about this being inefficient (maybe use
# evaluate() on a block and then take elements instead?).  To be clear, tensornet is smart
# enough not evaluate the whole tensor to obtain a single element; it first restricts the
# relevant indices in the primitives before contraction is done.  So, in this particular case,
# it actually is not inefficent.
print("relative error in check  8:", (E[0,0,0,0] - A_[0,0]*D_[0,0]) / (A_[0,0]*D_[0,0]))
# Note that it does not warn you a second time.
print("relative error in check  9:", (E[1,1,1,1] - A_[1,1]*D_[1,1]) / (A_[1,1]*D_[1,1]))

# Slicing also works, ...
print("relative error in check 10:", check(raw(A[0, 1:5]), A_[0, 1:5]))
# ... within contractions ...
AC = A[:,:8](0,p) @ C(p,1)
test = einsum("pq,qr->pr", A_[:,:8], C_)
print("relative error in check 11:", check(raw(AC), test))
# ... and on tensor networks (where again it first acts to slice the primitives in the network
# before evaluation is done.
test = einsum("p,q->pq", A_[:,0], D_[:,1])
print("relative error in check 12:", check(raw(E[:,0,:,1]), test))



    ### Some Finer Points

# One thing that is good to be aware of is that, if one contracts a sum (such as F below),
# the contractions of the terms are always done first and then added (most likely desired
# and far simpler dispatching algorithm) ...
F = A + 5*B(0,p)@C(p,1)
test = einsum("p,pq,q->", V_, A_, V_) +  5 * einsum("p,pq,qr,r->", V_, B_, C_, V_)
print("relative error in check 13:", check(raw(V(p) @ F(p,q) @ V(q)), test))
# ... but the user has control over this by forcing evaluation of the intermediate first,
# if desired.  The evaluate() function performs all internal contractions and index reductions
# and returns a primitive (not a raw) tensor.
print("relative error in check 14:", check(raw(V(p) @ evaluate(F)(p,q) @ V(q)), test))

# Should one need to loop over terms to build a sum, an empty accumulator can be started as so:
F = tensor_sum()
for term in [A, 5*B(0,p)@C(p,1)]:    # pretend these would be generated by some algorithm:
    F += term
print("relative error in check 15:", check(raw(V(p) @ F(p,q) @ V(q)), test))

# Here is another way to do this, if you know you want a concrete evaluation immediately.
# This does not yet save memory with an in-place build, but the point is to leave the
# door open for that to be implemented later (presently just a thin wrapper around
# +=raw(...)).
F_ = zeros(A.shape)
for term in [A, 5*B(0,p)@C(p,1)]:
    increment(F_, term)
print("relative error in check 16:", check(einsum("pq,p,q->", F_, V_, V_), test)), 

# And now a comment on a subtlty of the @ operator.  By necessity (since there is no 
# trigger to tell it to stop), it generates "incomplete" contractions.  Consider "ABCD" from
# above.  Done in the way below, the contraction labels q and r must be shared across two
# lines.  This is a potential danger since one might not realize a letter has been shared if
# lines like these are farther apart.  Whereas there might actually be a niche use case for
# this behavior it is mostly just a dark corner.
ABC4  = A(0,p) @ B(p,q) @ C(r,1)
ABCD4 = 2 * ABC4 @ D(q,r)
test = 2 * einsum("pr,rt,sq,ts->pq", A_, B_, C_, D_)
print("relative error in check 17:", check(raw(ABCD4), test))
# However, under normal circumstances, one would not write code in this way, so the danger
# is limited.  If one really means that, in ABC4, q and r should be summed over as individual
# indices, leaveing only two free indices in the product, then future contractions should
# be written more like the following (note the call with indices attached to ABC4), where we
# have to choose something other than D just for dimension matching.
ABCB1 = 2 * ABC4(p,0) @ B(p,1)
test = 2 * einsum("pq,qr,st,pu->tu", A_, B_, C_, B_)
print("relative error in check 18:", check(raw(ABCB1), test))
# And as a one-liner, we could have
ABCB2 = 2 * (A(0,p) @ B(p,q) @ C(r,1))(p,0) @ B(p,1)
print("relative error in check 19:", check(raw(ABCB2), test))
# And, just for clarity, this is of course the same (in every way) as
ABCB3 = 2 * A(s,p) @ B(p,q) @ C(r,0) @ B(s,1)
print("relative error in check 20:", check(raw(ABCB3), test))



    ### About Efficiency

# Here is a quick test that shows this code can be smarter than einsum in dispataching contractions.
# Setup ...
factor = 1
if len(sys.argv)==2:
    factor = float(sys.argv[1])
dim = math.floor(factor * 10)
M_ = random_tensor(dim, dim, dim, dim)
M  = prim_tensor(M_)
T_ = random_tensor(dim, dim)
T  = prim_tensor(T_)
# ... and test:
TTMTT = raw(M(p,q,r,s) @ T(p,0) @ T(q,1) @ T(r,2) @ T(s,3))
print("Tensornet done with 4-index transformation.  Waiting on einsum ... ")
test = einsum("pqrs,pw,qx,ry,sz->wxyz", M_, T_, T_, T_, T_)    # See note about opt_einsum in to-do (default in pytorch and accessible with numpy, also via tensorly)
print("... einsum done.")
print("relative error in check 21:", check(TTMTT, test))
print("If the timing difference was not dramatic enough, run the script with a float argument >~1")
print("to scale the dimensions.  But do not go nuts.  The einsum call scales with the 8th power.")



    ### About Backends

# Finally, the concept of a backend is local to the tensors.  Rather than have a global setting,
# any two tensors with the same backend may enter into tensornet expressions with each other.
# The long way of doing this for a numpy backend would be as below.  The default is not to make
# a copy of the original tensor data, but simply to wrap it.  As long as the user has no intention
# of changing the original data (tensornet will not either), it is more efficient, but if you 
# are getting weird results and want to be sure, or you just know you cannot promise not to
# modify the original, then use the copy_data flag.  The np_tensor function used above, is
# equivalent to the first option below and is provided for convenience. 

# from qode.math.tensornet import primitive_tensor_wrapper, numpy_backend
# np_tensor_ = primitive_tensor_wrapper(numpy_backend)                    # this ...
# np_tensor_ = primitive_tensor_wrapper(numpy_backend, copy_data=True)    # ... or this
# A = np_tensor_(A_)

# One can see now that providing a substitute for numpy_backend above is all one needs for a custom
# backend.  This should be a module that has the following functions defined in it.
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
#
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
