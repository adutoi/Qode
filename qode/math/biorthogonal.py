#    (C) Copyright 2019 Anthony D. Dutoi
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

def precise_numpy_inverse(M):
        M_inv  = numpy.linalg.inv(M)
        refine = numpy.linalg.inv(M_inv @ M)    # inverting something very close the the identity matrix
        return refine @ M_inv

# Let {|v^j>} be the full set of BO complements to the set {|v_j>}.  The resolution of a vector |u> in terms of {|v_j>}
# is the sum of "projections", |v_j><v^j|u>, into the 1-D spaces spanned by each of the |v_j>.  However, if what
# we want it the "projection" of |u> into the orthogonal complement of the 1-D space spanned by a given |v_j>,
# then we need
#                 |u'> = ( 1 - |v^j><v_j| ) |u>
# whereby it is easily verified that <v_j|u'> = 0.  The algorithm here works by orthogonalizing |u> to each vector in list
# V by sequentially projecting into the complement space of each.  Since the BO complements have overlap only with
# their associated vector of the original space, we have (for distinct ..., j, k, ...)
#      |u'>    =   ... ( 1 - |v^j><v_j| )( 1 - |v^k><v_k| ) ... |u>  =   ( 1 ... - |v^j><v_j| - |v^k><v_k| - ... ) |u>
# where it is easily verified that ... = <v_j|u'> = <v_k|u'> = ... = 0.
#
# The function below does the above orthogonalization for each vector |u>, in U to all of the vectors in V (with complements in Vc)
# excepting that vector in V which is supposed to be the complement of |u>, as indicated by the index that it comes in paired with.
# |u> is then normalized by projecting onto its purported complement.
#
# Importantly, U and V need not have the same size, one can orthogonalize the vectors in U to a larger set.  Also, it is permissible
# for U and Vc to be the same set or for U to be a subset of Vc (though then the Vc will changed, presumably slightly) over the course
# of this call.
#
# The call returns the matrix of overlaps, which will be a (non-square) block of the Id matrix at convergence
#
def biorthog_iteration(U, V, Vc=None):
    if Vc is None:  Vc = V	# assume ON set
    if len(Vc)!=len(V):  raise AssertionError
    overlaps = numpy.zeros((len(U),len(V)))
    for i,(k,u) in enumerate(U):
        for j,(v,vc) in enumerate(zip(V,Vc)):
            if j!=k:
                tmp = vc*(v|u)	# until u -= vc*(v|u) is allowed
                u -= tmp
        u /= (u|V[k])
        for j,v in enumerate(V):  overlaps[i,j] = (u|v)
    return overlaps

def iterative_biorthog(V, in_place=None, thresh=1e-12):      # V should be an iterable containing members of a space
    if in_place is not None:  U = in_place
    else:                     U = [v+0 for v in V]	# v+0 makes a copy
    for i,u in enumerate(U):  u /= (u|V[i])
    diff = float("inf")
    while diff>thresh:
        overlaps = biorthog_iteration(list(enumerate(U)), V, Vc=U)
        diff = numpy.linalg.norm( overlaps - numpy.identity(len(V)) ) / len(V)
        print(diff)
    return U
