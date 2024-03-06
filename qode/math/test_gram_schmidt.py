#    (C) Copyright 2024 Anthony D. Dutoi
# 
#    This file is part of QodeApplications.
# 
#    QodeApplications is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    QodeApplications is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with QodeApplications.  If not, see <http://www.gnu.org/licenses/>.
#

# If this is working, then the second block printed should be all near machine zero.

import numpy
from qode.math import linear_inner_product_space, numpy_space, gram_schmidt



def random_vec(dim):
    return 2 * numpy.random.random((dim,)) - 1    # all elements between -1 and 1

def print_overlap_dev(vecs):
    for i,bra in enumerate(vecs):
        for j,ket in enumerate(vecs):
            overlap = (bra|ket)
            if i==j:  overlap -= 1
            print(" {: 9.1e}".format(overlap), end="")
        print()
    print()



num = 10
dim = 10000

S = linear_inner_product_space(numpy_space.real_traits(dim))
vecs = [S.member(random_vec(dim)) for _ in range(num)]

print_overlap_dev(vecs)
gram_schmidt(vecs, normalize=True, n_times=2)
print_overlap_dev(vecs)


