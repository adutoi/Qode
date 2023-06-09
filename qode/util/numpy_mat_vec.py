#    (C) Copyright 2018 Anthony D. Dutoi
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
import random
import numpy

def random_matrix(dim, seed=12011976):
    random.seed(seed)
    mat = numpy.array([[0.]*dim]*dim)
    for i in range(dim):
        for j in range(i+1,dim):  mat[j,i] = mat[i,j] = random.uniform(-0.1,0.1)
        mat[i,i] = i+1 + random.uniform(-0.1,0.1)
    return mat

def zero_vector(dim):
    return numpy.array([0.]*dim)

def basis_vector(n,dim):
    vec = zero_vector(dim)
    vec[n] = 1.
    return vec
