#    (C) Copyright 2023 Anthony D. Dutoi
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
import numpy
from . import tensornet

def _SVD(M):
    return numpy.linalg.svd(M)

def svd_decomposition(nparray_Nd, indices_A, indices_B=None, thresh=1e-6, wrapper=tensornet.np_tensor):
    if indices_B is None:  indices_B = []
    all_free_indices = list(indices_A) + list(indices_B)
    if list(sorted(all_free_indices))!=list(range(len(nparray_Nd.shape))):
        raise RuntimeError("dimension mismatch")
    if len(indices_A)==0 or len(indices_B)==0:
        return wrapper(nparray_Nd)
    #
    nparray_Nd = nparray_Nd.transpose(all_free_indices)
    shape_A = nparray_Nd.shape[:len(indices_A)]
    shape_B = nparray_Nd.shape[len(indices_A):]
    M = nparray_Nd.reshape((numpy.prod(shape_A), numpy.prod(shape_B)))
    #
    U, s, Vh = _SVD(M)
    d, D = 0, len(s)
    for i in range(D):
        if abs(s[i])>thresh:  d += 1
    print("compression factor =", d/D)
    Shalf = numpy.sqrt(s[:d])
    A = (   U[:,:d] * Shalf).T
    B = (Vh.T[:,:d] * Shalf).T
    #
    A = wrapper(A.reshape([d] + list(shape_A)))
    contract_A = ["i"] + list(indices_A)
    B = wrapper(B.reshape([d] + list(shape_B)))
    contract_B = ["i"] + list(indices_B)
    return A(*contract_A) @ B(*contract_B)
