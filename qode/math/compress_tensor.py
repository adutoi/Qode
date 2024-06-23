#    (C) Copyright 2023, 2024 Anthony D. Dutoi
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
from . import tensornet, space, numpy_space, vector_set, gram_schmidt

def _SVD(M, big, print_info):
    A = M
    T_left, T_right = None, None
    n_left, n_right = M.shape
    if n_left>n_right and n_left>big:
        S = space.linear_inner_product_space(numpy_space.real_traits(n_left))
        T_left = numpy.array(M.T)                                 # copy M -> T
        vecs = vector_set(S, [S.member(T_left[i,:]) for i in range(n_right)])    # work with rows of T in place
        gram_schmidt.gram_schmidt(vecs, normalize=True, n_times=n_right)
        if print_info:
            print("ON deviation = ", vecs.orthonormality())
        M = T_left @ M
    if n_right>n_left and n_right>big:
        S = space.linear_inner_product_space(numpy_space.real_traits(n_right))
        T_right = numpy.array(M.T)                                # copy M -> T
        vecs = vector_set(S, [S.member(T_right[:,i]) for i in range(n_left)])    # work with columns of T in place
        gram_schmidt.gram_schmidt(vecs, normalize=True, n_times=n_left)
        if print_info:
            print("ON deviation = ", vecs.orthonormality())
        M = M @ T_right
    U, s, Vh = numpy.linalg.svd(M)

    if print_info:
        Z = U[:,:len(s)] @ numpy.diag(s) @ Vh[:len(s),:]
        print("(0) Error =", numpy.linalg.norm(M), numpy.linalg.norm(M-Z))

    if T_left is not None:
        U  = T_left.T @ U
    if T_right is not None:
        Vh = Vh @ T_right.T

    if print_info:
        Z = U[:,:len(s)] @ numpy.diag(s) @ Vh[:len(s),:]
        print("(1) Error =", numpy.linalg.norm(A), numpy.linalg.norm(A-Z))

    return U, s, Vh

def svd_decomposition(nparray_Nd, indices_A, indices_B=None, thresh=1e-6, big=10**5, wrapper=tensornet.np_tensor, print_info=False):
    if print_info:
        print("size =", numpy.prod(nparray_Nd.shape))
    if indices_B is None:  indices_B = []
    all_free_indices = list(indices_A) + list(indices_B)
    if list(sorted(all_free_indices))!=list(range(len(nparray_Nd.shape))):
        raise RuntimeError("dimension mismatch")
    if len(indices_A)==0 or len(indices_B)==0:
        return wrapper(nparray_Nd)
    #
    Y = nparray_Nd
    nparray_Nd = nparray_Nd.transpose(all_free_indices)
    shape_A = nparray_Nd.shape[:len(indices_A)]
    shape_B = nparray_Nd.shape[len(indices_A):]
    M = nparray_Nd.reshape((numpy.prod(shape_A), numpy.prod(shape_B)))
    #
    U, s, Vh = _SVD(M, big, print_info)
    d, D = 0, len(s)
    for i in range(D):
        if print_info:
            print("{:.1e} ".format(s[i]), end="")
        if abs(s[i])>thresh:  d += 1
    if print_info:
        print("\ncompression factor =", d/D)

    if d==0:
        d = 1    # a bit of a hack, but otherwise get some crazy errors

    Shalf = numpy.sqrt(s[:d])
    A = (   U[:,:d] * Shalf).T
    B = (Vh.T[:,:d] * Shalf).T
    #
    A = wrapper(A.reshape([d] + list(shape_A)))
    contract_A = ["i"] + list(indices_A)
    B = wrapper(B.reshape([d] + list(shape_B)))
    contract_B = ["i"] + list(indices_B)

    Z = tensornet.raw(A(*contract_A) @ B(*contract_B))
    if print_info:
        print("(2) Error =", numpy.linalg.norm(Y), numpy.linalg.norm(Y-Z))

    return A(*contract_A) @ B(*contract_B)
