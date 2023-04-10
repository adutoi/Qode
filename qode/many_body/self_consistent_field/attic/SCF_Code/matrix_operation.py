#    (C) Copyright 2016 Yuhong Liu
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
import numpy as np

def cp_arbit_np_mat(matrix):
    'DEEP copy arbitrary numpy matrix'
    I = [[x,y] for x,y in np.ndindex(matrix.shape)]
    lda = len(matrix)
    wid = len(I)//lda
    new_mat = np.matrix(np.zeros((lda,wid)))
    for i,j in I:
        new_mat[i,j] = matrix[i,j]
    return new_mat



def mat_replicate(M,N):
    'M, N: same size matrices'
    x = 0
    for i in M:
        y = 0
        for j in i:
            M[x][y] = N[x][y]
            y += 1
        x += 1
    return M

def mat_np_rep(M,N):
    'M, N: same size numpy matrices'
    size = len(M)
    for i in range(size):
        for j in range(size):
            M[i,j] = N[i,j]
    return M

def build_mat(input_list):
    size = len(input_list)
    return [[0.0 for i in range(size)] for j in range(size)]

def build_mat_rep(input_list):
    size = pow(len(input_list),2)
    return [[0.0 for i in range(size)] for j in range(size)]
    
def build_size_mat(size):
    return [[0.0 for i in range(size)] for j in range(size)]


def mat_add(M1,M2):
    'Two Python matrices, return addition'
    n = len(M1)
    M = build_size_mat(n)
    for i in range(len(M1)):
        for j in range(len(M1[i])):
            M[i][j] = M1[i][j] + M2[i][j]
    return M


def mat_multiply(M1,M2):
    'Two matrices, M1:contraction coefficents, M2:integrals, multiply accordingly'
    n = len(M1)
    M = build_size_mat(n)
    for i in range(len(M1)):
        for j in range(len(M1[i])):
            M[i][j] = M1[i][j] * M2[i][j]
    return M

def printline(matrix):
    for line in matrix:
        print(line)

def printnpline(matrix):
    I = [[x,y] for x,y in np.ndindex(matrix.shape)]
    output=''
    t = 0
    for i,j in I:
        if t < i:
            t = i
            output += '\n'
        if t==i:
            output += "%e   " %(matrix[i,j])
    print(output)

def writenpline(matrix,filename):
    I = [[x,y] for x,y in np.ndindex(matrix.shape)]
    output=''
    t = 0
    for i,j in I:
        if t < i:
            t = i
            output += '\n'
        if t==i:
            output += str(matrix[i,j])+'  '
    f = open(filename,'w')
    f.write(output)
    f.close()
    

def build_same_py_mat(py_matrix):
    'Input:Matrix; Output:same matrix at different location'
    # square python matrix only
    size = len(py_matrix)
    A = build_size_mat(size)
    A = mat_replicate(A,py_matrix)
    return A
    
    
    
def build_same_np_mat(np_matrix):
    'Input:Matrix; Output:same matrix at different location'
    # square numpy matrix only
    size = len(np_matrix)
    A = build_size_mat(size)
    B = np.matrix(A)
    B = mat_np_rep(B,np_matrix)
    return B
    
def empty_same_np_mat(np_matrix):
    'Input:Matrix; Output:same-size empty matrix at different location'
    A = build_same_np_mat(np_matrix)
    B = A - A
    return B
