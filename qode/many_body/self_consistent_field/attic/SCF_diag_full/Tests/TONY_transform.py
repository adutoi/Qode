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


def transform2idx(A,U):
    # Assume A is sqaure and U has same leading dimension, which is rectangular
    dim = len(A)
    dimU = len(U[0])
    B = []
    for i in dim:
        row = []
        B.append(row)
        for j in dimU:
            element = 0
            for k in dim:  element = element + A[i][k]*U[k][j]
            row.append(element)
    C = []
    for i in dimU:
        row = []
        C.append(row)
        for j in dimU:
            element = 0
            for k in dim:  element = element + B[k][j]*U[k][i]
            row.append(element)
    return C

