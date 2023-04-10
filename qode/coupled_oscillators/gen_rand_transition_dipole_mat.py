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
# This module is for TEST PURPOSE
# It generates a NUMPY Transition Dipole Matrix in such form:
#    <i'|<j'| V |i >|j > = <i'|<j'| xi xj |i >|j >
#
# Non-Zero Condition:  i' = i +/- 1  AND j' = j +/- 1
#
# The multiplier used is   1.0/ Leading Dimension    to scale off-diagonal numbers down.
#



import numpy as np

def gen_two_body_x_matrix( ld1, ld2 ):
    lda = ld1 * ld2
    A = np.matrix( np.zeros((lda,lda)) )
    for i in range(ld1):
        for j in range(ld2):
            for k in range(ld1):
                for l in range(ld2):
                    if abs(i-k) == 1 and abs(j-l) == 1:
                        A[i*ld2 + j ,k*ld2 + l] = (i+j+k+l)*(1.0/lda)
    return A



if __name__ == "__main__":
    n1 = 5
    n2 = 6
    A = gen_two_body_x_matrix(n1,n2)
    print(A)
    print("show one non-zero element =",A[0,7])