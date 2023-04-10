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

A = [[1,10],
     [50,100]]
B = [[1,0],
     [0,1]]

def multiply(A,B):
    'A:i--k  B:k--j'
    leni = len(A)
    lenj = len(B[0])
    lenk = len(B)

    C = [[0.0 for m in range(leni)]for n in range(lenj)]
    
    for i in range(leni):
        for j in range(lenj):
            for k in range(lenk):
                C[i][j] += A[i][k]*B[k][j]

    return C


M = multiply(A,B)
for line in M:
    print(line)
