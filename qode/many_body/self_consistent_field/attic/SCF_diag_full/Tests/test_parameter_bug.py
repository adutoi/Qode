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


def build(size):
    return [[ 1+i+j*size for i in range(size)] for j in range(size)]

def test_func(A):
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] += 10
    return A


M = build(4)

for i in M:
    print(i)

B = test_func(M)

for j in M:
    print(j)
for k in B:
    print(k)
