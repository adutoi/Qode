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


A = np.matrix([[1,2,3,4],
               [5,6,7,8],
               [9,10,11,12]])

print(len(A))
print(len(A[0]))
print(A[0,0])
for x in np.nditer(A):
    print(x)

for (x,y),value in np.ndenumerate(A):
    print(x,y)

I = [[x,y] for x,y in np.ndindex(A.shape)]
print(I)


output=''
t = 0
for i,j in I:
    if t < i:
        t = i
        output += '\n'
    if t==i:
        output += str(A[i,j])+'  '

print(output)

