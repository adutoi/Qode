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
'''
M = [[1,2],
     [3,4]]
N = [[10,20],
     [30,40]]

def mat_add(M1,M2):
    'Too matrices, return addition'
    M = M1
    for i in range(len(M1)):
        for j in range(len(M1[i])):
            M[i][j] = M1[i][j]+M2[i][j]
    return M


print(mat_add(M,N))
#print(len(M))
'''
def centerf(p,q,x,y):
    return (p*x+q*y)/(p+q)

def center_coor(l1,l2):
    p = l1[3]
    q = l2[3]
    x = centerf(p,q,l1[0],l2[0])
    y = centerf(p,q,l1[1],l2[1])
    z = centerf(p,q,l1[2],l2[2])
    return [x,y,z,0.0,0.0]

l1=[0.0,0.0,0.0,2.0,1.2]
l2=[1.0,2.0,3.0,8.0,1.2]
print(center_coor(l1,l2))
