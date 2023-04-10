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
def index(M,i,j,mat_type):
    if mat_type == 1:
        M[i][j] = -M[i][j]
    else:
        M[i,j] = -M[i,j]


def rev_sign(M):
    #print(M)
    #print('new M:')
    n = len(M)
    "mat_type = 1: python matrix"
    "mat_type = 2: numpy matrix"
    '''try:
        c = M[0,0]
        mat_type = 2
        print('Type 2')
    except:
        mat_type = 1
        print('Type 1')
    '''
    for i in range(n):
        for j in range(n):
            index(M,i,j,1)

    print(M)
            

