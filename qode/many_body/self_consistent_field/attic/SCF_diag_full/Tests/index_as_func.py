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

def index(i,j,number):
    if number == 1:
        x = i
        y = j
    elif number == 2:
        x = j
        y = i
    return x,y


def loop(value):
    size = 4
    A = [[ 1+i+j*size for i in range(size)] for j in range(size)]
    A = np.matrix(A)
    print(A)
    for i in range(size):
        for j in range(size):
            print(A[index(i,j,value)])

loop(1)
loop(2)
