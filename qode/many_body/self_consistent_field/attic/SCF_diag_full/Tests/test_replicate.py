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
from matrix_operation import *

def build(size):
    return [[ 1+i+j*size for i in range(size)] for j in range(size)]

X = build(4)
M = build_size_mat(4)
printline(X)
print(' ')
printline(M)

Y = build_same_py_mat(X)
print(' ')
printline(Y)

print(' ')
Z = mat_replicate(M,X)
printline(Z)
print(' ')
printline(M)
