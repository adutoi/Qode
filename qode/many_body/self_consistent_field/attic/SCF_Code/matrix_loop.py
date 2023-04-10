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
def mat_loop(input_list,operator):
    # INPUT: _list_: (x, y, z, exponential, contraction coefficient)
    # OUTPUT: Unitary Matrix
    # Normalized Gaussian Functions.    
    # BUILD a Matrix to hold data:
    # _list_ indicies: i, j.
    # S indicies: x,y
    #
    # x = 0
    # Loop i:
    #   y = 0
    #   Loop j:
    #      Compute each eqn: M[x][y] = ***
    #      y += 1
    #   x += 1
    M = []  
    for i in input_list:
        row = []
        for j in input_list:
            row += [operator(i,j)]
        M += [row]
    return M


def rep_loop(input_list,operator):
    M =[]
    for p in input_list:
        for q in input_list:
            row = []
            for r in input_list:
                for s in input_list:
                    row += [operator(p,q,r,s)]
            M += [row]                   
    return M
