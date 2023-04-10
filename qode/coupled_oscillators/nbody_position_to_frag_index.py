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
############  n-body arbitary index break down into fragment state index  ################

# INPUT:     x (or y) position in n-body Hamiltonian (or n-body coupling matrix, which may not be built explicitly
#                                                    ( x or y must be less than total dimension )
#            size_list:  number of states in each fragment make a list ( dimensions of each fragment )
#
# OUTPUT:    The breakdown individual index of each fragment for any valid given position 
#
#

def position_to_index( position, size_list ):
    num_frag = len(size_list)
    remainders = []
    for i in range(num_frag):
        if i < num_frag - 1:
            remainders += [ position % size_list[-1-i] ]
            position = position // size_list[-1-i]
        elif i == num_frag - 1:
            remainders += [ position ] 
    remainders.reverse()
    return remainders


#
# This is a test, size_list is arbitary
#
if __name__ == "__main__":
    size_list = [ 1,2,3,4,5 ]
    mul = 1
    for i in size_list:
        mul *= i
    print("size =" ,mul)

    for i in range(mul):
        index = position_to_index( i, size_list )
        print("position(",i,")=" ,index)



