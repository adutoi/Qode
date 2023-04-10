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
from math import sqrt


# Bohr Radius: 0.52917721092 angstrom
# Intake nuclei list, nuclei charge
# Output Nuclear Repulsion Energy


def distance(t1,t2):   
    # INPUT: _list_ elements t1 and t2
    # OUTPUT: First three x,y,z coordinates to give distance squared
    R2 = pow(t1[0]-t2[0],2) + pow(t1[1]-t2[1],2) + pow(t1[2]-t2[2],2)
    return sqrt(R2)/0.52917721092



def nurepul_main(_nulist_,A_num_list):

    
    NRE = 0.0
    n = len(_nulist_)
    #print("DEBUG")
    #print("n",n)
    P = _nulist_
    C = A_num_list
    for i in range(n):
        for j in range(i+1,n):
            d = distance(P[i],P[j])
            E = C[i]*C[j]/d
            NRE += E
    
    #print("NRE:",NRE)

    return NRE
