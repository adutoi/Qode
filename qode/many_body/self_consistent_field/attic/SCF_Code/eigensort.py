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
#IMPUTE: Numpy 1-D Array E
#OUTPUT: eigenvalue order and corresponding indices.
import numpy as np


def findN(M):
    #Return number of elements.
    ct = 0
    for x in np.nditer(M):
        ct += 1
    return ct

def sort(I,N,n):
    'I:Index, N:Numbers, n:dimension'
    t = 0.0
    indt = 0
    for i in range(n-1):
        for j in range(n-1-i):
            if N[j] > N[j+1]:
                t = N[j+1]
                N[j+1] = N[j]
                N[j] = t
                indt = I[j+1]
                I[j+1] = I[j]
                I[j] = indt
    
    return I,N

def sort_main(E):
    n = findN(E)//2
    Na = np.array([E[i] for i in range(n)])
    Ia = np.array([i for i in range(n)])
    Nb = np.array([E[i+n] for i in range(n)])
    Ib = np.array([i+n for i in range(n)])
    
    print("UNSORTED Guess Vector Energies")
    print(Ia,Na)
    print(Ib,Nb)
    
    Ia,Na = sort(Ia,Na,n)
    Ib,Nb = sort(Ib,Nb,n)
    
    print("SORTED   Guess Vector Energies ")
    print(Ia,Na)
    print(Ib,Nb)
    
    return Ia,Ib








