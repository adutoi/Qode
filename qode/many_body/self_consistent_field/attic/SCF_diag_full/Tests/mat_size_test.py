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
def build_empty(size):
    return [[ 0.0 for i in range(size)] for j in range(size)]

def printline(matrix):
    for i in matrix:
        print(i)

def loop_n(M,ng):
    N = build_empty(4)
    #The Loop:
    acc = []
    x = 0
    for m in range(len(ng)):
        for i in range(ng[m]):
            y = 0
            for n in range(len(ng)):
                for j in range(ng[n]):
                    N[m][n] += M[x][y]
                    acc += [[x,y,m,n]]
                    y += 1
            x += 1
    f=open('check_small.txt','w')
    to_write = ''
    for line in acc:
        to_write += str(line)+'\n'
    f.write(to_write)
    f.close()        
    return N
    

def loop_r(M,ng):
    #N = build_mat a different size one
    #The Loop:
    N = build_empty(16)
    size = len(ng)

    acc = []
    x = 0
    for p in range(len(ng)):
        for i in range(ng[p]):        
            for q in range(len(ng)):
                for j in range(ng[q]):
                    y = 0
                    for r in range(len(ng)):
                        for k in range(ng[r]):
                            for s in range(len(ng)):
                                for l in range(ng[s]):
                                    N[p*size + q][r*size + s] += M[x][y]
                                    acc+= [[x,y,p*size + q,r*size + s]]
                                    y += 1
                    x += 1
    f=open('check_cont.txt','w')
    to_write = ''
    for line in acc:
        to_write += str(line)+'\n'
    f.write(to_write)
    f.close()
    
    return N

ng = [3,1,3,1]
A = build(8)
B = build(64)

printline(A)
C = loop_n(A,ng)
printline(C)

#printline(B)
D = loop_r(B,ng)
printline(D)



