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
import sys
from math import exp,erf,sqrt,pi

def F0(T):
    # For atoms with only s orbitals (H,He)
    if T <= 15.0:
        t = 1.0
        s = t*exp(-T)
        n = 1       
        while n<=100:
            t *= 2.0*T/(2.0*n+1.0)
            s += t*exp(-T)
            n += 1
        #prt += str(round(T,2))+'\t'
        #correct = erf(sqrt(T))/(2*sqrt(T/pi))
        #prt += str(fabs((s-correct)/correct))+'\n'
        return s
    else:
        print('T out of range!')
        sys.exit(1)






def dfact(m):
    "Double Factorial"
    result = 1
    if (2*m-1) < 2:
        result = 1
    else:
        for i in range(3,2*m,2):
            result *= i
    return result




        

def F(m,T):
    # For all atoms 
    if T <= 15.0:
        c = dfact(m)
        tm = 2*m
        # Iteration Part:
        t = 1.0/dfact(tm + 1)
        print("m,tm,t",m,tm,t)
        s = t
        n = 1
        while n<=100:
            t *= 2.0*T/(tm + 2.0*n + 1.0)
            s += t
            n += 1
        #prt += str(round(T,2))+'\t'
        #correct = erf(sqrt(T))/(2*sqrt(T/pi))
        #prt += str(fabs((s-correct)/correct))+'\n'
        s *= c*exp(-T)  
        return s




maxdiff = -1.0


m = 4
T = 0.1
while T <= 15.0:
    a = F(m,T)
    b = erf(sqrt(T))/(2*sqrt(T/pi))
    print(a,b,abs(a-b)/b)
    if abs(a-b)> maxdiff:
        maxdiff = abs(a-b)
    T += 1.0
print(maxdiff)















