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



def readall(file):
    f = open(file,'r')
    cont = f.readlines()
    f.close()
    return cont

def chk_dif(x,y):
    x = float(x)
    y = float(y)
    return abs(x-y)

def main():
    a = readall('V_Direct.txt')
    b = readall('V_Success.txt')

    for i in range(len(a)):
        diff = chk_dif(a[i],b[i])
        if diff >= 1e-6:
            print(i,diff,a[i],b[i])
    print('Checking finished!')

        
if __name__=='__main__':
    main()
