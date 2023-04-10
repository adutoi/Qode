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

def textgen(m,n,o,p,i,j,k,l):
    v = str(m)+str(n)+str(o)+str(p)
    u1='U'+str(i)+str(m)
    u2='U'+str(j)+str(n)
    u3='U'+str(k)+str(o)
    u4='U'+str(l)+str(p)
    return 'V'+v+u1+u2+u3+u4

def loop(i,j,k,l):
    size = 2

    text = ''
    for m in range(size):
        for n in range(size):
            for o in range(size):
                for p in range(size):
                    text += textgen(m,n,o,p,i,j,k,l)
                    text += '  '
            text += '\n'
    text += '\n'
    return text

def main():
    size = 2
    output = ''
    for i in range(size):
        for j in range(size):
            for k in range(size):
                for l in range(size):
                    output += loop(i,j,k,l)

                    
    f = open('test_ijkl.txt','w')
    f.write(output)
    f.close()
    print('Job Done!')

main()
