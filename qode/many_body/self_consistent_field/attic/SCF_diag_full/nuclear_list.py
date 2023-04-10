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
def readline(file):
    f = open(file,'r')
    lines = f.readlines()
    f.close()
    return lines


################### MAIN FUNCTION  ###################################
# RETURN A LIST OF NUCLEUS COORDINATES ONLY FOR THE NUCLEAR ATTRATION
# MATRIX.

def numain(inputfile):
    nulist = []
    tmp = readline(inputfile)
    for i in tmp:
        line = i.split()
        if len(line) == 4:
            nulist += [[float(line[1]),float(line[2]),float(line[3])]]
    return nulist
    
    
