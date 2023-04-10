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
import static

lines = open("hamiltonian_lines.txt","r").readlines()

buffer = []

for line in lines:
	if line.rstrip():
		buffer += [line]
	else:
		filename = ""
		i = 0
		n = len(buffer[-1])
		while i<n:
			if buffer[-1][i:i+3]=="Occ":  filename += "Occ"
			if buffer[-1][i:i+3]=="Vrt":  filename += "Vrt"
			i += 1
		file = open("output/"+filename+".py","w")
		file.write(static.top)
		file.write("restrict_blocks(False)\n")
		for ln in buffer:  file.write(ln)
		file.write(static.bottom.format(filename+".tex"))
		file.close()
		buffer = []

print("Don\'t forget to add \"expression = add(expression)\" immediately before writing in OccVrt and VrtOcc!")
