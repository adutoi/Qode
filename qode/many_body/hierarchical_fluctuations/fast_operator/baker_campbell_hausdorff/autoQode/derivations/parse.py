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

lines = open("human-readable.txt","r").readlines()

spacings = []
for n,c in enumerate(lines[0]):
	if c=='.':  spacings += [n]


for line in lines[1:]:
	fields = []
	was = 0
	if line.rstrip():
		for spacing in spacings:
			fields += [ line[was:spacing] ]
			was = spacing
		fields += [ line[was:] ]
		fields = [field.rstrip() for field in fields]
		#
		occs = fields[-2].split()
		occs = [occ[0] for occ in occs]
		vrts = fields[-1].split()
		vrts = [vrt[0] for vrt in vrts]
		#
		i = 0
		n = len(fields[0])
		abbrev1 = ""
		while i<n-1:
			if fields[0][i:i+2] in ["Ex", "Fo", "Fv", "Dx"]:  abbrev1 += fields[0][i:i+2]
			i += 1
		abbrev1 = abbrev1.rjust(6,"_")
		i = 0
		n = len(fields[1])
		abbrev2 = ""
		while i<n-1:
			if fields[1][i:i+2] in ["Ex", "Fo", "Fv", "Dx"]:  abbrev2 += fields[1][i:i+2]
			i += 1
		filename = "output/commute_" +abbrev1 + "_" +abbrev2 + ".py"
		outfile = "commute_" +abbrev1 + "_" +abbrev2 + ".tex"
		output = open(filename, "w")
		output.write(static.top)
		#
		summand1 = "mult(" + fields[0] +  ")"
		summand2 = "mult(" + fields[1] +  ")"
		summand  = "comm( " + summand1 + " , " + summand2 + " )"
		p = fields[2].index(")(")
		fields[2] = fields[2][:p] + "),(" + fields[2][p+2:]
		summand1 = "Coeff(\'" + fields[2][0] + "\'," + fields[2][1:] + ")"
		p = fields[3].index(")(")
		fields[3] = fields[3][:p] + "),(" + fields[3][p+2:]
		summand2 = "Coeff(\'" + fields[3][0] + "\'," + fields[3][1:] + ")"
		summand = "mult(" + summand1 + ", " + summand2 + ", " + summand + ")"
		output.write("IDX_{} = Indices(\"{}\", bound=\"No\")\n".format("".join(occs), fields[-2]))
		output.write("{}, = IDX_{}()\n".format(",".join(occs), "".join(occs)))		# trailing comma needed for 1-tuples
		output.write("\n")
		output.write("IDX_{} = Indices(\"{}\", bound=\"Nv\")\n".format("".join(vrts), fields[-1]))	# trailing comma needed for 1-tuples
		output.write("{}, = IDX_{}()\n".format(",".join(vrts), "".join(vrts)))
		output.write("\nexpression = Sum(IDX_{},IDX_{})({})\n".format("".join(occs),"".join(vrts),summand))
		output.write(static.bottom.format(outfile))
		output.close()
