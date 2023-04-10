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
"""\
Though the language of the comments betrays the fact that I was thinking about a nested_operator implementation when I wrote this,
the logic behind these truncations is more general, so this stuff stays here.
"""

from ...util.output import no_print
#from qode.util.output import no_print		# if run locally as script, prints human-readable parse of logic to the screen

def _is_boolean_int(n):
    # This replaces the expression 'n is 0 or n is 1' in if statements, which, as of python 3.8 began resulting in the following warning
    # SyntaxWarning: "is" with a literal. Did you mean "=="?
    # This is explained in the following article
    # https://adamj.eu/tech/2020/01/21/why-does-python-3-8-syntaxwarning-for-is-literal/
    # Summarized, some implementations of python do not guarantee that multiple instances of the integer 0 or 1 are all the same object.
    # This seems like an oversight in the language spec to me, but maybe its just really hard to enforce at the implementation level?
    return (isinstance(n,int) and (n==0 or n==1))

# Just a utility to move elsewhere
def bool_matmul(mat1,mat2):
	""" Takes two matrices of ones(true) and zeros(false) and performs matrix multiplication, interpreting * as AND and + as OR """
	mat = []
	for i in range(len(mat1)):
		if len(mat1[i])!=len(mat2):  raise Exception("mismatched contraction dimensions and/or non-rectangular first argument")
		row = []
		for j in range(len(mat2[0])):
			Aij = 0
			for k in range(len(mat2)):
				if len(mat2[k])!=len(mat2[0]):  raise Exception("second argument is not rectangular matrix")
				Mik = mat1[i][k]
				Nkj = mat2[k][j]
				if not _is_boolean_int(Mik):  raise Exception("first argument is not a boolean matrix")
				if not _is_boolean_int(Nkj):  raise Exception("second argument is not a boolean matrix")
				Aij += Mik * Nkj
			if Aij>0:  Aij = 1
			row += [Aij]
		mat += [row]
	return mat

def vec_to_column_matrix(vec):
	return [[v] for v in vec]

def column_matrix_to_vec(mat):
	return [r[0] for r in mat]

def mask(element):
	if _is_boolean_int(element[0]) and element[0]==0:  return ' '
	else:                                              return '@'

def and_or_arrays(arr1, arr2, op):
	try:
		if len(arr1)!=len(arr2):  raise Exception("arrays should have same dimensions")
	except:
		textlog("one argument is not array-like")
		raise
	dim = len(arr1)
	new_array = [0]*dim
	for i in range(dim):
		a1 = arr1[i]
		a2 = arr2[i]
		try:
			iter(a1)
		except:
			if not _is_boolean_int(a1):  raise Exception("first argument is not a boolean array")
			if not _is_boolean_int(a2):  raise Exception("second argument is not a boolean array")
			if op=="and":  a = a1 * a2
			if op=="or":   a = a1 + a2
			if op=="xor":  a = a1 - a2
			if a!=0:  new_array[i] = 1
		else:
			new_array[i] = and_or_arrays(a1, a2, op)
	return new_array
def and_arrays(arr1, arr2):  return and_or_arrays(arr1, arr2, op="and")
def or_arrays( arr1, arr2):  return and_or_arrays(arr1, arr2, op="or")
def xor_arrays(arr1, arr2):  return and_or_arrays(arr1, arr2, op="xor")

			


"""\
The matrix below, encoded in text form, was written by hand.  It is for a general 2-particle operator (Hamiltonian, H) commuting with an operator (T) that contains single and double excitations.
This could be generated algorithmically, which would be the prefered way of generalizing to CCSDT or to Hamiltonians with higher than 3-body action.  For this, the rules appended to this comment will
be necessary.

The matrix can be thought of as a kind of boolean scattering matrix ('.' and ' ' mean False(0) and '@' means True(1).  A column matrix (acted upon from the left by this matrix) is thought of as containing
boolean values as to the presense or absense of a constant term (K), single excitations (E), single-particle "flat" OO or VV rearrangements (F), single deexcitations (D), or subsequent ordered
strings of such operators, in an operator (X) which is to be commuted with T (X might be H or the result of a previous such commutation).  In this notation, E4FD=EEEEFD, for example.  
Only a subset of terms higher than 2-body are needed because only these will be possibly produced by nested commutation of H with T.
The action of the matrix (multiplication is AND, and addition is OR) then indicates what orders of terms will be present after the operator with structure given by the original
column matrix is commuted once with T.  Therefore, repeated application of the matrix shows what terms are generated by repeated (nested) commutation of H with T.

Here are the rules needed to generate this algorithmically.
It should be made clear that any many-body operator can, with the application of commutator rules be put in a form such that all operator strings have excitations (if present) followed
by flat operators of the same type, OO or VV (if present), followed by deexciations (if present), where all indices of each type (occ creation, virt creation, occ annihilation, virt annihilation)
appear in descending order.  The rules below rely on enforcement of such restrictions to understand what kinds of such terms might arise from commutations.
Finally, where as F is generic (either Fo or Fv), in the final result, intermediate stages will need to keep track.  For example, FoFv never occurs in the final result because it rearranges to ED, 
but it will occur in intermediates.  If a generic F occurs in multiple places in an equation, it is presumed to be either Fo or Fv in each instance.

	Upon rearrangement

EE = EE		FoFv = ED	FE = EF + E
EF = EF		FvFo = ED	FF = FF + F
ED = ED				DE = ED + Fo + Fv + K
FD = FD				DF = FD + D
DD = DD

	Base commutators

[E,E] = 0
[F,E] = E
[D,E] = Fo + Fv

	Recursive commutator rules (can be derived from [AB,C]=[A,C]B+A[B,C] and [A,BC]=[A,B]C+B[A,C])

[abc...z, E]     = [a,E]bc...z     + a[b,E]c...z     + ab[c,E]...z     + ... + abc...[z,E]
[X, EaEbEc...Ez] = [X,Ea]EbEc...Ez + Ea[X,Eb]Ec...Ez + EaEb[X,Ec]...Ez + ... + EaEbEc...[X,Ez]

The general strategy:  For a given type of Hamiltonian term, use commutation and rearrangement with *all* possible excitation levels, along with subsequent rearrangement rules, in order to find out
what kinds of terms result from the overall commutator.  This presence or absence of a certain result for a certain input develops the boolean column of a matrix like the one below.  In general,
after the first commutation, you will find that terms of the form of the original hamiltonian are present, along with terms of that form multipled from the left by *excitations only*.  Such terms
must be added to the basis, even though they are not in the original Hamiltonian, but resolving there commutatation products can be done by recursion, since you can treat only the "inner" part of such
components (already a solved problem), multipled from the left by the excitation order of the "outer" part.
"""



commute_scatter = """\
key:       .         .    .    .         .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .




           K         E    F    D        EE   EF   ED   FF   FD   DD      E1EE E1EF E1ED E1FF E1FD E1DD      E2EE E2EF E2ED E2FF E2FD E2DD      E3EE E3EF E3ED E3FF E3FD E3DD      E4EE E4EF E4ED E4FF E4FD E4DD




 K         .         .    .    @         .    .    .    .    .    @




 E         .         .    @    @         .    .    @    @    @    .         .    .    .    .    .    @

 F         .         .    .    @         .    .    .    .    @    @         .    .    .    .    .    .

 D         .         .    .    .         .    .    .    .    .    @         .    .    .    .    .    .




EE         .         .    @    .         .    @    @    @    .    .         .    .    @    @    @    .         .    .    .    .    .    @

EF         .         .    .    @         .    .    @    @    @    .         .    .    .    .    @    @         .    .    .    .    .    .

ED         .         .    .    .         .    .    .    .    @    @         .    .    .    .    .    @         .    .    .    .    .    .

FF         .         .    .    .         .    .    .    .    @    @         .    .    .    .    .    .         .    .    .    .    .    .

FD         .         .    .    .         .    .    .    .    .    @         .    .    .    .    .    .         .    .    .    .    .    .

DD         .         .    .    .         .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .




E1EE                                     .    @    .    .    .    .         .    @    @    @    .    .         .    .    @    @    @    .         .    .    .    .    .    @

E1EF                                     .    .    @    @    .    .         .    .    @    @    @    .         .    .    .    .    @    @         .    .    .    .    .    .

E1ED                                     .    .    .    .    @    .         .    .    .    .    @    @         .    .    .    .    .    @         .    .    .    .    .    .

E1FF                                     .    .    .    .    @    .         .    .    .    .    @    @         .    .    .    .    .    .         .    .    .    .    .    .

E1FD                                     .    .    .    .    .    @         .    .    .    .    .    @         .    .    .    .    .    .         .    .    .    .    .    .

E1DD                                     .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .




E2EE                                                                        .    @    .    .    .    .         .    @    @    @    .    .         .    .    @    @    @    .         .    .    .    .    .    @

E2EF                                                                        .    .    @    @    .    .         .    .    @    @    @    .         .    .    .    .    @    @         .    .    .    .    .    .

E2ED                                                                        .    .    .    .    @    .         .    .    .    .    @    @         .    .    .    .    .    @         .    .    .    .    .    .

E2FF                                                                        .    .    .    .    @    .         .    .    .    .    @    @         .    .    .    .    .    .         .    .    .    .    .    .

E2FD                                                                        .    .    .    .    .    @         .    .    .    .    .    @         .    .    .    .    .    .         .    .    .    .    .    .

E2DD                                                                        .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .




E3EE                                                                                                           .    @    .    .    .    .         .    @    @    @    .    .         .    .    @    @    @    .

E3EF                                                                                                           .    .    @    @    .    .         .    .    @    @    @    .         .    .    .    .    @    @

E3ED                                                                                                           .    .    .    .    @    .         .    .    .    .    @    @         .    .    .    .    .    @

E3FF                                                                                                           .    .    .    .    @    .         .    .    .    .    @    @         .    .    .    .    .    .

E3FD                                                                                                           .    .    .    .    .    @         .    .    .    .    .    @         .    .    .    .    .    .

E3DD                                                                                                           .    .    .    .    .    .         .    .    .    .    .    .         .    .    .    .    .    .




E4EE                                                                                                                                              .    @    .    .    .    .         .    @    @    @    .    .

E4EF                                                                                                                                              .    .    @    @    .    .         .    .    @    @    @    .

E4ED                                                                                                                                              .    .    .    .    @    .         .    .    .    .    @    @

E4FF                                                                                                                                              .    .    .    .    @    .         .    .    .    .    @    @

E4FD                                                                                                                                              .    .    .    .    .    @         .    .    .    .    .    @

E4DD                                                                                                                                              .    .    .    .    .    .         .    .    .    .    .    .


"""

class load_truncation_info(object):
	def __init__(self, textlog=no_print):

		# Load the matrix as list of text lines and getrid of "spacer" lines
		alllines = commute_scatter.split("\n")
		datalines = []
		for line in alllines:
			if line!="":  datalines += [line]
		# separate the key from the data
		keyline   = datalines[0]
		# line "1" contains titles
		datalines = datalines[2:]

		# The keyline tells us the locations where to expect a "number" (' ', '.', or '@')
		key = []
		for i in range(len(keyline)):
			if keyline[i]=='.':  key += [i]

		# Resolve the matrix as a dense fill of "boolean" values (0 and 1 to support straightforward multiplication(AND) and addition(OR))
		matrix = []
		for line in datalines:
			row = []
			paddedline = line + " "*1000	# pad with a lot of "zeros"
			for i in key:
				if paddedline[i]=='@':  row += [1]
				else:                   row += [0]
			matrix += [row]

		# The dimension of this thing is going to be pretty important
		dim = len(matrix)



		# For completeness, make an ID matrix of the same dimension
		id_mat = []
		for i in range(dim):
			row = [0 for j in range(dim)]
			row[i] = 1
			id_mat += [row]

		# Convenient map of labels to integers
		K, E, F, D, EE, EF, ED, FF, FD, DD, E1EE, E1EF, E1ED, E1FF, E1FD, E1DD, E2EE, E2EF, E2ED, E2FF, E2FD, E2DD, E3EE, E3EF, E3ED, E3FF, E3FD, E3DD, E4EE, E4EF, E4ED, E4FF, E4FD, E4DD = range(dim)

		# Convenient map of integers to labels
		labels = ["K   ", "E   ", "F   ", "D   ", "EE  ", "EF  ", "ED  ", "FF  ", "FD  ", "DD  ", "E1EE", "E1EF", "E1ED", "E1FF", "E1FD", "E1DD", "E2EE", "E2EF", "E2ED", "E2FF", "E2FD", "E2DD", "E3EE", "E3EF", "E3ED", "E3FF", "E3FD", "E3DD", "E4EE", "E4EF", "E4ED", "E4FF", "E4FD", "E4DD"]
		if len(labels)!=dim:  raise Exception("dimension mismatch!")

		# Convenient for instantiating and up-converting vectors of booleans and printing the result of matrix multiplies on them
		def vector():  return [0]*dim



		# Structure of the 2e- Hamiltonian
		H = vector()
		for block in [K, E, F, D, EE, EF, ED, FF, FD, DD]:  H[block] = 1

		textlog()
		textlog("Structure of the operator produced after 0, 1, 2, 3, 4, and 5 nested commutations if no truncations are ever made")
		textlog()
		textlog("          ", 0, 1, 2, 3, 4, 5, sep="   ")
		textlog()
		H0 = vec_to_column_matrix(H)
		H1 = bool_matmul(matrix,H0)
		H2 = bool_matmul(matrix,H1)
		H3 = bool_matmul(matrix,H2)
		H4 = bool_matmul(matrix,H3)
		H5 = bool_matmul(matrix,H4)
		for i in range(dim):
			textlog(labels[i], "   ",
			mask(H0[i]), mask(H1[i]), mask(H2[i]), mask(H3[i]), mask(H4[i]), mask(H5[i]),
			sep="   ")
		textlog()



		# Action of the 0th-, 1st-, 2nd-, 3rd- and 4th-order nested commutators
		C0 = id_mat
		C1 = bool_matmul(matrix,C0)
		C2 = bool_matmul(matrix,C1)
		C3 = bool_matmul(matrix,C2)
		C4 = bool_matmul(matrix,C3)

		# Accumulation of terms from 1, 2, 3, or 4 commutations
		a1 = C1
		a2 = or_arrays(C1,C2)
		a3 = or_arrays(a2,C3)
		a4 = or_arrays(a3,C4)

		# Specifically which kinds of terms will produce zero, single or double excitations after 1, 2, 3, or 4 commutations
		A1 = vec_to_column_matrix(or_arrays(a1[K],or_arrays(a1[E],a1[EE])))
		A2 = vec_to_column_matrix(or_arrays(a2[K],or_arrays(a2[E],a2[EE])))
		A3 = vec_to_column_matrix(or_arrays(a3[K],or_arrays(a3[E],a3[EE])))
		A4 = vec_to_column_matrix(or_arrays(a4[K],or_arrays(a4[E],a4[EE])))

		textlog("Those terms that are present after 0, 1, 2, or 3 commutations that may lead to zero, single, or double excitations during the remaining commutations")
		textlog()
		textlog("          ", 0, 1, 2, 3, sep="   ")
		textlog()
		Z0 = and_arrays(H0,A4)
		Z1 = and_arrays(H1,A3)
		Z2 = and_arrays(H2,A2)
		Z3 = and_arrays(H3,A1)
		for i in range(dim):
			textlog(labels[i], "   ",
			mask(Z0[i]), mask(Z1[i]), mask(Z2[i]), mask(Z3[i]),
			sep="   ")
		textlog()




		t0 = H0
		P0 = and_arrays(t0,A4)
		D0 = xor_arrays(t0,P0)

		t1 = bool_matmul(matrix,P0)
		P1 = and_arrays(t1,A3)
		D1 = xor_arrays(t1,P1)

		t2 = bool_matmul(matrix,P1)
		P2 = and_arrays(t2,A2)
		D2 = xor_arrays(t2,P2)

		t3 = bool_matmul(matrix,P2)
		P3 = and_arrays(t3,A1)
		D3 = xor_arrays(t3,P3)

		textlog("Those terms that will be present after 0, 1, 2, or 3 commutations (accounting for prior truncations of things that will not matter), which can be discarded for future commutators")
		textlog("Clearly, excitations and the constants should be accumulated and saved elsewhere!")
		textlog()
		textlog("          ", 0, 1, 2, 3, sep="   ")
		textlog()
		for i in range(dim):
			textlog(labels[i], "   ",
			mask(D0[i]), mask(D1[i]), mask(D2[i]), mask(D3[i]),
			sep="   ")
		textlog()

		######################################################## printing done --- for export ########################################################

		self.excitations = vector()
		for block in [K, E, EE]:  self.excitations[block] = 1
		D0 = column_matrix_to_vec(D0)
		D1 = column_matrix_to_vec(D1)
		D2 = column_matrix_to_vec(D2)
		D3 = column_matrix_to_vec(D3)
		for block in [K, E, EE]:
			D0[block] = 0
			D1[block] = 0
			D2[block] = 0
			D3[block] = 0
		Z0 = column_matrix_to_vec(Z0)
		Z1 = column_matrix_to_vec(Z1)
		Z2 = column_matrix_to_vec(Z2)
		Z3 = column_matrix_to_vec(Z3)
		self.keep  = [Z0, Z1, Z2, Z3]
		self.trash = [D0, D1, D2, D3]


if __name__ == "__main__":
	load_truncation_info(print)
