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
import sys
import numpy as np
from qode.fermion_field          import occ_strings
from qode.fermion_field.state    import configuration, state, dot
from qode.fermion_field.reorder2 import *
from hamiltonian2 import hamiltonian



occ = [0,4]
vrt = [1,2,3, 5,6,7]

all_occ_strings = occ_strings.CISD(occ,vrt)	# Works because CISD=FCI for 2 e- (else basis projections might make equivalent ops look different)
num_config = len(all_occ_strings)

h_mat = np.load('data/h_mat.npy')
V_mat = np.load('data/V_mat.npy')
H = hamiltonian(0.0, h_mat, V_mat, all_occ_strings, occ, vrt)





#all_occ_strings = occ_strings.all_occ_strings(8,4)
#occ = [0,2,4,6]
#vrt = [1,3,5,7]


all_e = []
for i in occ:
	for a in vrt:
		all_e += [excitation(a,i)]
all_d = []
for i in occ:
	for a in vrt:
		all_d += [deexcitation(a,i)]
all_fo = []
for i in occ:
	for j in occ:
		all_fo += [occ_rearrange(i,j)]
all_fv = []
for a in vrt:
	for b in vrt:
		all_fv += [vrt_rearrange(a,b)]



def test_mult(op, trans_list, all_occ_strings):
	error = 0.
	print("LEFT MULTIPLY")
	for t in trans_list:
		print(t)
		for i in range(len(all_occ_strings)):
			printable = str(configuration(all_occ_strings[i]))
			psi_in = state(all_occ_strings,i)
			control = t(op(psi_in))
			experiment = op.left_multiply(t)(psi_in)
			#
			control.increment(experiment,-1)
			err = dot(control,control)
			printable += "         error = {}".format(err)
			print(printable)
			error += err
			if err>1e-10:  raise Exception("inequivalent operators")
			sys.stdout.flush()
	"""\
	print("RIGHT MULTIPLY")
	for t in trans_list:
		print(t)
		for i in range(len(all_occ_strings)):
			printable = str(configuration(all_occ_strings[i]))
			psi_in = state(all_occ_strings,i)
			control = op(t(psi_in))
			experiment = op.right_multiply(t)(psi_in)
			control.increment(experiment,-1)
			err = dot(control,control)
			printable += "         error = {}".format(err)
			print(printable)
			error += err
			if err>1e-10:  raise Exception("inequivalent operators")
	"""
	return error



print("##################   E    ##################")
error = test_mult(H, all_e, all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   D    ##################")
error = test_mult(H, all_d, all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   Fo   ##################")
error = test_mult(H, all_fo, all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   Fv   ##################")
error = test_mult(H, all_fv, all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
