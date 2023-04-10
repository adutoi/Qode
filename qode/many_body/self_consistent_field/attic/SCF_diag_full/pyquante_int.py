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
# Run under python2.7
#
# This Module is called by system calls in this format:
#
#     python2.7 pyquante_int.py Units Charge Multiplicity basis_set
#
# There must be two files in this folder:
#   1) Atomic_Num_List.npy  (format: [1,1,8] for example for H2O ) 
#   2) Coordinates_List.npy (format: [[0.0,0.0,0.0], [0.0,0.0,1.0], ...] )             
#
import sys
import os
import numpy
from PyQuante.Molecule import Molecule
from PyQuante.PyQuante2 import SCF
from PyQuante.cints import ijkl2intindex as intindex



def make_PyQuante_Mol_Obj( atom_num_list, coordinates_list, Units='Angstrom', Charge='0', Multiplicity='1' ):
	# print( atom_num_list )
	# print( coordinates_list )
	if len(atom_num_list) == len( coordinates_list ):
		pass
	else:
		raise AssertionError
	n_atom = len(atom_num_list)
	# Build PyQuante Molecule Object Dependencies...
	mol_list = []
	for i in range(n_atom):
		mol_list += [( atom_num_list[i], ( coordinates_list[i,0], coordinates_list[i,1], coordinates_list[i,2] ) )]
	# Initialization of PyQuante Molecule Object...
	mol_obj = Molecule( 'SCF_Calc', mol_list, units=Units, charge=Charge, multiplicity=Multiplicity )
	return mol_obj

def normalize_U(U,E):
	from math import sqrt
	ld1,ld2 = U.shape
	# print ld1,ld2, E.shape[0]
	if ld1 == E.shape[0] and ld1 == ld2:
		pass
	else:
		raise AssertionError
	for i in range(ld1):
		for j in range(ld1):
			U[j,i] /= sqrt(E[i])
	return U

# Example of PyQuante Molecule Object Construction:
# h2o = Molecule('H2O',
#                [(8,  ( 0.00000000,     0.00000000,     0.04851804)),
#                 (1,  ( 0.75300223,     0.00000000,    -0.51923377)),
#                 (1,  (-0.75300223,     0.00000000,    -0.51923377))],
#                 units='Angstrom', charge='0', multiplicity='1')

def Integral_main( atom_num_list, coordinates_list, Units, Charge, Multiplicity, basis_set ):
	Mol_Obj = make_PyQuante_Mol_Obj( atom_num_list, coordinates_list, Units, Charge, Multiplicity )
	# PyQuante 2.x
	solver = SCF(Mol_Obj, method="HF", basis=basis_set)
	nbfs = len(solver.basis_set.get()) # Num of Basis Function
	#
	# Overlap          Matrix: solver.integrals.S
	# Core Hamiltonian Matrix: solver.integrals.h 
	# Electronic Repulsion Matrix: solver.integrals.ERI
	#
	# E,U = numpy.linalg.eigh( numpy.matrix(solver.integrals.S) )
	# U = normalize_U(U,E)
	s_mat = numpy.matrix(solver.integrals.S)
	h_mat = numpy.matrix(solver.integrals.h) 
	v_pyquante = solver.integrals.ERI
	v_mat = numpy.zeros((nbfs**2,nbfs**2), dtype=numpy.float64)
	# Conversion from PyQuante to my v_mat indexings...
	for i in range(nbfs):
		for j in range(nbfs):
			ij = i*(i+1)/2+j
			for k in range(nbfs):
				for l in range(nbfs):
					kl = k*(k+1)/2+l
					v_mat[i*nbfs+k, j*nbfs+l] = v_pyquante[intindex(i,j,k,l)] 
	# solver.iterate(etol=1e-6)
	# print "\n\nHF Result = ", solver.energy
	# numpy.save('Eggen_Values.npy',E)
	# numpy.save('Normalized_Eigen_Vectors.npy',U)
	numpy.save('Overlap_Integral_Mat.npy',s_mat)
	numpy.save('Core_Hamiltonian_Mat.npy',h_mat)
	numpy.save('Elec_Repul_Mat.npy',v_mat)
	# print "Temporary Files Saved: Eggen_Values.npy, Normalized_Eigen_Vectors.npy, Core_Hamiltonian_Mat.npy, Elec_Repul_Mat.npy"

if __name__ == "__main__":
	Units  = sys.argv[1]
	Charge = sys.argv[2]
	Multiplicity = sys.argv[3]
	basis_set = sys.argv[4]
	atom_num_list    = numpy.load('Atomic_Num_List.npy')
	coordinates_list = numpy.load('Coordinates_List.npy')
	# Calling Integral Engine...
	Integral_main( atom_num_list, coordinates_list, Units, Charge, Multiplicity, basis_set ) 
	



