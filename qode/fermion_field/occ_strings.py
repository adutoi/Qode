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
from copy import copy



# The recursive workhorse and its user-friendly wrapper

def _make_occ_strings(occ_string,n_elec,preface=[]):
	"""\
	This is a workhorse function, not likely to be of use to a caller outside of this file.

	Given an input occ_string (boolean list of occupancies), this will add n_elec more electrons
	to only the empty slots, in all possible configruations in those empty slots.  The results
	are returned as a list with the preface occupation string pre-pended to it.  (The preface is
	needed for the recursive manner in which this function does its work, and is not meant to
	be specified from the "outside.")
	"""
	the_list = []
	if n_elec==0:
		the_list += [preface+occ_string]
	else:
		n_elec -= 1
		for i in range(len(occ_string)):
			if not occ_string[i]:
				occ_string_temp = copy(occ_string)
				occ_string_temp[i] = True
				the_list += _make_occ_strings(occ_string_temp[i+1:],n_elec,preface+occ_string_temp[:i+1])
	return the_list

def all_occ_strings(N,n):
	"""\
	This wrapper is suitable for generating all FCI occupation strings for n electrons in N (spin) orbitals.
	It is also useful for generating the separate particle and hole configurations for a given excitation order
	above a reference configuration.  It returns a list of occupation strings, each of which is a list of
	boolean occupancies.
	"""
	return _make_occ_strings([False]*N,n)



# Helpers for generating occupation strings of a certain excitation order, given a reference implied by the occupied set

def _sanity_check(occ_orbs,vrt_orbs):
	"""\
	occ_orbs and vrt_orbs are lists of integer indices of occupied and virtual orbitals, respectively.
	Make sure orbitals all belong to contiguous index block and that the sets are non-overapping.
	"""
	n_orbs = len(occ_orbs) + len(vrt_orbs)
	for i in occ_orbs:
		if i not in range(n_orbs):  raise Exception("occupied index out of range")
	for a in vrt_orbs:
		if a not in range(n_orbs):  raise Exception("virtual index out of range")
	for i in occ_orbs:
		if i in vrt_orbs:  raise Exception("overlapping occupied and virtual sets")
	return n_orbs

def _invert_ph(occ_strings):
	""" Takes a set of occupation strings and turns each hole (False) into a particle (True) and each particle into a hole.  Returns a new list. """
	new_occ_strings = []
	for occ_string in occ_strings:
		new_occ_string = [not orb for orb in occ_string]
		new_occ_strings += [new_occ_string]
	return new_occ_strings

def _n_tuple_substitutions(n,occ_orbs,vrt_orbs):
	"""\
	Generates all occupation strings associated with given exciation order n above a reference
	implied by the occupied set, in an orbital space that is the union of the occupied and virtual
	sets (both lists of integer orbital indices).  Returns the result as a list of occupation
	strings, each of which is a list of boolean occupancies.
	"""
	n_orbs = _sanity_check(occ_orbs,vrt_orbs)
	particle_occ_strings = all_occ_strings(len(vrt_orbs),n)
	hole_occ_strings     = all_occ_strings(len(occ_orbs),n)
	hole_occ_strings = _invert_ph(hole_occ_strings)
	occ_strings = []
	for hole_occ_string in hole_occ_strings:
		for particle_occ_string in particle_occ_strings:
			occ_string = [None]*n_orbs
			for i,TF in enumerate(    hole_occ_string):  occ_string[occ_orbs[i]] = TF
			for a,TF in enumerate(particle_occ_string):  occ_string[vrt_orbs[a]] = TF
			occ_strings += [occ_string]
	return occ_strings



# CI occupation string sets

def CI(orders,occ_orbs,vrt_orbs):
	"""\
	occ_orbs and vrt_orbs are lists containing the indices of occupied and virtual spin orbitals in any order.

	This returns a list of occupation strings (each of which is a list of boolean occupancies) representing all excitations
	of all orders in orders (a list of integers) above a reference implied by the occupied set, in an orbital space that is
	the union of the occupied and virtual sets.

	The blocking by excitation order will be in the order given in the input list, but the ordering of occupation
	strings associated with a given excitation order is defined internally by the algorithm.

	This basically functions as a more sophisticated wrapper for _n_tuple_substitutions().
	"""
	occ_strings = []
	for n in orders:
		if n>len(occ_orbs):  raise Exception("CI excitation order exceeds available number of electrons")
		occ_strings += _n_tuple_substitutions(n,occ_orbs,vrt_orbs)
	return occ_strings

def CIS(   occ_orbs,vrt_orbs):  return CI([0,1],      occ_orbs,vrt_orbs)
def CISD(  occ_orbs,vrt_orbs):  return CI([0,1,2],    occ_orbs,vrt_orbs)
def CISDT( occ_orbs,vrt_orbs):  return CI([0,1,2,3],  occ_orbs,vrt_orbs)
def CISDTQ(occ_orbs,vrt_orbs):  return CI([0,1,2,3,4],occ_orbs,vrt_orbs)

def FCI(occ_orbs,vrt_orbs):	# Not blocked by excitation order, and pays no heed to orbital ordering!
		n_orbs = _sanity_check(occ_orbs,vrt_orbs)
		return all_occ_strings(n_orbs,len(occ_orbs))
