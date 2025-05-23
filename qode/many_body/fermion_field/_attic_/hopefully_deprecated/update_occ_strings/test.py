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
from qode.fermion_field import config, occ_strings



# Test raw occ_string code, in particular generation of all configurations

print("FCI occupation strings for 2 electrons in 4 spin orbitals")
strings = occ_strings.all_occ_strings(4,2)
for string in strings:  print(strings.index(string),string)
print("LENGTH =",len(strings))



# This uses the config module but really only tests the __init__ function and the __str__ function
# The configurations mimuc CISD on a 2e- system (so FCI) where there are 4 spatial orbitals blocked into alpha and beta spin.

occ_orbs = [0, 4]
vrt_orbs = [1,2,3, 5,6,7]

print("CISD occupation strings for 2 electrons in 8 spin orbitals blocked by alpha and beta")
cisd_strings = occ_strings.CISD(occ_orbs,vrt_orbs)
for configuration in cisd_strings:  print(config.configuration(configuration))
print("LENGTH =",len(cisd_strings))



# Some tests making many-e- configurations where the locations of the electrons are restricted to subsets that might be overlapping

strings = occ_strings.all_occ_strings(8,2,[[0,1,2,3],[4,5,6,7]])
for configuration in strings:  print(config.configuration(configuration))
print("LENGTH =",len(strings))

strings = occ_strings.all_occ_strings(8,2,[[0,2,3,4,6],[1,3,4,5,7]])
for configuration in strings:  print(config.configuration(configuration))
print("LENGTH =",len(strings))


# This test breaks!  

strings = occ_strings.all_occ_strings(9,3,[[2,3,4,6],[1,3,4,5,7],[0,3,5,6,8]])
for configuration in strings:  print(config.configuration(configuration))
print("LENGTH =",len(strings))

"""\
here is the output.  Note that   | . . . 3 4 5 . . . >   occurs twice because
when 3 gets eliminated from the possible choices for the 2nd electron (when
the 1st is in state 4), there is no record indicating that the 3rd electron
cannot be in 3.

Need to be careful in fixing this so as not to eliminate configurations like
| . 1 . 3 4 . . . . >  where electron 3 is in state 3 even though electron 1 is
in state 4 (because electron 2 is in state 1 which is not overlapping other sets).

+1.000   x   | 0 1 2 . . . . . . >
+1.000   x   | . 1 2 3 . . . . . >
+1.000   x   | . 1 2 . . 5 . . . >
+1.000   x   | . 1 2 . . . 6 . . >
+1.000   x   | . 1 2 . . . . . 8 >
+1.000   x   | 0 . 2 3 . . . . . >
+1.000   x   | . . 2 3 . 5 . . . >
+1.000   x   | . . 2 3 . . 6 . . >
+1.000   x   | . . 2 3 . . . . 8 >
+1.000   x   | 0 . 2 . 4 . . . . >
+1.000   x   | . . 2 3 4 . . . . >
+1.000   x   | . . 2 . 4 5 . . . >
+1.000   x   | . . 2 . 4 . 6 . . >
+1.000   x   | . . 2 . 4 . . . 8 >
+1.000   x   | 0 . 2 . . 5 . . . >
+1.000   x   | . . 2 . . 5 6 . . >
+1.000   x   | . . 2 . . 5 . . 8 >
+1.000   x   | 0 . 2 . . . . 7 . >
+1.000   x   | . . 2 3 . . . 7 . >
+1.000   x   | . . 2 . . 5 . 7 . >
+1.000   x   | . . 2 . . . 6 7 . >
+1.000   x   | . . 2 . . . . 7 8 >
+1.000   x   | 0 1 . 3 . . . . . >
+1.000   x   | . 1 . 3 . 5 . . . >
+1.000   x   | . 1 . 3 . . 6 . . >
+1.000   x   | . 1 . 3 . . . . 8 >
+1.000   x   | 0 . . 3 4 . . . . >
+1.000   x   | . . . 3 4 5 . . . >
+1.000   x   | . . . 3 4 . 6 . . >
+1.000   x   | . . . 3 4 . . . 8 >
+1.000   x   | 0 . . 3 . 5 . . . >
+1.000   x   | . . . 3 . 5 6 . . >
+1.000   x   | . . . 3 . 5 . . 8 >
+1.000   x   | 0 . . 3 . . . 7 . >
+1.000   x   | . . . 3 . 5 . 7 . >
+1.000   x   | . . . 3 . . 6 7 . >
+1.000   x   | . . . 3 . . . 7 8 >
+1.000   x   | 0 1 . . 4 . . . . >
+1.000   x   | . 1 . 3 4 . . . . >
+1.000   x   | . 1 . . 4 5 . . . >
+1.000   x   | . 1 . . 4 . 6 . . >
+1.000   x   | . 1 . . 4 . . . 8 >
+1.000   x   | 0 . . . 4 5 . . . >
+1.000   x   | . . . 3 4 5 . . . >
+1.000   x   | . . . . 4 5 6 . . >
+1.000   x   | . . . . 4 5 . . 8 >
+1.000   x   | 0 . . . 4 . . 7 . >
+1.000   x   | . . . 3 4 . . 7 . >
+1.000   x   | . . . . 4 5 . 7 . >
+1.000   x   | . . . . 4 . 6 7 . >
+1.000   x   | . . . . 4 . . 7 8 >
+1.000   x   | 0 1 . . . . 6 . . >
+1.000   x   | . 1 . . . 5 6 . . >
+1.000   x   | . 1 . . . . 6 . 8 >
+1.000   x   | 0 . . 3 . . 6 . . >
+1.000   x   | . . . 3 . 5 6 . . >
+1.000   x   | . . . 3 . . 6 . 8 >
+1.000   x   | 0 . . . 4 . 6 . . >
+1.000   x   | . . . . 4 5 6 . . >
+1.000   x   | . . . . 4 . 6 . 8 >
+1.000   x   | 0 . . . . 5 6 . . >
+1.000   x   | . . . . . 5 6 . 8 >
+1.000   x   | 0 . . . . . 6 7 . >
+1.000   x   | . . . . . 5 6 7 . >
+1.000   x   | . . . . . . 6 7 8 >
"""
