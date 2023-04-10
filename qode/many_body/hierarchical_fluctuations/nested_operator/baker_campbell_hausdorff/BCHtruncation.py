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
The stuff in this file translates the generic concepts behind the BCH truncation into a form that cooperates well with nested_operator
"""

from ..... import util
from ....coupled_cluster.commutator_scatter_2_2 import load_truncation_info
from .. import mask



class BCHtruncation(object):
	def __init__(self, CC_order, H_order):
		if CC_order!=2 or H_order!=2:  raise Exception("Sorry, BCH code is only for CCSD on a 2-body operator right now!")
		self.max_commutator = 4
	def masks(self, commutator_order):
		#info = load_truncation_info(print)
		info = load_truncation_info()
		excitations = mask(encode_nested(info.excitations))
		if commutator_order<4:
			to_keep     = mask(encode_nested(info.keep[commutator_order]))
			#print_nested_mask(to_keep)
			to_trash    = mask(encode_nested(info.trash[commutator_order]))
		else:
			to_keep     = mask(False)
			to_trash    = mask(False)
		return util.output(excitations=excitations, to_keep=to_keep, to_trash=to_trash)



def empty():  return [False, False, False, False]

def true():   return [True,  False, False, False]

def encode_nested(inp):
	K, E, F, D, EE, EF, ED, FF, FD, DD, E1EE, E1EF, E1ED, E1FF, E1FD, E1DD, E2EE, E2EF, E2ED, E2FF, E2FD, E2DD, E3EE, E3EF, E3ED, E3FF, E3FD, E3DD, E4EE, E4EF, E4ED, E4FF, E4FD, E4DD = range(34)
	out = empty()
	if inp[K]:  out[K] = True

	if inp[F] or inp[FF] or inp[FD]:  out[F] = empty()
	if inp[ F]:  out[F][K] = True
	if inp[FF]:  out[F][F] = true()
	if inp[FD]:  out[F][D] = true()

	if inp[D] or inp[DD]:             out[D] = empty()
	if inp[ D]:  out[D][K] = True
	if inp[DD]:  out[D][D] = true()

	if (inp[E] or
	    inp[EE] or inp[EF] or inp[ED] or
	    inp[E1EE] or inp[E1EF] or inp[E1ED] or inp[E1FF] or inp[E1FD] or inp[E1DD] or
	    inp[E2EE] or inp[E2EF] or inp[E2ED] or inp[E2FF] or inp[E2FD] or inp[E2DD] or
	    inp[E3EE] or inp[E3EF] or inp[E3ED] or inp[E3FF] or inp[E3FD] or inp[E3DD] or
	    inp[E4EE] or inp[E4EF] or inp[E4ED] or inp[E4FF] or inp[E4FD] or inp[E4DD]):
		out[E] = empty()

	if (inp[EE] or
	    inp[E1EE] or inp[E1EF] or inp[E1ED] or
	    inp[E2EE] or inp[E2EF] or inp[E2ED] or inp[E2FF] or inp[E2FD] or inp[E2DD] or
	    inp[E3EE] or inp[E3EF] or inp[E3ED] or inp[E3FF] or inp[E3FD] or inp[E3DD] or
	    inp[E4EE] or inp[E4EF] or inp[E4ED] or inp[E4FF] or inp[E4FD] or inp[E4DD]):
		out[E][E] = empty()

	if (inp[E1EE] or
	    inp[E2EE] or inp[E2EF] or inp[E2ED] or
	    inp[E3EE] or inp[E3EF] or inp[E3ED] or inp[E3FF] or inp[E3FD] or inp[E3DD] or
	    inp[E4EE] or inp[E4EF] or inp[E4ED] or inp[E4FF] or inp[E4FD] or inp[E4DD]):
		out[E][E][E] = empty()

	if (inp[E2EE] or
	    inp[E3EE] or inp[E3EF] or inp[E3ED] or
	    inp[E4EE] or inp[E4EF] or inp[E4ED] or inp[E4FF] or inp[E4FD] or inp[E4DD]):
		out[E][E][E][E] = empty()

	if (inp[E3EE] or
	    inp[E4EE] or inp[E4EF] or inp[E4ED]):
		out[E][E][E][E][E] = empty()

	if (inp[E4EE]):
		out[E][E][E][E][E][E] = empty()

	if inp[   E]:  out[E][K] = True
	if inp[  EE]:  out[E][E][K] = True
	if inp[E1EE]:  out[E][E][E][K] = True
	if inp[E2EE]:  out[E][E][E][E][K] = True
	if inp[E3EE]:  out[E][E][E][E][E][K] = True
	if inp[E4EE]:  out[E][E][E][E][E][E][K] = True

	if inp[  EF] or inp[E1FF] or inp[E1FD]:  out[E][F] = empty()
	if inp[  EF]:  out[E][F][K] = True
	if inp[E1FF]:  out[E][F][F] = true()
	if inp[E1FD]:  out[E][F][D] = true()
	if inp[  ED] or inp[E1DD]:  out[E][D] = empty()
	if inp[  ED]:  out[E][D][K] = True
	if inp[E1DD]:  out[E][D][D] = true()

	if inp[E1EF] or inp[E2FF] or inp[E2FD]:  out[E][E][F] = empty()
	if inp[E1EF]:  out[E][E][F][K] = True
	if inp[E2FF]:  out[E][E][F][F] = true()
	if inp[E2FD]:  out[E][E][F][D] = true()
	if inp[E1ED] or inp[E2DD]:  out[E][E][D] = empty()
	if inp[E1ED]:  out[E][E][D][K] = True
	if inp[E2DD]:  out[E][E][D][D] = true()

	if inp[E2EF] or inp[E3FF] or inp[E3FD]:  out[E][E][E][F] = empty()
	if inp[E2EF]:  out[E][E][E][F][K] = True
	if inp[E3FF]:  out[E][E][E][F][F] = true()
	if inp[E3FD]:  out[E][E][E][F][D] = true()
	if inp[E2ED] or inp[E3DD]:  out[E][E][E][D] = empty()
	if inp[E2ED]:  out[E][E][E][D][K] = True
	if inp[E3DD]:  out[E][E][E][D][D] = true()

	if inp[E3EF] or inp[E4FF] or inp[E4FD]:  out[E][E][E][E][F] = empty()
	if inp[E3EF]:  out[E][E][E][E][F][K] = True
	if inp[E4FF]:  out[E][E][E][E][F][F] = true()
	if inp[E4FD]:  out[E][E][E][E][F][D] = true()
	if inp[E3ED] or inp[E4DD]:  out[E][E][E][E][D] = empty()
	if inp[E3ED]:  out[E][E][E][E][D][K] = True
	if inp[E4DD]:  out[E][E][E][E][D][D] = true()

	if inp[E4EF]:  out[E][E][E][E][E][F] = true()
	if inp[E4ED]:  out[E][E][E][E][E][D] = true()

	return out



def print_nested(encoding, prestring=""):
	K, E, F, D = range(4)
	chars = ["K", "E", "F", "D"]
	if    encoding[K] and prestring=="":  print(chars[K], "  ", end="")
	elif  encoding[K]:  print(prestring, "  ", end="")
	for i in [E,F,D]:
		if encoding[i] is not False:  print_nested(encoding[i], prestring+chars[i])
	if prestring=="":  print("\n")


def print_nested_mask(the_mask, prestring=""):
	if prestring=="":  print("PRINTING FROM MASK\n")
	if   the_mask.C_on and prestring=="":  print("K  ", end="")
	elif the_mask.C_on:                    print(prestring, "  ", end="")
	if the_mask.Ex_on:  print_nested_mask(the_mask.Ex_mask, prestring+"Ex")
	if the_mask.Fo_on:  print_nested_mask(the_mask.Fo_mask, prestring+"Fo")
	if the_mask.Fv_on:  print_nested_mask(the_mask.Fv_mask, prestring+"Fv")
	if the_mask.Dx_on:  print_nested_mask(the_mask.Dx_mask, prestring+"Dx")
	if prestring=="":  print("\n")
