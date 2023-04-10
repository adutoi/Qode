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
top= """\
from scalars import Coeff
from containers import add, mult, comm, do_commutators, reorder_ops, reorder_pdts, nest_summations, resolve_deltas, expand, sum_case, trim_high_order
from operators import Ex, Fo, Fv, Dx, CrtOcc, CrtVrt, DstOcc, DstVrt
from indices import Indices
from summation import Sum



"""


bottom = """\



out = open("{}","w")

#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = expand(expression)				# commutators only work with primitives, so expand first
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = do_commutators(expression)			# take primitive commutators
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = sum_case(expression)			# need to do this before reordering, since reordering depends on case
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = expand(expression)				# first "flatten" products in summands (had a good reason for not doing this before sum_case ... maybe just more useful printing?)
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = reorder_ops(expression)			# primary purpose here is to collect together operators and reorder, also collects deltas and tensors together (just an easy by-product)
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = expand(expression)				# operator reordering could have resulted in multi-term summands, break out into many summations
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = trim_high_order(expression)		# eliminate unnecessary components, also collects deltas and tensors together (useful for next step)
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = resolve_deltas(expression)			# finally, resolve the deltas
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = expand(expression)				# this mainly strips zeros out of an otherwise very long addition of zero terms with scattered things we want
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = reorder_pdts(expression)			# collect together tensors and already-ordered operator strings (untouched)
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

expression = nest_summations(expression)		# separate math operations from storage operations
#out.write(str(expression))
#out.write("\\n\\\\\\\\\\n")

out.write(str(expression.multiline(self_contained=True)))	# only has an effect on top most "add"
#out.write("\\n")
out.close()
"""
