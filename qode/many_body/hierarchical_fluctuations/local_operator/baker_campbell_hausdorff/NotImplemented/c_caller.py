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
# This could be done better, starting with embedding this stuff into the autoQode folder by changing import paths etc.
# Even better though would be to force autogeneration of the code (maybe . . . it takes a long time)

instructions = """\
The module you are calling relies on some automatically generated C code.
Got to:
    Qode/qode/many_body/local_operator/baker_campbell_hausdorff
and run:
    ./run_autoqode.sh
Sorry it takes quite some time (~1 hr).  Then retry what you are doing.  This error should go away.
"""

def compute_next_X( prev_X_obj, prev_T_obj, new_X_obj, rec_num_states, n_threads ):
	print(instructions)
	raise NotImplmentedError
