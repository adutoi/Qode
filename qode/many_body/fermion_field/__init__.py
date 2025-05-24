#    (C) Copyright 2025 Anthony Dutoi
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
from . import field_op
from .field_op_ham    import Hamiltonian
from .CI_space_traits import CI_space_traits
from .configurations  import dn_up_elec, combine_orb_lists, all_configs, Sz_configs, decompose_configs, recompose_configs, config_combination, tensor_product_configs, print_configs
from . import CI_methods

# Contains direct interfaces to low-level (bit-fiddling) code for operating on configutations strings
#   field_op.py
# The aforementioned low-level (bit-fiddling) code for operating on configutations strings
#   field_op.c
#   antisymm.c

# The next layer up builds on field_op.py (etc) but do not know about each other
#   field_op_ham.py       given appropriate sets if integrals resolves the action of the Hamiltonian on a list of coefficients given low-level config reps
#   CI_space_traits.py    connects arrays of coefficients to low-level reps of configurations for building state vectors

# An isolated module (no dependencies on foregoing) that generates lists of configs represented as integers.
#   configurations.py

# Built directly on all but lowest layer above, plus high-level numerics (ie lanczos) and an understanding of fragment/molecule data structures
#   CI_methods.py

# An important change of perspective happens between the low-level routines in field_op.py (etc) and configurations.py (and somewhat CI_methods).
#   In the low-level routines, we imagine integers written in their usual way, with low-order bits to the right because we will be constantly 
# manipulating them explicitly.  To match this, we can therefore imagine arrays as written right-to-left (opposite of the usual indo-european language
# language convention) so that indices and orders match in spatial representation as well and numerical values.  This is fine because we never explicitly
# write out any arrays at this level, but rather, just call elements by index.
#   In the high-level routines, we imagine arrays written in their usual way, with low-index elements to the left because we will constantly be
# writing arrays out explicitly with the [] notation.  To match this, we can therefore imagine integers written with their low-order bits to the 
# left (opposite of the usual now-worldwide arabic number convention) so that indices and orders match in spatial representation as well and numerical
# values.  This is fine because we never expicitly manipulate the bits of these integers.  This will be helpful in conceptualizing the order of
# some loop nesting and the order in which some arrays are written (when the user knows that they come from the decompostion of an integer).
