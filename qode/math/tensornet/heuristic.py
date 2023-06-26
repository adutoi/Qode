#    (C) Copyright 2023 Anthony D. Dutoi
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

import math



# This is separated out from the calling code to encapsulate the potentially fragile
# algorithm (and make clear its dependencies), away from the more robust logic surrounding it.
def heuristic(scalar, contraction_groups, free_indices_groups, shapes):
    all_groups = {group for group in contraction_groups} | {group for group in free_indices_groups}
    reductions = {}
    for group in all_groups:
        total_size, product_size = 0, 1
        for tens in group:
            tens_size = math.prod(shapes[tens])
            total_size   += tens_size
            product_size *= tens_size    # first get size of the outer product before reductions
        if group in contraction_groups:
            for contraction in contraction_groups[group]:
                tens0, pos0 = contraction[0]    # always instantiated by function that checks congruence
                axis_length = shapes[tens0][pos0]
                product_size = product_size // axis_length**len(contraction)
        if group in free_indices_groups:
            for free_index in free_indices_groups[group]:
                tens0, pos0 = free_index[0]     # always instantiated by function that checks congruence
                axis_length = shapes[tens0][pos0]
                product_size = product_size // axis_length**(len(free_index)-1)
        reductions[group] = product_size / total_size
    do_scalar_mult, do_reduction, target = False, True, None
    if len(reductions)>0:
        target = min(reductions, key=reductions.get)
        if reductions[target]>1 and scalar!=1:
            do_scalar_mult, do_reduction = True, False
    elif scalar!=1 and len(shapes)>0:
        do_scalar_mult, do_reduction = True, False
    else:
        do_reduction = False    # only thing left is to take outer product or return scalar
    if do_scalar_mult:
        target = min(shapes, key=lambda tens: math.prod(shapes[tens]))    # the smallest tensor
    return do_scalar_mult, do_reduction, target
