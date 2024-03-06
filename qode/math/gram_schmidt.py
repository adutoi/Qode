#    (C) Copyright 2024 Anthony D. Dutoi
# 
#    This file is part of QodeApplications.
# 
#    QodeApplications is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    QodeApplications is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with QodeApplications.  If not, see <http://www.gnu.org/licenses/>.
#
from .field_traits import machine_epsilon, sqrt

def gram_schmidt(vectors, normalize=False, normalized=False, n_times=1):
    if len(vectors)<2:                           # If one (or fewer) ...
        if normalize:                            # ... vectors, ...
            for v in vectors:  v /= sqrt(v|v)    # ... renormalize(?) ...
        return                                   # ... and leave.

    eps = machine_epsilon[vectors[0].field]    # to avoid division by zero (so result not properly normalized if linear dependency)

    half = len(vectors) // 2              # If we do the GS ...
    order  = list(range(len(vectors)))    # ... multiple times ...
    orderA = order[:half]                 # ... we will keep ...
    orderB = order[half:]                 # ... shuffling the order.

    count = 0                             # Just to keep track of shuffling
    for _ in range(n_times):
        # Precompute squared norms for efficiency (unless we know they are normalized)
        if normalized:
            norms2 = [1.] * len(vectors)
        else:
            norms2 = [(v|v) for v in vectors]

        # The main Gram-Schmidt loop (but perhaps in shuffled order after the first time)
        for n,i in enumerate(order):
            for j in order[:n]:
                vectors[i] -= vectors[j] * ((vectors[j]|vectors[i]) / (norms2[j]+eps))    # extra parentheses so division is on scalar

        # Normalize the now (approximately) orthogonal vectors, if requested
        normalized = False    # This is the default assumption, since vectors changed
        if normalize:
            for v in vectors:  v /= (sqrt(v|v)+eps)
            normalized = True

        # Shuffle the ordering
        if count%2:                               # If count is odd, prepare an even-numbered loop.
            orderA = orderA[-1:] + orderA[:-1]    # ... Each even-numbered loop brings ...
            orderB = orderB[1:]  + orderB[:1]     # ... vectors from the center to ...
            order  = orderA + orderB              # ... the outsides of the list.
        else:
            order = list(reversed(order))         # Odd loops do the even orderings in reverse.
        count += 1

    return
