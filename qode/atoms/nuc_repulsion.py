#    (C) Copyright 2019 Anthony D. Dutoi
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
from ..util import constants
from .elements import atomic_number

energy_unit = constants.E_H	# In case the user needs it, this is the unit of energy for the return value (atomic units)



def energy(atoms, distance_unit=constants.Ang):	# by default, geometries assumed to be in Angstrom
    in_AU = float(distance_unit/constants.a_0)
    repulsion = 0
    for i,atom_i in enumerate(atoms):
        X_i, (xi,yi,zi) = atom_i  #("element position")    # atom_i is a struct
        Z_i = atomic_number(X_i)
        for atom_j in atoms[i+1:]:
            X_j, (xj,yj,zj) = atom_j  #("element position")    # atom_j is a struct
            Z_j = atomic_number(X_j)
            R = in_AU * math.sqrt((xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2)
            repulsion += (Z_i * Z_j / R)
    return repulsion
