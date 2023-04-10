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
from math import sqrt

def dimensionless_oscillator_energy( state_number ):
    # E_dimensionless = ( n + 1/2 )
    return state_number + 0.5

def unitless_oscillator_energy( state_number, effect_k ):
    # E_unitless = sqrt( Eii ) * E_dimensionless
    return sqrt( effect_k ) * dimensionless_oscillator_energy( state_number )

def unitless_HO_system_energy( state_config, effect_k_list ):
    # E( n oscillators ) = Sum_i E_unitless
    # Each HO has a state read from state_config[ i ]
    # Each HO has a effective k read from effect_k_list[ i ]
    num_ho = len( effect_k_list )
    sys_energy = 0.0
    for i in range( num_ho ):
        sys_energy += unitless_oscillator_energy( state_config[ i ], effect_k_list[ i ] )
    return sys_energy



