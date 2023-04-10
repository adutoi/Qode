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
import sys
import numpy as np
import numpy.linalg as LA
from time import ctime
from math import sqrt
import copy
#from . import analytical_index
from qode.coupled_oscillators.analytical import analytical_lowest_states

from qode.coupled_oscillators.analytical.analytical_routines \
          import dimensionless_oscillator_energy,            \
                 unitless_oscillator_energy,                 \
                 unitless_HO_system_energy

from qode.coupled_oscillators.analytical.analytical_recursive_states import get_lowest_n_states
#
# Analytical Solution to Harmonic Oscillator Systems 
#
#


# A set of n harmonic oscillators has the following Hamiltonian:
# 
# H =  Sum_i h(i) + Sum_i<j kij xi xj
#
#   =  Sum_i    hb^2 / mi * d^2/dxi^2   +    "Kinetic   Energy"
#      Sum_i    1/2 kii xi^2            +    "Potential Energy"
#      Sum_i<j  kij xi xj                    "Linear    Coupling"
#


# 1) Transformations to Unitless Hamiltonian:
# Let  xi    =  sqrt( m0 / mi ) * L0 * xi~
#      xi^2  =  m0/mi * L0^2 * xi~
#      xj    =  sqrt( m0 / mj ) * L0 * xj~
#
#
# H = Sum_i     hb**2 / mi * d^2 / (m0/mi * L0^2 )dxi~^2        +
#     Sum_i     1/2 * kii * m0/mi * L0^2 * xi~^2                +
#     Sum_i!=j  1/2 * kij * m0*L0^2 / sqrt( mi*mj ) xi~ * xj~
#
#   = hb^2 / (m0 L0^2) * [ Sum_i    (-1/2 * d^2/dxi~^2 + 1/2 * kii~ xi~^2 ) + 
#                          Sum_i,j    1/2 * kij~ xi~ xj~                      ]
#
# hb^2/(m0 L0^2) is the energy prefactor.
#
#
# Therefore,
# H = hb^2 / (m0 L0^2) * [ Sum_i    (-1/2 * d^2/dxi~^2) +  
#                          Sum_i,j  ( 1/2 * kij~ xi~ xj~ )    ]    
# where, 
#       kii~ = kii * (m0L0^2)^2 / ( hb^2 * mi )
#       kij~ = kij * (m0L0^2)^2 / ( hb^2 * sqrt( mi * mj ) )


# 2) Diagonalization of Coupling Matrix (K):
# ( Transformations to effective force constants coordinates )
# 
# E, U = diag( K )
# 
# Therefore,
# H = hb^2 / (m0 L0^2) * Sum_i [-1/2 * d^2/dyi^2 + 1/2 Eii yi^2 ] 
# 
# where, Eii are eigenvalues of coupling matrix K.
#


# 3) Transformations to Dimensionless Hamiltonian:
# Let  yi = zi * Eii^(-1/4) 
#
# H = hb^2 / (m0 L0^2) * Sum_i sqrt(Eii) * { -1/2 * d^2/dzi^2 + 1/2 zi^2 }
# 
# Eigenvalues are:
# Eni(i) = hb^2 / (m0 L0^2) * sqrt(Eii) * ( ni + 1/2 ) * hb
# 
# For each configuration, the total energy = Sum_i Eni(i) 
#


def build_k_mat(
                 primitive_k_mat,  
                 mass_list,  
                 unit_length, 
                 unit_mass,  
                 unit_hbar 
               ):
    # INPUT:   primitive k (coupling) matrix , masses, unit mass, unit length, unit hbar
    #          k matrix consists of primitive k values 
    # OUTPUT:  tranformed k matrix ready to diagonalize

    # new k matrix formulas:
    #   kii~ = kii * (m0L0^2)^2 / ( hb^2 * mi ) 
    #        = (m0L0^2)^2/hb^2 * kii / mi 
    #        = shared_factor   * kii / mi 
    #
    #   kij~ = kij * (m0L0^2)^2 / ( hb^2 * sqrt( mi * mj ) )    
    #        = (m0L0^2)^2/hb^2 * kij / sqrt( mi*mj ) 
    #        = shared_factor   * kij / sqrt( mi*mj )

    # make a deepcopy of primitive k matrix
    primitive_k_mat = copy.deepcopy( primitive_k_mat )
    primitive_k_mat = np.matrix( primitive_k_mat )

    # Initialize new k matrix to be an empty numpy matrix:
    lda = len( primitive_k_mat )
    k_mat_new  = np.zeros(( lda, lda ))
    k_mat_new = np.matrix( k_mat_new )

    # Compute this shared_fatcor = (m0L0^2)^2/hb^2 
    shared_factor = pow( unit_mass * unit_length**2, 2 ) / unit_hbar**2

    # Fill up new k matrix elements:
    # Looping over i, j in the range of lda.
    for i in range( lda ):
        for j in range( lda ):
            if i != j:
                k_mat_new[i,j] = shared_factor * primitive_k_mat[i,j] / sqrt( mass_list[i] * mass_list[j] ) 
            else:
                k_mat_new[i,i] = shared_factor * primitive_k_mat[i,i] / mass_list[i]
    return k_mat_new






#
# Analytical Solution Code:
#
# import module "analytical_index"
# initialize object analytical_calc,
# pass energy upper bound,
#      energy function and its dependencies
#      to the analytical_calc()
#
'''def dimensionless_oscillator_energy( state_number ):
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
    return sys_energy'''


def unitless_analytical_solution(
                                 effective_k_list, 
                                 lowest_num_states
                                 ):
    # should return both state numbers and eigenvalues.
    num_oscillators = len( effective_k_list )
    
    #states, eigen_energy = analytical_lowest_states.get_lowest_n_states(
    #                                                                    effective_k_list,
    #                                                                lowest_num_states
    #                                                                )
    states, eigen_energy = get_lowest_n_states(                                            
                                               lowest_num_states,
                                               effective_k_list
                                               )
    
    
    # states and eigen_energy are sorted ( according to energy )
    return states, eigen_energy








#
#
# Module main function:
#
#

def analytical_main(   
                      mass_list,
                      primitive_k_mat,
                      lowest_num_states,
                      unit_length, 
                      unit_mass, 
                      unit_hbar
                    ):
    # INPUT:   masses, primitive k matrix
    #          Unit length, unit mass, unit hbar 
    # 
    # OUTPUT:  Sorted analytical eigenvalues. 
    #          State Numbers are available but not returned. 
    #
    #print("\n\n==================================================================")
    #print("ANALYTICAL SOLUTION:")
    #print("##########################################################")
    # Build Coupling Matrix for diagonalization
    K = build_k_mat( primitive_k_mat,  mass_list , 
                     unit_length    ,  unit_mass ,  unit_hbar )

    #print("\nK:\n",K)
    #
    # Diagonalization of K matrix.
    E,U = LA.eigh(K, UPLO='U')
    vectors = U
 
    # eigenvalues are effective k:
    k_eff = E.tolist()

    #print('\nEffective Force Constant:\n',k_eff )
    #print('\nUnitary Matirx:\n',U )


    # Sanity Check. Effective k should all be positive (if not, couplings are too strong)
    pass_eigenvalue_check = True
    
    n = 0
    num_k = len( k_eff )
    while pass_eigenvalue_check and n < num_k:
        if k_eff[n] < 0.0:
            print("Couplings Too Strong, Eigenvalue(s) Turn Negative!")
            pass_eigenvalue_check = False
        else:
            n += 1

    #print("Passed Effective Force Constants Checking: They are all positive.")

    #
    # Energy prefactor ( Unit of Energy ):
    pfactor = unit_hbar**2 / (unit_mass * unit_length**2)
    #print("\n\nEnergy Pre-Factor =",pfactor,"for Analytical Solution.")
    #
 
    # 
    # Carry Out Analytical Calculations
    # Return Sorted Energies
    #
    if pass_eigenvalue_check == True:
        #print("ANALYTICAL SOLUTION STARTS AT", ctime() )
        states_configure, E_analytical = unitless_analytical_solution( k_eff, lowest_num_states )
        #E_ana_sorted = np.sort( E_analytical )
        if pfactor != 1.0:
            print("WARNING! ENERGY PREFACTOR IS NOT EQUAL TO ONE!\n\n")
        #print("\n\n##########################################################")
        #print("ANALYTICAL SOLUTION DONE AT", ctime() )
        #print("==================================================================\n\n")
        # Return list_of_k, all indecies, and eigen energyies.
        return k_eff, states_configure, E_analytical, vectors
    else:
        errors  = np.zeros((1,1)) 
        print("\n\n==================================================================")
        print("ANALYTICAL SOLUTION FAILED!")
        print("==================================================================\n\n")
        return errors, errors, errors, errors 





