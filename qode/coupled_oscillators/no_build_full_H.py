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
import analytical_solution
import numpy as np
import numpy.linalg as LA
import hamiltonian_class





def hamiltonian_multiply_vec( eigen_energy1, eigen_energy2, dipole_mat, vector ):
    if len( vector ) != len( dipole_mat ):
        print("ERROR in Hamiltonian and Vector Sizes")
        sys.exit(1)
    new_vector = [ 0.0 for i in range( len( vector ) ) ]

    x = 0
    while x < len(eigen_energy1) * len(eigen_energy2):
        y = 0
        for i in range(len(eigen_energy1)):
            for j in range(len(eigen_energy2)):
                new_vector[ x ] += dipole_mat[ x ][ y ]* vector[ y ]
                if x == y:
                    new_vector[ x ] += ( eigen_energy1[ i ] + eigen_energy2[ j ] ) * vector[y] 
                y += 1  
        x += 1             
    return new_vector        


if __name__ == "__main__":
    eigen1 = [ 1.0, 2.0, 1.4, 1.6]
    eigen2 = [ 1.1, 1.5, 2.1, 0.7]
    len1 = len(eigen1)
    len2 = len(eigen2)
    lda = len1 * len2
    vector = [ float(i) for i in range(1,1+lda) ]
    mat    = [ [float(i+j) for i in range(1,1+lda)] for j in range(1,1+lda) ]
    new_vec = hamiltonian_multiply_vec( eigen1, eigen2, mat, vector )

    print("Original Vector =", vector)
    print("Original Dipole Matrix =")
    for line in mat:
        print(line)

    mat_np = np.matrix(mat)
    vec_np = np.matrix(vector)
    vec_np = vec_np.H
    x = 0
    for i in range(len1):
        for j in range(len2):
            mat_np[x,x] += eigen1[i] + eigen2[j]
            x += 1
    print("REAL Hamiltonian =")
    for i in range(lda):
        for j in range(lda):
            print(mat_np[i,j], end="   ")
        print(" ")
    new_vec_np = mat_np*vec_np

    print("Compare Vectors (New vs Ori)")
    for i in range(len(new_vec)):
        print( new_vec[i]," <==> ",new_vec_np[i,0], "    diff =", abs(new_vec[i] - new_vec_np[i,0]) )









'''
def get_analytical_sln( oscillator_system_object, energy_upper_bound, 
                        unit_length, unit_mass, unit_hbar ):
    k_eff, states_configure, E_analytical, vectors = \
                analytical_solution.analytical_main( 
                                         oscillator_system_object.list_of_masses() ,
                                         oscillator_system_object.get_k_mat()      ,
                                         energy_upper_bound                        , 
                                         unit_length,
                                         unit_mass, 
                                         unit_hbar   )

    return np.sort( E_analytical )



def get_two_group_recouple_sln( recouple_hamiltonian_object ):
    dipole_mat        =  recouple_hamiltonian_object.get_pair_dipole_mat()[0]
    energy1, energy2  =  recouple_hamiltonian_object.get_sorted_eigen_energy()

    x = 0
    for i in range(len(energy1)):
        for j in range(len(energy2)):
            dipole_mat[ x ][ x ] += energy1[i] + energy2[j]
            x += 1

    energies, vectors = LA.eigh( dipole_mat, UPLO='U' )
    sorted_energies = np.sort( energies ) 
    return sorted_energies 







def driver_main( oscillator_system_object , recouple_hamiltonian_object, energy_upper_bound ):
    unit_length        = recouple_hamiltonian_object.unit_length
    unit_mass          = recouple_hamiltonian_object.unit_mass
    unit_hbar          = recouple_hamiltonian_object.unit_hbar    



    E_analytical = get_analytical_sln( oscillator_system_object, 
                                       energy_upper_bound, 
                                       unit_length,
                                       unit_mass,
                                       unit_hbar    )
                                                         

    E_recouple   = get_two_group_recouple_sln( recouple_hamiltonian_object )

    return np.sort(E_analytical) , np.sort(E_recouple) 
'''











