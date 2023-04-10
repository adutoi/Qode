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
#
#   THIS MODULE TAKES recsys_HO_hamiltonian class TO BUILD three-body Hamiltonian 
#   FROM paired transition dipole matrices 
#
#





def two_to_three_main( recsys_hamil_obj ):
    # Unpack three eigenvalue 1D arrays
    assert( len( recsys_hamil_obj.get_sorted_eigen_energy() ) == 3 )
    eigenenergy1 , eigenenergy2 , eigenenergy3 = \
                   recsys_hamil_obj.get_sorted_eigen_energy()

    # Unpack Three 2D Transition Dipole Matrices.
    assert( len( recsys_hamil_obj.get_pair_dipole_mat() ) == 3 )
    dipole_mat1 , dipole_mat2 , dipole_mat3 = \
                   recsys_hamil_obj.get_pair_dipole_mat()

    # Get first dimensions of the three
    lda1 = len( eigenenergy1 )
    lda2 = len( eigenenergy2 )
    lda3 = len( eigenenergy3 )

    print("Start Building Three-Body Hamiltonian")
    print("dimensions = ", lda1*lda2*lda3 , "*" , lda1*lda2*lda3)

    # build THREE-BODY HAMILTONIAN:
    three_body_hamiltonian = [ [ 0.0 for i in range( lda1*lda2*lda3 ) ] \
                                      for j in range( lda1*lda2*lda3 ) ]

    # Generate 6-fold nested loops to build the THREE-BODY MATRIX:
    # first  dimension index: i, k, m
    # second dimension index: j, l, n
    # first  index = i*lda2*lda3 + k*lda3 + m
    # second index = j*lda2*lda3 + l*lda3 + n
    for i in range( lda1 ):
        for k in range( lda2 ):
            for m in range( lda3 ):
                # Do something
                for j in range( lda1 ):
                    for l in range( lda2 ):
                        for n in range( lda3 ):
                           if ( m == n ):
                               three_body_hamiltonian[ i*lda2*lda3 + k*lda3 + m ][ j*lda2*lda3 + l*lda3 + n ] +=\
                                     dipole_mat1[i*lda2 + k][j*lda2 + l]
                           elif ( k == l ):
                               three_body_hamiltonian[ i*lda2*lda3 + k*lda3 + m ][ j*lda2*lda3 + l*lda3 + n ] +=\
                                     dipole_mat2[i*lda3 + m][j*lda3 + n]
                           elif ( i == j ):
                               three_body_hamiltonian[ i*lda2*lda3 + k*lda3 + m ][ j*lda2*lda3 + l*lda3 + n ] +=\
                                     dipole_mat3[k*lda3 + m][l*lda3 + n]

    for i in range( lda1 ):
        for k in range( lda2 ):
            for m in range( lda3 ):
                three_body_hamiltonian[ i*lda2*lda3 + k*lda3 + m ][ i*lda2*lda3 + k*lda3 + m ] += \
                        eigenenergy1[ i ] + eigenenergy2[ k ] + eigenenergy3[ m ]


    return three_body_hamiltonian





