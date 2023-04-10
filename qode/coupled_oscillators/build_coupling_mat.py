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

def correlate( 
               mat_list,
               coupling
	      ):
    "Take any numbers of matrices and their coupling matrix to build big coupling matrix"
    # total LD and build big matrix
    ld_new = 0
    for mat in mat_list:
        ld_new += len(mat)
    new_mat = [ [ 0.0 for i in range(ld_new) ] for j in range(ld_new) ]
    # Fill in diagonal matrices
    shift = 0
    for mat in mat_list:
        loop_range = len(mat)
        for x in range( loop_range ):
            for y in range( loop_range ):
                new_mat[ x+shift ][ y+shift ] = mat[ x ][ y ]
        shift += loop_range
    
    # Fill in coupling blocks
    x_shift = 0
    loop_range = len(coupling)
    for i in range( loop_range ):         # looping over i for kij
        y_shift = 0                          # get y_shift 
        for k in range(i+1):                 # get y_shift
            y_shift += len( mat_list[ k ] )  # get y_shift
        for j in range( i+1, loop_range ):   # looping over j for kij
            x_range = len(mat_list[i])          # find mat[i] size
            y_range = len(mat_list[j])          # find mat[j] size
            for x in range(x_range):            # fill in kij into that block( loop over x )
                for y in range(y_range):        # fill in kij into that block( loop over y ) 
                    new_mat[ x + x_shift  ][ y + y_shift  ] = coupling[ i ][ j ]
            y_shift +=  y_range 
        x_shift += len(mat_list[i])
        
    return new_mat


if __name__ == '__main__':
    from qode.util.printfunc import printline
    m = 1
    n = 1
    mat1 = [[ 10.0*(i+1) for  i in range(m)] for j in range(m) ]
    mat2 = [[ 20.0*(i+1) for  i in range(n)] for j in range(n) ]
    mat3 = [[ 30.0*(i+1) for  i in range(m)] for j in range(m) ]
    mat4 = [[ 40.0*(i+1) for  i in range(m)] for j in range(m) ]
    printline(mat1,"    ")
    printline(mat2,"    ")
    printline(mat3,"    ")
    printline(mat4,"    ")
    
    cp2 = [[ i+j+1.0 for i in range(2) ] for j in range(2) ]
    cp3 = [[ i+j+1.0 for i in range(3) ] for j in range(3) ]
    cp4 = [[ i+j+1.0 for i in range(4) ] for j in range(4) ]
    
    print("COUPLING MAT FOR TWO")
    printline(cp2,"    ")
    print("COUPLING MAT FOR THREE")
    printline(cp3,"    ")
    print("COUPLING MAT FOR FOUR")
    printline(cp4,"    ")
    
    x = correlate( [mat1,mat2], cp2 )
    print("TWO BODY COUPLING TEST")
    printline(x,"  ")
    
    y = correlate( [mat1,mat2,mat2], cp3 )
    print("THREE BODY COUPLING TEST")
    printline(y,"  ")
    
    z = correlate( [mat4,mat1,mat2,mat2], cp4 )
    print("FOUR BODY COUPLING TEST")
    printline(z,"  ")

    sing_num_mat = [ [ [i * 111.] ] for i in range(1, 5) ]
    print(sing_num_mat)

    test = correlate( sing_num_mat, cp4 )
    printline(test,"    ")





