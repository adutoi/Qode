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
import copy
from math import factorial


def combination(N,M):
    return factorial(N) // factorial(M) // factorial(N-M)



class cisd_amplitude(object):
    def __init__(self,
                 num_alpha_elec,
                 num_beta_elec,
                 num_basis_func
                 ):
        # num_basis_func is the number of spatial orbitals
        self.num_alpha_elec = num_alpha_elec
        self.num_beta_elec  = num_beta_elec
        self.num_basis_func = num_basis_func
        self.num_alpha_virt_orb = self.num_basis_func - self.num_alpha_elec
        self.num_beta_virt_orb  = self.num_basis_func - self.num_beta_elec
        
        # CHECK FIRST IF ALL KINDS OF DOUBLE EXCITATIONS ARE POSSIBLE!!
        if self.num_alpha_virt_orb >= 0 and self.num_beta_virt_orb >= 0:
            pass
        else:
            print("Something really weird went wrong")
            raise RuntimeError
                
        # build three amplitude holders, Ref is for one number,
        # Singles is a matrix completely filled with numbers,
        # Doubles is a list of upper submatrices ( Please refer to the png file in current folder )
        
        # Reference Amplitude:
        self.ref_amplitude = np.matrix(np.zeros((1,1)))
        self.ref_amplitude[0,0] = 1.0
    
        # Singles Amplitudes: ( ROWS: m_alpha , m_beta ; COLUMNS: e_alpha , e_beta )
        self.single_amp_mat = np.matrix(np.zeros(( self.num_alpha_elec + self.num_beta_elec , self.num_alpha_virt_orb \
                                        + self.num_beta_virt_orb )) )
        
        # Doubles Amplitudes: ( ROWS: e_alpha, e_beta; COLUMNS: f_alpha, f_beta )
        # Just make a super buge matrix to hold zeros.
        total_ld = ( self.num_alpha_elec + self.num_beta_elec ) * ( self.num_alpha_virt_orb + self.num_beta_virt_orb )
        self.double_amp_mat = np.matrix( np.zeros(( total_ld, total_ld )) )
    
    # Functions to retrieve electronic configurations
    def get_num_alpha_elec(self):
        return self.num_alpha_elec
    
    def get_num_beta_elec(self):
        return self.num_beta_elec
    
    def get_num_spatial_orb(self):
        return self.num_basis_func
    
    def get_num_alpha_virt_orb(self):
        return self.num_alpha_virt_orb
    
    def get_num_bete_virt_orb(self):
        return self.num_beta_virt_orb
    
    # Functions to retrieve CISD amplitudes
    def get_ref_amplitude(self):
        return self.ref_amplitude[0,0]
           
    def get_single_amplitude(self):
        return self.single_amp_mat
    
    def get_double_amplitude(self):
        return self.double_amp_mat

    def get_num_singles(self):
        num_occ = self.get_num_alpha_elec() + self.get_num_beta_elec()
        num_vrt = self.get_num_spatial_orb() * 2 - num_occ
        return num_occ*num_vrt

    def get_num_doubles(self):
        num_occ = self.get_num_alpha_elec() + self.get_num_beta_elec()
        num_vrt = self.get_num_spatial_orb() * 2 - num_occ
        return combination(num_occ,2) * combination(num_vrt,2)

    # Functions to update CISD amplitudes
    def update_ref_amplitude(self, a_number ):
        self.ref_amplitude[0,0] = a_number
        
    def update_single_amplitude(self, a_np_matrix ):
        if self.single_amp_mat.shape == a_np_matrix.shape:
            self.single_amp_mat = np.matrix( copy.deepcopy( a_np_matrix ) )
        else:
            print("Incompatible Matrix Size assigned to Single Excitation Amplitude Matrix.")
            sys.exit(1)
    
    def update_double_amplitude(self, a_np_matrix ):
        if self.double_amp_mat.shape == a_np_matrix.shape:
            self.double_amp_mat = np.matrix( copy.deepcopy( a_np_matrix ) )
        else:
            print("Incompatible Matrix Size assigned to Double Excitation Amplitude Matrix.")
            sys.exit(1)

    def get_vec_dimension(self):
        return 1 + self.get_num_singles() + self.get_num_doubles()

    def amp_obj_to_long_vec(self):
        num_occ = self.get_num_alpha_elec() + self.get_num_beta_elec()
        num_vrt = self.get_num_spatial_orb() * 2 - num_occ
        long_vec = [ self.ref_amplitude[0,0] ]
        for i in range(num_occ):
            for a in range(num_vrt):
                long_vec += [ self.single_amp_mat[i,a] ]
        for i in range(num_occ):
            for j in range(i+1, num_occ):
                for a in range(num_vrt):
                    for b in range(a+1, num_vrt):
                        long_vec += [self.double_amp_mat[i*num_vrt+a, j*num_vrt+b] ]
        size = len(long_vec)
        return copy.deepcopy( np.array(long_vec).reshape((size,1)) )

    def update_amp_from_long_vec(self, long_vec):
        num_occ = self.get_num_alpha_elec() + self.get_num_beta_elec()
        num_vrt = self.get_num_spatial_orb() * 2 - num_occ
        idx = 0
        self.ref_amplitude[0,0] = long_vec[idx,0]
        idx += 1
        for i in range(num_occ):
            for a in range(num_vrt):
                self.single_amp_mat[i,a] = long_vec[idx,0]
                idx += 1
        for i in range(num_occ):
            for j in range(i+1, num_occ):
                for a in range(num_vrt):
                    for b in range(a+1, num_vrt):
                        self.double_amp_mat[i*num_vrt+a, j*num_vrt+b] = long_vec[idx,0]
                        idx += 1
        if idx != long_vec.shape[0]:
            raise ValueError
        

    def clean_all_amplitude(self):
        self.ref_amplitude[0,0] = 0.0
        self.single_amp_mat.fill(0.0)
        self.double_amp_mat.fill(0.0)

    def return_info(self):
        output_buf  = "CISD Vector Object\n"
        output_buf += "-----------------------------------------------------\n"
        output_buf += str(self.get_ref_amplitude()) + '\n'
        output_buf += "-----------------------------------------------------\n"
        output_buf += str(self.get_single_amplitude()) + '\n'
        output_buf += "-----------------------------------------------------\n"
        output_buf += str(self.get_double_amplitude()) + '\n'
        output_buf += "-----------------------------------------------------\n"
        return output_buf

    def print_info(self):
        print(self.return_info())

def add_and_update_cisd_amp_obj( destination_amp_obj, incremental_amp_obj  ):
    # They are ALL NUMPY 2D ARRAYS. Straight Addition Is OK!
    destination_amp_obj.update_ref_amplitude( destination_amp_obj.get_ref_amplitude() + incremental_amp_obj.get_ref_amplitude() )
    destination_amp_obj.update_single_amplitude( destination_amp_obj.get_single_amplitude() + incremental_amp_obj.get_single_amplitude() )
    destination_amp_obj.update_double_amplitude( destination_amp_obj.get_double_amplitude() + incremental_amp_obj.get_double_amplitude() )
    return destination_amp_obj


if __name__ == "__main__":
    cisd_vec = cisd_amplitude(1,1,4)
    cisd_vec.print_info()

    