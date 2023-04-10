#    (C) Copyright 2016 Yuhong Liu
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
import numpy as np
from qode.SCF_Code import repul_transform
from qode.SCF_Code import sum_repulsion
from qode.SCF_Code.matrix_operation import cp_arbit_np_mat

# build Fock matrix for electrons
def build_Fock(U1,T0,N0,V0,alpha_e,beta_e):
    '''
    E1: eigenvalues; U1:eigenvectors; T0,N0,V0: spin matrices for
    Kin, Nucl Attr, and Elec Repul Energy;
    return new Fock matrix: F1; transformed matrices: T1,N1,V1.
    '''
    
    U1 = cp_arbit_np_mat(U1)
    T0 = cp_arbit_np_mat(T0)
    N0 = cp_arbit_np_mat(N0)

    T1 = U1.H * T0 * U1
    N1 = U1.H * N0 * U1
    
    if (alpha_e + beta_e) > 1:
        V0 = cp_arbit_np_mat(V0)
        V1 = repul_transform.transform_V(V0,U1)
        vHF = sum_repulsion.sum_repul(V1,alpha_e,beta_e)
        F1 = T1 + N1 + vHF 
    else:
        F1 = T1 + N1
        V1 = V0

    return F1,T1,N1,V1
