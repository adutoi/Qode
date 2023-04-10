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
import sys
import numpy as np
from qode.SCF_Code.matrix_operation import printnpline, cp_arbit_np_mat
import qode.SCF_Code.elec_fock_builder as e_fock

'''
convergence class, is initialized before SCF cycles. It took initial guess to
calculate first SCF energy.

In each cycle, __call__(***) is called to check SCF energy difference and the
products of corresponding guess vectors.

If both converged, __call__ return True and terminates the SCF loop.
'''
def get_low_evec(F):
    E,V = np.linalg.eigh(F)                  # diag(Fock) w/ evals sorted
    vec = np.matrix( np.zeros((len(F),1)) )  # vec to hold the lowest vector
    for i in range(len(F)):
        vec[i,0] = V[i,0]    
    return vec, E[0]

def replace_vec(guess, index, vector):
    for i in range(len(guess)):
        guess[i,index] = vector[i,0]

class convergence(object):
    def __init__(self, num_occ, thresh):        
        self.num_occ = num_occ
        self.thresh = thresh
        self.dS = np.zeros((self.num_occ))
        self.dE = np.zeros((self.num_occ))
        self.count = 0
        
    def __call__(self, old_vec, new_vec, old_E, new_E):
        k = self.count % self.num_occ
        self.dS[k] = e_fock.col_dot(old_vec, 0, new_vec, 0)
        self.dE[k] = abs( old_E - new_E )
        # check vector
        vec_conv = 0
        for i in self.dS:
            if abs(abs(i) - 1.0) < self.thresh:
                vec_conv += 1
        
        # check energy
        E_conv = 0
        for i in self.dE:
            if abs(i/old_E) < self.thresh:
                E_conv += 1
        print("conv: vec =",vec_conv,"energy =",E_conv)
        if vec_conv == self.num_occ and E_conv == self.num_occ:
            return True
        else:
            self.count += 1
        


def SCF(fock, guess, thresh):
    converged = False
    scf_count = 0
    conv = convergence(fock.num_e, thresh)
    
    while not converged:
        i = 0
        while i < fock.num_e and (not converged):
            j = fock.occ[i]         # To get the index
            scf_count += 1
            print('Cycle',scf_count)
            #print(i,"in guess",j,"in Matrix")
            F = fock(guess, j)
            old_vec = np.matrix( np.zeros((fock.n,1)) )
            e_fock.cp_col(old_vec, 0, guess, i)
            old_E = old_vec.H * F * old_vec
            new_vec, new_E = get_low_evec(F)
            replace_vec(guess, i, new_vec)
            converged = conv(old_vec, new_vec, old_E, new_E)
            i += 1
            
    return guess














