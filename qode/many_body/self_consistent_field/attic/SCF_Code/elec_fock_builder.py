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
import copy
from qode.SCF_Code.matrix_operation import printnpline, cp_arbit_np_mat
from qode.SCF_Code import elec_mat2fock
from qode.SCF_Code import elec_HF_energy


# A big number to fill out diagonal elements in Fock matrix that is
# not of interest.
big_number = 1.0e6



def col_dot(M1, col1, M2, col2):
    'Return the product of column 1 in M1 and column 2 in M2'
    product = 0.0
    for i in range(len(M1)):
        product += M1[i,col1] * M2[i,col2]
    return product

def cp_col(M1, col1, M2, col2):
    'copy column 2 in M2 to column 1 in M1'
    for i in range(len(M1)):
        M1[i,col1] = M2[i,col2]

def reduce_fock(F, i, occ, alpha_e, beta_e):
    "Wipe Out OTHER occupies' matrix elements before diagonalization"
    n = len(F)
    # Wipe Out other occupies.
    if i in occ:
        for k in occ:
            if k != i:              # Not the one that is interested
                for x in range(n):  # Then wipe out
                    F[x,k] = 0.0
                    F[k,x] = 0.0
                F[k,k] = big_number
    '''
    print("Reduced Fock")
    printnpline(F)'''
    return F

def build_basis(guess, ae, be):
    'build full vector (unitary) matrix from occupies'
    n = len(guess)
    spat_n = n//2   # number of spatial orbitals
    

    I1 = np.matrix( np.zeros((n, spat_n)) ) # alpha Identity part
    I2 = np.matrix( np.zeros((n, spat_n)) ) # beta  Identity part    
    # copy Identity into I1 and I2
    for i in range(spat_n):
        I1[i,i] = 1.0
        I2[i+spat_n,i] = 1.0
    '''
    print("I1")
    print(I1)
    print("I2")
    print(I2)
    '''
    # Subtract Guess from I1 and I2
    # |i1> = |i1> - |j><j|i1>; |i2> = |i2> - |j><j|i2>
    for i in range(spat_n):
        for j in range(ae + be):
            pro1 = col_dot(I1, i, guess, j) # pro1 = <j|i1>
            pro2 = col_dot(I2, i, guess, j) # pro2 = <j|i2>
            for l in range(n):
                I1[l,i] -= guess[l,j] * pro1     # |i1> = |i1> - |j> * pro1
                I2[l,i] -= guess[l,j] * pro2     # |i2> = |i2> - |j> * pro2
                
    J1 = cp_arbit_np_mat(I1)
    J2 = cp_arbit_np_mat(I2)

    '''
    print("J1")
    printnpline(J1)
    print("J2")
    printnpline(J2)
    '''
    
    # Overlap Matrix for I1 and I2
    S1 = np.matrix( np.zeros((spat_n, spat_n)) )
    S2 = np.matrix( np.zeros((spat_n, spat_n)) )
    
    for i in range(spat_n):
        for j in range(spat_n):
            S1[i,j] = col_dot(J1, i, J1, j) # column vector products
            S2[i,j] = col_dot(J2, i, J2, j) # column vector products
    
    '''
    print("S1")
    printnpline(S1)
    print("S2")
    printnpline(S2)
    '''
    
    # diag(overlap matrix) for S1 and S2
    # eigh(Matix) sorts eigenvalues so here it is used,
    # all the virtual eigenvectors have higher energy on the right side of V1, V2
    E1, V1 = np.linalg.eigh(S1)
    E2, V2 = np.linalg.eigh(S2)
    '''
    print("E1",E1)
    print("E2",E2)
    print("V1")
    printnpline(V1)
    print("V2")
    printnpline(V2)'''
    
    # L1 = J1 * K1; L2 = J2 * K2
    L1 = J1 * V1
    L2 = J2 * V2
    '''
    print("L1")
    printnpline(L1)
    print("L2")
    printnpline(L2)'''

    
    # Combine guess, alpha and beta to the whole U matrix:
    U = np.matrix( np.zeros((n,n)) )
    
    for i in range(spat_n):    
        cp_col(U, i, L1, i)
        cp_col(U, i + spat_n, L2, i)

    # Copy the guess into self.U:
    for i in range(ae):
        cp_col(U, i, guess, i)
    for i in range(be):
        cp_col(U, i + spat_n, guess, i + ae)
    '''
    print("U")
    printnpline(U)'''
    return U




class fock_builder(object):
    'Build Fock Matrices with Spin Separated'
    def __init__(self, T_sp, N_sp, V_sp ,alpha_e, beta_e, NRE):
        self.n = len(T_sp)
        self.spat_n = self.n//2    # number of spatial orbitals
        self.n_elec = alpha_e + beta_e
        self.T0 = cp_arbit_np_mat(T_sp)
        self.N0 = cp_arbit_np_mat(N_sp)
        self.V0 = cp_arbit_np_mat(V_sp)
        self.ae = alpha_e
        self.be = beta_e
        self.num_e = alpha_e + beta_e
        self.NRE = NRE
        self.occa, self.occb = self.get_scf_index()
        self.occ = self.occa + self.occb
        print("Occupied Fock Vector Index =",self.occ)
        
        
    def get_scf_index(self):
        "Sort eigenvalues of diagonalizing initial guess, find the indecies of occupied orbitals"
        occa = []
        occb = []
        for i in range(self.ae):
            occa += [i]
        for i in range(self.spat_n, self.spat_n+self.be):
            occb += [i]
        print("occupied Alpha electron Index =",occa)
        print("occupied Beta  electron Index =",occb)
        return occa, occb

        
    def __call__(self, guess, i):
        "input guess and orbital index i,solve for virtuals, build and reduce new Fock, diag(F) and update guess"
        self.guess = cp_arbit_np_mat(guess)
        # solve for virtuals to build basis Ufull:
        self.U = build_basis(self.guess, self.ae, self.be)
        
        # Build new Fock matrix
        self.F1, self.T1, self.N1, self.V1 = elec_mat2fock.build_Fock(\
            self.U, self.T0, self.N0, self.V0, self.ae, self.be)

        # Reduce Fock matrix to remove other occupied orbitals' information
        self.F1 = reduce_fock(self.F1, i, self.occ, self.ae, self.be)
        self.F1 = self.U * self.F1 * (self.U.H)
        return self.F1
    





