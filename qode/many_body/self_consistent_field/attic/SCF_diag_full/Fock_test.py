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
from qode.SCF_diag_full.matrix_operation import printline,printnpline
from math import exp,pi

def get_wf(col,_list_,ng):
    wavef = []
    count = 0
    for i in range(len(ng)):
        for j in range(ng[i]):
            wavef += [[_list_[count][3],_list_[count][4],_list_[count][2],\
                       col[i]]]
            count += 1
    return wavef



def get_column(U,col_num):
    'col_num: starts with one.' 
    cn = col_num - 1
    size = len(U)
    c_wf = []
    for i in range(size):
        c_wf += [U[i,cn]]
    return c_wf

def NFactor(expn):
    'Normalization Factor'
    return pow(2.0*expn/pi,0.75)

def comp_wf(x,wf):
    'Equation: c * N * Uji * exp( -expn[x-z]^2 )'
    func = 0.0
    for i in wf:
        func += i[1]*NFactor(i[0])*i[3]*exp(-i[0]*pow(x-i[2],2))
    return func

def plot(wf,xrange,dx):
    new_plot = ''
    x = -abs(xrange)
    while x <= abs(xrange):
        tmp = comp_wf(x,wf)
        new_plot += str(x)+'\t'+str(tmp)+'\n'
        x += dx
    
    return new_plot

#============ Main Fuction ============
def fock_main(U_on,T_t,N_t,V_t,_list_,ng):
    U0 = U_on
    F0 = T_t + N_t
    E1,U1 = np.linalg.eig(F0)
    print(E1)
    U_wf = U0 * U1   # U_wf: Unitary for WaveFunction

    col1 = get_column(U_wf,1)
    #print(col1)
    col2 = get_column(U_wf,2)
    col3 = get_column(U_wf,3)
    col4 = get_column(U_wf,4)

    ############### CHECK HERE  ####################
    ################################################
    xrange = 5.0
    dx = 0.1
    ct = 1
    for n in [col1,col2,col3,col4]:
        wf = get_wf(n,_list_,ng)
        to_write = plot(wf,xrange,dx)

        file = 'wavef'+str(ct)+'.txt'
        f = open(file,'w')
        f.write(to_write)
        f.close()

        ct += 1

