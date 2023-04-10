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
# RIGHT NOW ONLY FOR Ground State MP2 Correlation Energy
from qode.SCF_diag_full import eigensort
import numpy as np


def mp2_E(n,a,b,r,s,V,E):
    V1 = V[a*n+b,r*n+s]
    V2 = V[a*n+b,s*n+r]
    dE = E[a] + E[b] - E[r] - E[s]
    return 0.5*(V1*V1-V1*V2)/dE

def mp2_main(EHF,E,V,alpha_e,beta_e):
    n = len(E)
    hn = n//2
    Ia,Ib = eigensort.sort_main(E)
    #print(Ia,Ib)
    #aE = np.zeros(n//2)
    #bE = np.zeros(n//2)
    #print(aE,bE)
    #print("HF ENERGY:",EHF)
    #print("EigenValues:\n",E)
    #print("Repulsion Mat:\n",V)

    aaaaE = 0.0
    bbbbE = 0.0
    ababE = 0.0
    
    if (alpha_e + beta_e) > 1:  # CLOSE SHELL
        # <aa|aa> Part:
        for i in range(alpha_e):
            a = Ia[i]
            for j in range(alpha_e):
                b = Ia[j]
                for k in range(alpha_e,hn):
                    r = Ia[k]
                    for l in range(alpha_e,hn):
                        s = Ia[l]
                        aaaaE += mp2_E(n,a,b,r,s,V,E)
        # <bb|bb> Part:
        for i in range(beta_e):
            a = Ib[i]
            for j in range(beta_e):
                b = Ib[j]
                for k in range(beta_e,hn):
                    r = Ib[k]
                    for l in range(beta_e,hn):
                        s = Ib[l]
                        bbbbE += mp2_E(n,a,b,r,s,V,E)
        # <ab|ab> Part:
        for i in range(alpha_e):
            a = Ia[i]
            for j in range(beta_e):
                b = Ib[j]
                for k in range(alpha_e,hn):
                    r = Ia[k]
                    for l in range(beta_e,hn):
                        s = Ib[l]
                        ababE += mp2_E(n,a,b,r,s,V,E)
        # <ba|ba> Part:
        for i in range(beta_e):
            a = Ib[i]
            for j in range(alpha_e):
                b = Ia[j]
                for k in range(beta_e,hn):
                    r = Ib[k]
                    for l in range(alpha_e,hn):
                        s = Ia[l]
                        ababE += mp2_E(n,a,b,r,s,V,E)

    else:
        pass

    EMP2 = aaaaE + bbbbE + ababE
    print("\nMP2 Correlation:\n")
    print("aaaa correlation E =",aaaaE)
    print("bbbb correlation E =",bbbbE)
    print("abab correlation E =",ababE)
    print("MP2 CORRELATION ENERGY =",EMP2)
    print("Total MP2 Energy =",EHF + EMP2)
    return EMP2



