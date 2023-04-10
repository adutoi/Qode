#    (C) Copyright 2019 Anthony D. Dutoi
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
import numpy



class operator(object):
    def __init__(self, order, n_crop, n_dsop, empty=False):
        self.order   = order
        self.dims    = (n_crop, n_dsop)
        self.scalars = [0.]	# The first element is not an array but a scalar
        dims = self.dims
        for o in range(1,order+1):
            if empty:  self.scalars += [None]
            else:      self.scalars += [numpy.zeros(dims)]
            dims = (n_crop,) + dims + (n_dsop,)



def hamiltonian(h, V, occupied):	# V in raw spin-orb rep, physicists ordering
    n_orb = h.shape[0]
    n_occ = len(occupied)
    n_vrt = n_orb - n_occ

    virtuals = [a for a in range(n_orb) if a not in occupied]
    orbitals = occupied + virtuals

    H = operator(2, n_orb, n_orb)

    for i in occupied:
        H.scalars[0] += h[i,i]
        for j in occupied:
            H.scalars[0] += (1/2.) * (V[i,j,i,j] - V[i,j,j,i])

    for pp,p in enumerate(orbitals):
        for qq,q in enumerate(orbitals):
            H.scalars[1][pp,qq] = h[p,q]
            for i in occupied:
                H.scalars[1][pp,qq] += V[p,i,q,i] - V[p,i,i,q]

    for pp,p in enumerate(orbitals):
        for qq,q in enumerate(orbitals):
            for rr,r in enumerate(orbitals):
                for ss,s in enumerate(orbitals):
                    H.scalars[2][pp,qq,rr,ss] += V[p,q,r,s] - V[p,q,s,r]

    E0 = H.scalars[0]
    H.scalars[0] = 0
    return E0, H



def dot(t1,t2):
    order = min(t1.order, t2.order)
    result = t1.scalars[0] * t2.scalars[0]	# order should not be less than zero
    denom = 1.
    axes = [0,1]
    Daxes = [2,3]
    for o in range(1,order+1):
        denom *= (o**2)
        result += numpy.tensordot(t1.scalars[o], t2.scalars[o], axes=(axes,axes)) / denom
        axes += Daxes
        Daxes = [Daxes[0]+2, Daxes[1]+2]
    return result

def scale(c, t):
    for o in range(t.order+1):
        t.scalars[o] *= c

def increment(t,Dt,c=1):
    for o in range(Dt.order+1):
        t.scalars[o] += c * Dt.scalars[o]

def copy(t):
    order = t.order
    n_crop, n_dsop = t.dims
    tNew = operator(order, n_crop, n_dsop, empty=True)
    for o in range(t.order+1):
        tNew.scalars[o] = numpy.copy(t.scalars[o])
    return tNew


def apply_denominators(energies, omega):
    e_occ, e_vrt = energies
    if omega.order>2:  raise RuntimeError("no denominators yet for general orders")
    result = copy(omega)
    result.scalars[0] = 0.
    for a,ea in enumerate(e_vrt):
        for i,ei in enumerate(e_occ):
            result.scalars[1][a,i] /= (ea-ei)
    for a,ea in enumerate(e_vrt):
        for b,eb in enumerate(e_vrt):
            for i,ei in enumerate(e_occ):
                for j,ej in enumerate(e_occ):
                    result.scalars[2][a,b,i,j] /= (ea+eb-ei-ej)
    return result

def check_member(v):  pass
def check_lin_op(op): return False
