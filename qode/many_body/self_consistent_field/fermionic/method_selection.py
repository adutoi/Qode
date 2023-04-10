from . import RoothanHall
from .mean_field import HartreeFock

def RHF_RoothanHall_Orthonormal(n_elec, integrals, C=None, damp=(5, 0.2), thresh=1e-8):
    if (n_elec % 2) != 0:  raise ValueError
    h, V = integrals
    return RoothanHall.compute(n_elec//2, RoothanHall.orthonormal_diagonalizer(), HartreeFock.restricted(h,V,damp), C, thresh)

def RHF_RoothanHall_Nonorthogonal(n_elec, integrals, C=None, damp=(5, 0.2), thresh=1e-8):
    if (n_elec % 2) != 0:  raise ValueError
    S, h, V = integrals
    return RoothanHall.compute(n_elec//2, RoothanHall.nonorthogonal_diagonalizer(S), HartreeFock.restricted(h,V,damp), C, thresh)

def RHF_RoothanHall_BiorthogonalS(n_elec, integrals, C=None, damp=(5, 0.2), thresh=1e-8):
    if (n_elec % 2) != 0:  raise ValueError
    S, h, V = integrals
    return RoothanHall.compute(n_elec//2, RoothanHall.biorthogonal_diagonalizer_Sdressed(S), HartreeFock.restricted(h,V,damp), C, thresh)

def RHF_RoothanHall_Biorthogonal(n_elec, integrals, C=None, damp=(5, 0.2), thresh=1e-8):
    if (n_elec % 2) != 0:  raise ValueError
    h, V = integrals
    return RoothanHall.compute(n_elec//2, RoothanHall.biorthogonal_diagonalizer_general(), HartreeFock.restricted(h,V,damp), C, thresh)

def UHF_RoothanHall_Orthonormal(n_elec, integrals, C=None, damp=(5, 0.2), thresh=1e-8):
    h, V = integrals
    return RoothanHall.compute(n_elec, RoothanHall.orthonormal_diagonalizer(), HartreeFock.unrestricted(h,V,damp), C, thresh)

def UHF_RoothanHall_Nonorthogonal(n_elec, integrals, C=None, damp=(5, 0.2), thresh=1e-8):
    S, h, V = integrals
    return RoothanHall.compute(n_elec, RoothanHall.nonorthogonal_diagonalizer(S), HartreeFock.unrestricted(h,V,damp), C, thresh)

def UHF_RoothanHall_BiorthogonalS(n_elec, integrals, C=None, damp=(5, 0.2), thresh=1e-8):
    S, h, V = integrals
    return RoothanHall.compute(n_elec, RoothanHall.biorthogonal_diagonalizer_Sdressed(S), HartreeFock.unrestricted(h,V,damp), C, thresh)

def UHF_RoothanHall_Biorthogonal(n_elec, integrals, C=None, damp=(5, 0.2), thresh=1e-8):
    h, V = integrals
    return RoothanHall.compute(n_elec, RoothanHall.biorthogonal_diagonalizer_general(), HartreeFock.unrestricted(h,V,damp), C, thresh)
