/*   (C) Copyright 2019 Anthony D. Dutoi and Yuhong Liu
 *
 *   This file is part of Qode.
 *
 *   Qode is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Qode is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Qode.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdlib.h>		// malloc pointer arrays and dimensions
#include <math.h>		// fabs
#include "PyC_types.h"
#include "partition_storage.h"	// Defines ptrs1 and ptrs2

// only works if output tensors are zero upon entry ... H2 needs to be a flat list with dummy arrays in the lower triangle

void load(PyInt Nfrag, BigInt* Nvrt, PyFloat thresh,
          Double** H1, Double** H2,
          Double* Hc,
          Double* HEx, Double* HFo, Double* HFv, Double* HDx,
          Double* HExEx, Double* HExFo, Double* HExFv, Double* HExDx, Double* HFoFo, Double* HFoFv, Double* HFoDx, Double* HFvFv, Double* HFvDx, Double* HDxDx)
    {
    Double tmp;
    Double** H    = (Double**)malloc(sizeof(Double*)*Nfrag*Nfrag);
    BigInt*  Ntot = (BigInt*)malloc(sizeof(Int)*Nfrag);
    for (int M=0;  M<Nfrag;  M++) {Ntot[M] = Nvrt[M]+1;}

    // constant
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        Hc[0] += H1[M1][0*Ntot[M1] + 0];
        for (int M2=0;  M2<Nfrag;  M2++)
            {
            if (M1<M2) {Hc[0] += (1/2.) * H2[M1*Nfrag + M2][0*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + 0];}
            if (M1>M2) {Hc[0] += (1/2.) * H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + 0*Ntot[M2]*Ntot[M1] + 0*Ntot[M1] + 0];}
            }
        }

    // Ex
    ptrs1(H, HEx, Nfrag, Nvrt, 1,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int u=0;  u<Nvrt[M1];  u++)
            {
            tmp = H1[M1][(u+1)*Ntot[M1] + 0];
            for (int M2=0;  M2<Nfrag;  M2++)
                {
                if (M1<M2) {tmp += H2[M1*Nfrag + M2][(u+1)*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + 0];}
                if (M1>M2) {tmp += H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M2]*Ntot[M1] + 0*Ntot[M1] + 0];}
                }
            if (fabs(tmp)>thresh) {H[M1][u] = tmp;}
            }
        }

    // Fo
    ptrs1(H, HFo, Nfrag, Nvrt, 0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        tmp = -H1[M1][0*Ntot[M1] + 0];
        for (int M2=0;  M2<Nfrag;  M2++)
            {
            if (M1<M2) {tmp -= H2[M1*Nfrag + M2][0*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + 0];}
            if (M1>M2) {tmp -= H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + 0*Ntot[M2]*Ntot[M1] + 0*Ntot[M1] + 0];}
            }
        if (fabs(tmp)>thresh) {H[M1][0] = tmp;}
        }

    // Fv
    ptrs1(H, HFv, Nfrag, Nvrt, 1,1);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int u=0;  u<Nvrt[M1];  u++)
            {
            for (int v=0;  v<Nvrt[M1];  v++)
                {
                tmp = H1[M1][(u+1)*Ntot[M1] + (v+1)];
                for (int M2=0;  M2<Nfrag;  M2++)
                    {
                    if (M1<M2) {tmp += H2[M1*Nfrag + M2][(u+1)*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + (v+1)*Ntot[M2] + 0];}
                    if (M1>M2) {tmp += H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M2]*Ntot[M1] + 0*Ntot[M1] + (v+1)];}
                    }
                if (fabs(tmp)>thresh) {H[M1][u*Nvrt[M1] + v] = tmp;}
                }
            }
        }

    // Dx
    ptrs1(H, HDx, Nfrag, Nvrt, 1,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int u=0;  u<Nvrt[M1];  u++)
            {
            tmp = H1[M1][0*Ntot[M1] + (u+1)];
            for (int M2=0;  M2<Nfrag;  M2++)
                {
                if (M1<M2) {tmp += H2[M1*Nfrag + M2][0*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + (u+1)*Ntot[M2] + 0];}
                if (M1>M2) {tmp += H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + 0*Ntot[M2]*Ntot[M1] + 0*Ntot[M1] + (u+1)];}
                }
            if (fabs(tmp)>thresh) {H[M1][u] = tmp;}
            }
        }

    // ExEx
    ptrs2(H, HExEx, Nfrag, Nvrt, 1,2,0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    tmp = (1/2.) * H2[M2*Nfrag + M1][(v+1)*Ntot[M1]*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M2]*Ntot[M1] + 0*Ntot[M1] + 0];
                    if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2] + v] = tmp;}
                    }
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    tmp = (1/2.) * H2[M1*Nfrag + M2][(u+1)*Ntot[M2]*Ntot[M1]*Ntot[M2] + (v+1)*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + 0];
                    if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2] + v] = tmp;}
                    }
                }
            }
        }

    // ExFo
    ptrs2(H, HExFo, Nfrag, Nvrt, 1,0,0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                tmp = -H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M2]*Ntot[M1] + 0*Ntot[M1] + 0];
                if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u] = tmp;}
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                tmp = -H2[M1*Nfrag + M2][(u+1)*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + 0];
                if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u] = tmp;}
                }
            }
        }

    // ExFv
    ptrs2(H, HExFv, Nfrag, Nvrt, 1,2,2,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    for (int w=0;  w<Nvrt[M2];  w++)
                        {
                        tmp = H2[M2*Nfrag + M1][(v+1)*Ntot[M1]*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M2]*Ntot[M1] + (w+1)*Ntot[M1] + 0];
                        if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2]*Nvrt[M2] + v*Nvrt[M2] + w] = tmp;}
                        }
                    }
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    for (int w=0;  w<Nvrt[M2];  w++)
                        {
                        tmp = H2[M1*Nfrag + M2][(u+1)*Ntot[M2]*Ntot[M1]*Ntot[M2] + (v+1)*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + (w+1)];
                        if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2]*Nvrt[M2] + v*Nvrt[M2] + w] = tmp;}
                        }
                    }
                }
            }
        }

    // ExDx
    ptrs2(H, HExDx, Nfrag, Nvrt, 1,2,0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    tmp = H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M2]*Ntot[M1] + (v+1)*Ntot[M1] + 0];
                    if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2] + v] = tmp;}
                    }
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    tmp = H2[M1*Nfrag + M2][(u+1)*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + (v+1)];
                    if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2] + v] = tmp;}
                    }
                }
            }
        }

    // FoFo
    ptrs2(H, HFoFo, Nfrag, Nvrt, 0,0,0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            tmp = (1/2.) * H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + 0*Ntot[M2]*Ntot[M1] + 0*Ntot[M1] + 0];
            if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][0] = tmp;}
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            tmp = (1/2.) * H2[M1*Nfrag + M2][0*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + 0];
            if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][0] = tmp;}
            }
        }

    // FoFv
    ptrs2(H, HFoFv, Nfrag, Nvrt, 2,2,0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M2];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    tmp = -H2[M2*Nfrag + M1][(u+1)*Ntot[M1]*Ntot[M2]*Ntot[M1] + 0*Ntot[M2]*Ntot[M1] + (v+1)*Ntot[M1] + 0];
                    if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2] + v] = tmp;}
                    }
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M2];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    tmp = -H2[M1*Nfrag + M2][0*Ntot[M2]*Ntot[M1]*Ntot[M2] + (u+1)*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + (v+1)];
                    if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2] + v] = tmp;}
                    }
                }
            }
        }

    // FoDx
    ptrs2(H, HFoDx, Nfrag, Nvrt, 2,0,0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M2];  u++)
                {
                tmp = -H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + 0*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M1] + 0];
                if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u] = tmp;}
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M2];  u++)
                {
                tmp = -H2[M1*Nfrag + M2][0*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + 0*Ntot[M2] + (u+1)];
                if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u] = tmp;}
                }
            }
        }

    // FvFv
    ptrs2(H, HFvFv, Nfrag, Nvrt, 1,2,1,2);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    for (int w=0;  w<Nvrt[M1];  w++)
                        {
                        for (int x=0;  x<Nvrt[M2];  x++)
                            {
                            tmp = (1/2.) * H2[M2*Nfrag + M1][(v+1)*Ntot[M1]*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M2]*Ntot[M1] + (x+1)*Ntot[M1] + (w+1)];
                            if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2]*Nvrt[M1]*Nvrt[M2] + v*Nvrt[M1]*Nvrt[M2] + w*Nvrt[M2] + x] = tmp;}
                            }
                        }
                    }
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    for (int w=0;  w<Nvrt[M1];  w++)
                        {
                        for (int x=0;  x<Nvrt[M2];  x++)
                            {
                            tmp = (1/2.) * H2[M1*Nfrag + M2][(u+1)*Ntot[M2]*Ntot[M1]*Ntot[M2] + (v+1)*Ntot[M1]*Ntot[M2] + (w+1)*Ntot[M2] + (x+1)];
                            if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2]*Nvrt[M1]*Nvrt[M2] + v*Nvrt[M1]*Nvrt[M2] + w*Nvrt[M2] + x] = tmp;}
                            }
                        }
                    }
                }
            }
        }

    // FvDx
    ptrs2(H, HFvDx, Nfrag, Nvrt, 1,1,2,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M1];  v++)
                    {
                    for (int w=0;  w<Nvrt[M2];  w++)
                        {
                        tmp = H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + (u+1)*Ntot[M2]*Ntot[M1] + (w+1)*Ntot[M1] + (v+1)];
                        if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M1]*Nvrt[M2] + v*Nvrt[M2] + w] = tmp;}
                        }
                    }
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M1];  v++)
                    {
                    for (int w=0;  w<Nvrt[M2];  w++)
                        {
                        tmp = H2[M1*Nfrag + M2][(u+1)*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + (v+1)*Ntot[M2] + (w+1)];
                        if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M1]*Nvrt[M2] + v*Nvrt[M2] + w] = tmp;}
                        }
                    }
                }
            }
        }

    // DxDx
    ptrs2(H, HDxDx, Nfrag, Nvrt, 1,2,0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<M1;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    tmp = (1/2.) * H2[M2*Nfrag + M1][0*Ntot[M1]*Ntot[M2]*Ntot[M1] + 0*Ntot[M2]*Ntot[M1] + (v+1)*Ntot[M1] + (u+1)];
                    if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2] + v] = tmp;}
                    }
                }
            }
        for (int M2=M1+1;  M2<Nfrag;  M2++)
            {
            for (int u=0;  u<Nvrt[M1];  u++)
                {
                for (int v=0;  v<Nvrt[M2];  v++)
                    {
                    tmp = (1/2.) * H2[M1*Nfrag + M2][0*Ntot[M2]*Ntot[M1]*Ntot[M2] + 0*Ntot[M1]*Ntot[M2] + (u+1)*Ntot[M2] + (v+1)];
                    if (fabs(tmp)>thresh) {H[M1*Nfrag + M2][u*Nvrt[M2] + v] = tmp;}
                    }
                }
            }
        }

    free(H);
    free(Ntot);

    return;
    }
