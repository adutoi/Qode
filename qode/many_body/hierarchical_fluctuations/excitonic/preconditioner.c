/*   (C) Copyright 2018, 2019 Anthony D. Dutoi and Yuhong Liu
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
#include <stdlib.h>		// malloc for array of pointers
#include "PyC_types.h"
#include "partition_storage.h"	// Defines ptrs1 and ptrs2

void divide_by_diagE(PyInt Nfrag, BigInt* Nvrt, Double* OEx, Double* OExEx, Double** diagE)
    {
    Double** O = (Double**)malloc(sizeof(Double*)*Nfrag*Nfrag);

    // Ex
    ptrs1(O, OEx, Nfrag, Nvrt, 1,0);
    for (int M=0;  M<Nfrag;  M++)
        {
        for (int u=0;  u<Nvrt[M];  u++)
            {
            O[M][u] /= (diagE[M][u+1] - diagE[M][0]);
            }
        }

    // ExEx
    ptrs2(O, OExEx, Nfrag, Nvrt, 1,2,0,0);
    for (int M1=0;  M1<Nfrag;  M1++)
        {
        for (int M2=0;  M2<Nfrag;  M2++)
            {
            if (M1!=M2)
                {
                for (int u=0;  u<Nvrt[M1];  u++)
                    {
                    for (int v=0;  v<Nvrt[M2];  v++)
                        {
                        O[M1*Nfrag + M2][u*Nvrt[M2] + v] /= (diagE[M1][u+1] + diagE[M2][v+1] - diagE[M1][0] - diagE[M2][0]);
                        }
                    }
                }
            }
        }

    free(O);

    return;
    }
