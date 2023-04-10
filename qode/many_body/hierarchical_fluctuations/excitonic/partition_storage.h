/*   (C) Copyright 2018, 2019 Anthony D. Dutoi
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
#include "PyC_types.h"



// Assign fragment blocked pointers based on how many indices there are per fragment(s)

Double** ptrs1(Double** Aptrs, Double* A, BigInt Nfrag, BigInt* Nv, BigInt i1, BigInt i2)
	{
	BigInt inc;
	BigInt ptr = 0;
	BigInt idx = 0;
	for (BigInt M=0;  M<Nfrag;  M++)
		{
		Aptrs[ptr++] = A + idx;
		inc = 1;
		inc *= (i1==0) ? 1 : Nv[M];
		inc *= (i2==0) ? 1 : Nv[M];
		idx += inc;
		}
	return Aptrs;
	}

Double** ptrs2(Double** Aptrs, Double* A, BigInt Nfrag, BigInt* Nv, BigInt i1, BigInt i2, BigInt i3, BigInt i4)
	{
	BigInt inc;
	BigInt ptr = 0;
	BigInt idx = 0;
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			Aptrs[ptr++] = A + idx;
			inc = 1;
			inc *= (i1==0) ? 1 : ((i1==1) ? Nv[M1] : Nv[M2]);
			inc *= (i2==0) ? 1 : ((i2==1) ? Nv[M1] : Nv[M2]);
			inc *= (i3==0) ? 1 : ((i3==1) ? Nv[M1] : Nv[M2]);
			inc *= (i4==0) ? 1 : ((i4==1) ? Nv[M1] : Nv[M2]);
			if (M1!=M2) {idx += inc;}	// reserve zero space for M1=M2 since those blocks never referenced
			}
		}
	return Aptrs;
	}
