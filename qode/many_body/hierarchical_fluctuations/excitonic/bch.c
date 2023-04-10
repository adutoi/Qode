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
#include <stdlib.h>		// malloc for workspace
#include "PyC_types.h"
#include "partition_storage.h"	// Defines ptrs1 and ptrs2



// Stores the (combined) lengths of the dimensions that both precede (rows) and follow (columns) the dimension to be contracted
typedef struct
	{
	BigInt  rows;
	BigInt  cols;
	Double* data;
	}
tensor;
tensor TNS(BigInt rows, Double* data, BigInt cols)
	{
	tensor value;
	value.rows = rows;
	value.cols = cols;
	value.data = data;
	return value;
	}

// Primitive operations on vectors (or tensors treated as vectors)
void zero(Double* A, BigInt dim)
	{
	for (BigInt i=0;  i<dim;  i++) {A[i] = 0.;}
	return;
	}
void add(Double* Z, Double a, Double* B, BigInt dim)
	{
	if      (a ==  1.)
		{
		for (BigInt i=0;  i<dim;  i++) {Z[i] +=   B[i];}
		}
	else if (a == -1.)
		{
		for (BigInt i=0;  i<dim;  i++) {Z[i] -=   B[i];}
		}
	else if (a !=  0.)
		{
		for (BigInt i=0;  i<dim;  i++) {Z[i] += a*B[i];}
		}
	return;
	}

// Contract two tensors, A and B, where the length of the contracted dimension is ABdot, and the result scaled by c is stored at location Z
void contract(BigInt c, Double* Z, tensor Atns, BigInt ABdot, tensor Btns)
	{
	Double* A     = Atns.data;
	BigInt  Arows = Atns.rows;
	BigInt  Acols = Atns.cols;
	Double* B     = Btns.data;
	BigInt  Brows = Btns.rows;
	BigInt  Bcols = Btns.cols;
	for (BigInt Ar=0;  Ar<Arows;  Ar++)
		{
		for (BigInt Ac=0;  Ac<Acols;  Ac++)
			{
			Double* a = A + Ar*ABdot*Acols + Ac;
			for (BigInt Br=0;  Br<Brows;  Br++)
				{
				Double* b = B + Br*ABdot*Bcols;
				Double* z = Z + Ar*Acols*Brows*Bcols + Ac*Brows*Bcols + Br*Bcols;
				BigInt  j = 0;
				if (Bcols==1)
					{
					if      (c==1)
						{
						for (BigInt i=0;  i<ABdot;  i++)
							{
							z[0] += a[j] * b[i];
							j += Acols;
							}
						}
					else if (c==-1)
						{
						for (BigInt i=0;  i<ABdot;  i++)
							{
							z[0] -= a[j] * b[i];
							j += Acols;
							}
						}
					else
						{
						for (BigInt i=0;  i<ABdot;  i++)
							{
							z[0] += c * a[j] * b[i];
							j += Acols;
							}
						}
					}
				else
					{
					for (BigInt i=0;  i<ABdot;  i++)
						{
						add(z, c*a[j], b, Bcols);
						b += Bcols;
						j += Acols;
						}
					}
				}
			}
		}
	return;
	}



void execute_Ex(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		// + \mathcal{H}^{u_m}_{o_m}
		add(O1[M1], 1., H[M1], Nv[M1]);
		}
	return;
	}

void execute_Fo(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		// + \mathcal{H}^{o_m}_{o_m} t_1^{u_m}
		add(O1[M1], H[M1][0], t1[M1], Nv[M1]);
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{o_{m_1}}_{o_{m_1}} t_2^{u_{m_1} u_{m_2}}
			if (M1!=M2) {add(O2[M1*Nfrag+M2], H[M1][0], t2[M1*Nfrag+M2], Nv[M1]*Nv[M2]);}
			}
		}
	return;
	}

void execute_Fv(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		// + \mathcal{H}^{u_m}_{\dot{x}_m} t_1^{\dot{x}_m}
		contract(+1, O1[M1], TNS(Nv[M1], H[M1], 1), Nv[M1], TNS(1, t1[M1], 1));
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{u_{m_1}}_{\dot{x}_{m_1}} t_2^{\dot{x}_{m_1} u_{m_2}}
			if (M1!=M2) {contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1], H[M1], 1), Nv[M1], TNS(1, t2[M1*Nfrag+M2], Nv[M2]));}
			}
		}
	return;
	}

void execute_Dx(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				zero(tmp, Nv[M1]);
				contract(+1, tmp, TNS(1, H[M2], 1), Nv[M2], TNS(1, t2[M2*Nfrag+M1], Nv[M1]));
				// - \mathcal{H}^{o_{m_2}}_{\dot{x}_{m_2}} t_2^{\dot{x}_{m_2} u_{m_1}} t_1^{u_{m_2}}
				contract(-1, O2[M1*Nfrag+M2], TNS(Nv[M1], tmp, 1), 1, TNS(1, t1[M2], Nv[M2]));
				// + \mathcal{H}^{o}_{\dot{x}} t_2^{\dot{x} u_{m}}
				add(O1[M1], 1., tmp, Nv[M1]);
				}
			}
		}
	Double E = 0.;
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		tmp[0] = 0.;
		contract(+1, tmp, TNS(1, H[M1], 1), Nv[M1], TNS(1, t1[M1], 1));
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// - \mathcal{H}^{o_{m_1}}_{\dot{x}_{m_1}} t_1^{\dot{x}_{m_1}} t_2^{u_{m_1} u_{m_2}}
			if (M1!=M2) {add(O2[M1*Nfrag+M2], -tmp[0], t2[M1*Nfrag+M2], Nv[M1]*Nv[M2]);}
			}
		// - \mathcal{H}^{o_m}_{\dot{x}_m} t_1^{\dot{x}_m} t_1^{u_m}
		add(O1[M1], -tmp[0], t1[M1], Nv[M1]);
		E += tmp[0];
		}
	// + \mathcal{H}^{o}_{\dot{x}} t_1^{\dot{x}}
	O0[0] += E;
	return;
	}

void execute_ExEx(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{u_{m_1} u_{m_2}}_{o_{m_1} o_{m_2}}
			if (M1!=M2) {add(O2[M1*Nfrag+M2], 1., H[M1*Nfrag+M2], Nv[M1]*Nv[M2]);}
			}
		}
	return;
	}

void execute_ExFo(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{u_{m_1} o_{m_2}}_{o_{m_1} o_{m_2}} t_1^{u_{m_2}}
			if (M1!=M2) {contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1], H[M1*Nfrag+M2], 1), 1, TNS(1, t1[M2], Nv[M2]));}
			}
		}
	return;
	}

void execute_ExFv(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{u_{m_1} u_{m_2}}_{o_{m_1} \dot{x}_{m_2}} t_1^{\dot{x}_{m_2}}
			if (M1!=M2) {contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1]*Nv[M2], H[M1*Nfrag+M2], 1), Nv[M2], TNS(1, t1[M2], 1));}
			}
		}
	return;
	}

void execute_ExDx(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				for (BigInt N=0;  N<Nfrag;  N++)
					{
					if (N!=M1 && N!=M2)
						{
						// + \mathcal{H}^{u_{m_1} o}_{o_{m_1} \dot{x}} t_2^{\dot{x} u_{m_2}}
						contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1], H[M1*Nfrag+N], 1), Nv[N], TNS(1, t2[N*Nfrag+M2], Nv[M2]));
						}
					}
				}
			}
		}
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				zero(tmp, Nv[M1]);
				contract(+1, tmp, TNS(Nv[M1], H[M1*Nfrag+M2], 1), Nv[M2], TNS(1, t1[M2], 1));
				// - \mathcal{H}^{u_{m_1} o_{m_2}}_{o_{m_1} \dot{x}_{m_2}} t_1^{\dot{x}_{m_2}} t_1^{u_{m_2}}
				contract(-1, O2[M1*Nfrag+M2], TNS(Nv[M1], tmp, 1), 1, TNS(1, t1[M2], Nv[M2]));
				// + \mathcal{H}^{u_{m} o}_{o_{m} \dot{x}} t_1^{\dot{x}}
				add(O1[M1], 1., tmp, Nv[M1]);
				}
			}
		}
	return;
	}

void execute_FoFo(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{o_{m_1} o_{m_2}}_{o_{m_1} o_{m_2}} t_\text{II}^{u_{m_1} u_{m_2}}
			if (M1!=M2) {add(O2[M1*Nfrag+M2], H[M1*Nfrag+M2][0], tII[M1*Nfrag+M2], Nv[M1]*Nv[M2]);}
			}
		}
	return;
	}

void execute_FoFv(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{o_{m_2} u_{m_1}}_{o_{m_2} \dot{x}_{m_1}} t_\text{II}^{\dot{x}_{m_1} u_{m_2}}
			if (M1!=M2) {contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1], H[M2*Nfrag+M1], 1), Nv[M1], TNS(1, tII[M1*Nfrag+M2], Nv[M2]));}
			}
		}
	return;
	}

void execute_FoDx(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	Double* tmp1 = tmp;
	Double* tmp2 = tmp1 + Nfrag*NvMax;
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		zero(tmp1, Nfrag*Nv[M1]);
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			for (BigInt N=0;  N<Nfrag;  N++)
				{
				if (N!=M1 && N!=M2)
					{
					if (M1==M2)
						{
						zero(tmp2, Nv[M1]);
						contract(+1, tmp2, TNS(1, H[M1*Nfrag+N], 1), Nv[N], TNS(1, t2[N*Nfrag+M1], Nv[M1]));
						add(tmp1+M2*Nv[M1], 1., tmp2, Nv[M1]);
						add(tmp1+N*Nv[M1], -1., tmp2, Nv[M1]);
						}
					else
						{
						contract(+1, tmp1+M2*Nv[M1], TNS(1, H[M2*Nfrag+N], 1), Nv[N], TNS(1, t2[N*Nfrag+M1], Nv[M1]));
						}
					}
				}
			}
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + ( \mathcal{H}^{o_{m_2} o}_{o_{m_2} \dot{x}} t_2^{\dot{x} u_{m_1}} - \mathcal{H}^{o_{m_1} o_{m_2}}_{o_{m_1} \dot{x}_{m_2}} t_2^{\dot{x}_{m_2} u_{m_1}} ) t_1^{u_{m_2}}
			if (M1!=M2) {contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1], tmp1+M2*Nv[M1], 1), 1, TNS(1, t1[M2], Nv[M2]));}
			}
		}
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		Double h = 0.;
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				tmp[0] = 0.;
				contract(+1, tmp, TNS(1, H[M1*Nfrag+M2], 1), Nv[M2], TNS(1, t1[M2], 1));
				// - \mathcal{H}^{o_{m_1} o_{m_2}}_{o_{m_1} \dot{x}_{m_2}} t_1^{\dot{x}_{m_2}} t_\text{II}^{u_{m_1} u_{m_2}}
				add(O2[M1*Nfrag+M2], -tmp[0], tII[M1*Nfrag+M2], Nv[M1]*Nv[M2]);
				h += tmp[0];
				}
			}
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{o_{m_1} o}_{o_{m_1} \dot{x}} t_1^{\dot{x}} t_2^{u_{m_1} u_{m_2}}
			if (M1!=M2) {add(O2[M1*Nfrag+M2], h, t2[M1*Nfrag+M2], Nv[M1]*Nv[M2]);}
			}
		}
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt N=0;  N<Nfrag;  N++)
			{
			// + \mathcal{H}^{o_{m} o}_{o_{m} \dot{x}} t_\text{II}^{\dot{x} u_{m}}
			if (M1!=N) {contract(+1, O1[M1], TNS(1, H[M1*Nfrag+N], 1), Nv[N], TNS(1, tII[N*Nfrag+M1], Nv[M1]));}
			}
		}
	return;
	}

void execute_FvFv(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// + \mathcal{H}^{u_{m_1} u_{m_2}}_{\dot{x}_{m_1} \dot{x}_{m_2}} t_\text{II}^{\dot{x}_{m_1} \dot{x}_{m_2}}
			if (M1!=M2) {contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1]*Nv[M2], H[M1*Nfrag+M2], 1), Nv[M1]*Nv[M2], TNS(1, tII[M1*Nfrag+M2], 1));}
			}
		}
	return;
	}

void execute_FvDx(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	Double* tmp1 = tmp;
	Double* tmp2 = tmp1 + NvMax*NvMax*Nfrag;
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				for (BigInt N=0;  N<Nfrag;  N++)
					{
					if (N!=M1 && N!=M2)
						{
						zero(tmp, Nv[M1]*Nv[N]);
						contract(+1, tmp, TNS(Nv[M1], H[M1*Nfrag+N], Nv[N]), Nv[M1], TNS(1, t1[M1], 1));
						// + \mathcal{H}^{u_{m_1} o}_{\dot{x}_{m_1} \dot{x}} t_1^{\dot{x}_{m_1}} t_2^{\dot{x} u_{m_2}}
						contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1], tmp, 1), Nv[N], TNS(1, t2[N*Nfrag+M2], Nv[M2]));
						}
					}
				}
			}
		zero(tmp1, Nv[M1]*Nv[M1]*Nfrag);
		zero(tmp2, Nv[M1]*Nv[M1]);
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				contract(-1, tmp1+M2*Nv[M1]*Nv[M1], TNS(Nv[M1]*Nv[M1], H[M1*Nfrag+M2], 1), Nv[M2], TNS(1, t1[M2], 1));
				add(tmp2, -1., tmp1+M2*Nv[M1]*Nv[M1], Nv[M1]*Nv[M1]);
				}
			}
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				add(tmp1+M2*Nv[M1]*Nv[M1], 1., tmp2, Nv[M1]*Nv[M1]);
				// + ( \mathcal{H}^{u_{m_1} o}_{\dot{x}_{m_1} \dot{x}} t_1^{\dot{x}} - \mathcal{H}^{u_{m_1} o_{m_2}}_{\dot{x}_{m_1} \dot{x}_{m_2}} t_1^{\dot{x}_{m_2}} ) t_2^{\dot{x}_{m_1} u_{m_2}}
				contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1], tmp1+M2*Nv[M1]*Nv[M1], 1), Nv[M1], TNS(1, t2[M1*Nfrag+M2], Nv[M2]));
				}
			}
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				zero(tmp, Nv[M1]);
				contract(+1, tmp, TNS(Nv[M1], H[M1*Nfrag+M2], 1), Nv[M1]*Nv[M2], TNS(1, tII[M1*Nfrag+M2], 1));
				// - \mathcal{H}^{u_{m_1} o_{m_2}}_{\dot{x}_{m_1} \dot{x}_{m_2}} t_\text{II}^{\dot{x}_{m_1} \dot{x}_{m_2}} t_1^{u_{m_2}}
				contract(-1, O2[M1*Nfrag+M2], TNS(Nv[M1], tmp, 1), 1, TNS(1, t1[M2], Nv[M2]));
				// + \mathcal{H}^{u_{m} o}_{\dot{x}_{m} \dot{x}} t_\text{II}^{\dot{x}_{m} \dot{x}}
				add(O1[M1], 1., tmp, Nv[M1]);
				}
			}
		}
	return;
	}

void execute_DxDx(BigInt nthd, BigInt Nfrag, BigInt* Nv, BigInt NvMax, Double* O0, Double** O1, Double** O2, Double** t1, Double** t2, Double** tII, Double* tmp, Double** H)
	{
	Double* tmp1 = tmp;
	Double* tmp2 = tmp1 + ((Nfrag>NvMax) ? Nfrag*NvMax : NvMax*NvMax);
	Double* tmp3 = tmp2 + NvMax*NvMax;
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt N=0;  N<Nfrag;  N++)
			{
			zero(tmp1, Nv[N]*Nv[M1]);
			for (BigInt M2=0;  M2<Nfrag;  M2++)
				{
				if (M1!=M2 && N!=M2)
					{
					zero(tmp2, Nv[N]*Nv[M1]);
					contract(+1, tmp2, TNS(Nv[N], H[N*Nfrag+M2], 1), Nv[M2], TNS(1, t2[M2*Nfrag+M1], Nv[M1]));
					// - 2 \mathcal{H}^{o o_{m_2}}_{\dot{x} \dot{x}_{m_2}} t_2^{\dot{x}_{m_2} u_{m_1}} t_2^{\dot{x} u_{m_2}}
					contract(-2, O2[M1*Nfrag+M2], TNS(1, tmp2, Nv[M1]), Nv[N], TNS(1, t2[N*Nfrag+M2], Nv[M2]));
					// + \mathcal{H}^{o_{m_1} o_{m_2}}_{\dot{x}_{m_1} \dot{x}_{m_2}} t_2^{\dot{x}_{m_2} u_{m_1}} t_2^{\dot{x}_{m_1} u_{m_2}}
					if (N==M1) {contract(+1, O2[M1*Nfrag+M2], TNS(1, tmp2, Nv[M1]), Nv[M1], TNS(1, t2[M1*Nfrag+M2], Nv[M2]));}
					add(tmp1, 1., tmp2, Nv[N]*Nv[M1]);
					}
				}
			for (BigInt M2=0;  M2<Nfrag;  M2++)
				{
				// + \mathcal{H}^{o o}_{\dot{x} \dot{x}'} t_2^{\dot{x}' u_{m_1}} t_2^{\dot{x} u_{m_2}}
				if (M1!=M2 && N!=M2) {contract(+1, O2[M1*Nfrag+M2], TNS(1, tmp1, Nv[M1]), Nv[N], TNS(1, t2[N*Nfrag+M2], Nv[M2]));}
				}
			}
		}
	for (BigInt M2=0;  M2<Nfrag;  M2++)
		{
		zero(tmp1, Nv[M2]*Nfrag);
		zero(tmp2, Nv[M2]);
		for (BigInt M1=0;  M1<Nfrag;  M1++)
			{
			if (M1!=M2)
				{
				contract(+2, tmp1+M1*Nv[M2], TNS(Nv[M2], H[M2*Nfrag+M1], 1), Nv[M1], TNS(1, t1[M1], 1));
				add(tmp2, 1., tmp1+M1*Nv[M2], Nv[M2]);
				}
			}
		for (BigInt M1=0;  M1<Nfrag;  M1++)
			{
			add(tmp1+M1*Nv[M2], -1., tmp2, Nv[M2]);
			}
		for (BigInt M1=0;  M1<Nfrag;  M1++)
			{
			zero(tmp2, Nv[M1]);
			if (M1!=M2)
				{
				contract(+1, tmp2, TNS(1, tmp1+M1*Nv[M2], 1), Nv[M2], TNS(1, t2[M2*Nfrag+M1], Nv[M1]));
				// + 2 ( \mathcal{H}^{o o}_{\dot{x} \dot{x}'} t_1^{\dot{x}'} - \mathcal{H}^{o o_{m}}_{\dot{x} \dot{x}_{m}} t_1^{\dot{x}_{m}} ) t_2^{\dot{x} u_{m}}
				add(O1[M1], -1., tmp2, Nv[M1]);
				for (BigInt N=0;  N<Nfrag;  N++)
					{
					if (N!=M1 && N!=M2)
						{
						zero(tmp3, Nv[N]);
						contract(+2, tmp3, TNS(Nv[N], H[N*Nfrag+M2], 1), Nv[M2], TNS(1, t1[M2], 1));
						contract(-1, tmp2, TNS(1, tmp3, 1), Nv[N], TNS(1, t2[N*Nfrag+M1], Nv[M1]));
						}
					}
				// - 2 (\mathcal{H}^{o o_{m_2}}_{\dot{x} \dot{x}_{m_2}} t_1^{\dot{x}_{m_2}} t_2^{\dot{x} u_{m_1}}
				//      + (\mathcal{H}^{o_{m_2} o}_{\dot{x}_{m_2} \dot{x}} t_1^{\dot{x}} - \mathcal{H}^{o_{m_2} o_{m_1}}_{\dot{x}_{m_2} \dot{x}_{m_1}} t_1^{\dot{x}_{m_1}} ) t_2^{\dot{x}_{m_2} u_{m_1}} ) t_1^{u_{m_2}}
				contract(+1, O2[M1*Nfrag+M2], TNS(Nv[M1], tmp2, 1), 1, TNS(1, t1[M2], Nv[M2]));
				}
			}
		}
	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		Double E = 0.;
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				tmp[0] = 0.;
				contract(+1, tmp, TNS(1, H[M1*Nfrag+M2], 1), Nv[M1]*Nv[M2], TNS(1, tII[M1*Nfrag+M2], 1));
				// + \mathcal{H}^{o_{m_1} o_{m_2}}_{\dot{x}_{m_1} \dot{x}_{m_2}} t_\text{II}^{\dot{x}_{m_1} \dot{x}_{m_2}} t_\text{II}^{u_{m_1} u_{m_2}}
				add(O2[M1*Nfrag+M2], tmp[0], tII[M1*Nfrag+M2], Nv[M1]*Nv[M2]);
				E += tmp[0];
				}
			}
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			// - 2 \, \mathcal{H}^{o_{m_1} o}_{\dot{x}_{m_1} \dot{x}} t_\text{II}^{\dot{x}_{m_1} \dot{x}} t_2^{u_{m_1} u_{m_2}}
			if (M1!=M2) {add(O2[M1*Nfrag+M2], -2*E, t2[M1*Nfrag+M2], Nv[M1]*Nv[M2]);}
			}
		// - 2 \, \mathcal{H}^{o_{m} o}_{\dot{x}_{m} \dot{x}} t_\text{II}^{\dot{x}_m \dot{x}} t_1^{u_{m}}
		add(O1[M1], -2*E, t1[M1], Nv[M1]);
		// + \mathcal{H}^{o o}_{\dot{x} \dot{x}'} t_\text{II}^{\dot{x} \dot{x}'}
		O0[0] += E;
		}
	return;
	}



void execute(
	PyInt nthd, PyInt Nfrag, BigInt* Nv,
	Double* Hc,
	Double* HEx, Double* HFo, Double* HFv, Double* HDx,
	Double* HExEx, Double* HExFo, Double* HExFv, Double* HExDx, Double* HFoFo, Double* HFoFv, Double* HFoDx, Double* HFvFv, Double* HFvDx, Double* HDxDx,
	Double* TEx, Double* TExEx,
	Double* OK, Double* OEx, Double* OExEx
	)
	{
	BigInt NvMax  = 0;
	BigInt NvTot  = 0;
	BigInt SumNv2 = 0;
        for (BigInt M=0;  M<Nfrag;  M++)
		{
		NvTot  += Nv[M];
		SumNv2 += Nv[M]*Nv[M];
		if (Nv[M]>NvMax) {NvMax = Nv[M];}
		}

	Double*  tmp       = (Double*)malloc(sizeof(Double)*(((Nfrag+1)*NvMax+1)*NvMax));
	Double*  tII_alloc = (Double*)malloc(sizeof(Double)*(NvTot*NvTot-SumNv2));	// do not need physical space for diagonal blocks
	Double** Hptrs     = (Double**)malloc(sizeof(Double*)*Nfrag*Nfrag);

	Double*  O0  = OK;
	Double** O1  = ptrs1((Double**)malloc(sizeof(Double*)*Nfrag),       OEx,       Nfrag, Nv, 1,0);
	Double** O2  = ptrs2((Double**)malloc(sizeof(Double*)*Nfrag*Nfrag), OExEx,     Nfrag, Nv, 1,2,0,0);
	Double** t1  = ptrs1((Double**)malloc(sizeof(Double*)*Nfrag),       TEx,       Nfrag, Nv, 1,0);
	Double** t2  = ptrs2((Double**)malloc(sizeof(Double*)*Nfrag*Nfrag), TExEx,     Nfrag, Nv, 1,2,0,0);
	Double** tII = ptrs2((Double**)malloc(sizeof(Double*)*Nfrag*Nfrag), tII_alloc, Nfrag, Nv, 1,2,0,0);

	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=0;  M2<Nfrag;  M2++)
			{
			if (M1!=M2)
				{
				zero(tII[M1*Nfrag+M2], Nv[M1]*Nv[M2]);
				add(tII[M1*Nfrag+M2], 1., t2[M1*Nfrag+M2], Nv[M1]*Nv[M2]);
				contract(+1, tII[M1*Nfrag+M2], TNS(Nv[M1], t1[M1], 1), 1, TNS(1, t1[M2], Nv[M2]));
				}
			}
		}

	// O0[0] += Hc[0];	// Take care of the constant part of the energy externally
	execute_Ex(  nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs1(Hptrs, HEx,   Nfrag, Nv, 1,0));
	execute_Fo(  nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs1(Hptrs, HFo,   Nfrag, Nv, 0,0));
	execute_Fv(  nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs1(Hptrs, HFv,   Nfrag, Nv, 1,1));
	execute_Dx(  nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs1(Hptrs, HDx,   Nfrag, Nv, 1,0));
	execute_ExEx(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HExEx, Nfrag, Nv, 1,2,0,0));
	execute_ExFo(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HExFo, Nfrag, Nv, 1,0,0,0));
	execute_ExFv(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HExFv, Nfrag, Nv, 1,2,2,0));
	execute_ExDx(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HExDx, Nfrag, Nv, 1,2,0,0));
	execute_FoFo(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HFoFo, Nfrag, Nv, 0,0,0,0));
	execute_FoFv(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HFoFv, Nfrag, Nv, 2,2,0,0));
	execute_FoDx(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HFoDx, Nfrag, Nv, 2,0,0,0));
	execute_FvFv(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HFvFv, Nfrag, Nv, 1,2,1,2));
	execute_FvDx(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HFvDx, Nfrag, Nv, 1,1,2,0));
	execute_DxDx(nthd, Nfrag, Nv, NvMax, O0, O1, O2, t1, t2, tII, tmp, ptrs2(Hptrs, HDxDx, Nfrag, Nv, 1,2,0,0));

	for (BigInt M1=0;  M1<Nfrag;  M1++)
		{
		for (BigInt M2=M1+1;  M2<Nfrag;  M2++)
			{
			for (BigInt u1=0;  u1<Nv[M1];  u1++)
				{
				for (BigInt u2=0;  u2<Nv[M2];  u2++)
					{
					tmp[0] = O2[M1*Nfrag+M2][u1*Nv[M2]+u2] + O2[M2*Nfrag+M1][u2*Nv[M1]+u1];
					O2[M1*Nfrag+M2][u1*Nv[M2]+u2] = tmp[0];
					O2[M2*Nfrag+M1][u2*Nv[M1]+u1] = tmp[0];
					}
				}
			}
		}

	free(tII);
	free(t2);
	free(t1);
	free(O2);
	free(O1);
	free(Hptrs);
	free(tII_alloc);
	free(tmp);

	return;
	}
