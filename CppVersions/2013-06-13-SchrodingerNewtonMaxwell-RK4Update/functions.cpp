// functions.cpp
#include "SN_3D.h"



void RK4 ( Schrodinger elc, Domain dom, Derivative D2, Manager boss, CX dt )
{
		
	CX* k1 = new CX [dom.N3];
	CX* k2 = new CX [dom.N3];
	CX* k3 = new CX [dom.N3];
	CX* k4 = new CX [dom.N3];
	CX* Utemp = new CX [dom.N3];
	
	CX coef1 = dt*1.0/6.0;
	CX coef2 = dt*2.0/6.0;
	CX coef3 = dt*2.0/6.0;
	CX coef4 = dt*1.0/6.0;

	int Nt;

	if ( real(dt) == 0)
		Nt = dom.Nt_im;
	else
		Nt = dom.Nt_re;

	for (int r = 0; r < Nt; r++)
	{
		
		feval ( 1, Utemp, elc, k1, k1, dom, D2 , dt);
		feval ( 2, Utemp, elc, k1, k2, dom, D2 , dt);
		feval ( 3, Utemp, elc, k2, k3, dom, D2 , dt);
		feval ( 4, Utemp, elc, k3, k4, dom, D2 , dt);

		for (int i = 0; i < dom.N3; i++)
			elc.U[i] = elc.U[i] + (coef1*k1[i] + coef2*k2[i] + coef3*k3[i] + coef4*k4[i] );

		elc.density();

		// Normalize
		if ( real(dt) == 0.0)
		{
			double sum = trapz3D(elc.rho, dom.Nx, dom.Ny, dom.Nz, dom.dx, dom.dy, dom.dz);
			for (int i = 0; i < dom.N3; i++)
				elc.U[i] = elc.U[i]/sqrt(sum);
		}

		boss.vidWrite(r, elc.rho, elc.U, dt);

	}


	delete [] k1;
	delete [] k2;
	delete [] k3;
	delete [] k4;
	delete [] Utemp;

}

void feval ( int stage, CX* Utemp,  Schrodinger elc, CX* km, CX* kn, Domain dom, Derivative D2, CX dt)
{
	CX coef;
	if (stage == 1)	{
		coef = 0.0;
	} else if (stage == 2)	{
		coef = dt/2.0;
	} else if (stage == 3)	{
		coef = dt/2.0;
	} else if (stage == 4)	{
		coef = dt;
	}

	for (int i = 0; i < dom.N3; i++)
		Utemp[i] = elc.U[i] + coef*km[i];

	RHS(dom, Utemp, kn, D2, elc);

}

void RHS (Domain dom, CX* Utemp,  CX * kn, Derivative D2, Schrodinger elc  )
{
		D2.derivCX(Utemp, kn);

		for (int i = 0; i < dom.N3; i++)
			kn[i] = CX(0.0, 1.0)* kn[i] - CX(0.0, 1.0)*elc.Vext[i]*Utemp[i];



		//for (int i = 0; i < dom.N3; i++)
			//kn[i] = elc.c1*kn[i] + elc.c2*elc.Vext[i]*Utemp[i];

		//for (int i = 0; i < dom.N3; i++)
			//kn[i] = elc.c1*kn[i] + elc.c2*elc.Vext[i]*Utemp[i];
		//for (int i = 0; i < dom.N3; i++)
			//kn[i] = CX(0,1)*Utemp[i];
}



