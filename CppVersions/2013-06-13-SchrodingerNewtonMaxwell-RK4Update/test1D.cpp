// test1D.cpp
#include "SN_3D.h"

void RK4_1D ( CX* U, int Nx, int Nt, CX dt, double dx, Derivative1D  D2, double* Vext)
{
		
	CX* k1 = new CX [Nx];
	CX* k2 = new CX [Nx];
	CX* k3 = new CX [Nx];
	CX* k4 = new CX [Nx];
	CX* Utemp = new CX [Nx];
	
	CX coef1 = dt*1.0/6.0;
	CX coef2 = dt*2.0/6.0;
	CX coef3 = dt*2.0/6.0;
	CX coef4 = dt*1.0/6.0;

	double* rho = new double [Nx];

	for (int r = 0; r < Nt; r++)
	{
		
		feval1D ( 1, U, Utemp, k1, k1, Nx, dt, D2, Vext );
		feval1D ( 2, U, Utemp, k1, k2, Nx, dt, D2, Vext );
		feval1D ( 3, U, Utemp, k2, k3, Nx, dt, D2, Vext );
		feval1D ( 4, U, Utemp, k3, k4, Nx, dt, D2, Vext );


		for (int i = 0; i < Nx; i++)
			U[i] = U[i] + (coef1*k1[i] + coef2*k2[i] + coef3*k3[i] + coef4*k4[i] );

		for (int k = 0; k < Nx; k++)
		{
			rho[k] = pow ( abs(U[k]), 2.0);
		}

		// Normalize
		if ( real(dt) == 0.0)
		{
			double sum = trapz( rho, Nx, dx);
			for (int i = 0; i < Nx; i++)
				U[i] = U[i]/sqrt(sum);
		}
		
		
		if ( mod(r,10) == 0 ){
		cout << r << "\t";
		ofstream data1("1Dvid.bin", ios::out | ios::binary | ios::app); 
		double temp;
		for (int k = 0; k < Nx; k++)
		{
			temp = real(U[k]);
			data1.write((char *) & temp, sizeof temp);
		}
		data1.close();

		}


	}


	delete [] k1;
	delete [] k2;
	delete [] k3;
	delete [] k4;
	delete [] Utemp;
	delete [] rho;
}

void feval1D ( int stage, CX* U, CX* Utemp, CX* km, CX* kn, int Nx, CX dt, Derivative1D  D2, double* Vext )
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

	for (int i = 0; i < Nx; i++)
		Utemp[i] = U[i] + coef*km[i];

	RHS_1D( Nx, Utemp, kn, D2, Vext);
}

void RHS_1D (int Nx, CX* Utemp,  CX * kn, Derivative1D D2, double* Vext )
{
		D2.derivCX(Utemp, kn);

		for (int i = 0; i < Nx; i++)
			kn[i] = CX(0.0, 1.0)* kn[i] - CX(0.0, 1.0)* Vext[i]*Utemp[i];

}

void run1D()
{
	Domain dom ;
	Manager boss(dom);
	boss.startTime();


	CX dt;
	int N, Nt;
	double tf, L;

	string boundary = "zero";
	N	= 50;
	Nt	= 800;
	tf	= 20;
	L	= 4*PI - 4*PI/(N);
	

	Derivative1D D2( 2.0, 2.0, boundary, N, L);


	CX* U = new CX [N];
	CX* dU = new CX[N];
	double* Vext = new double [N];

	for (int k = 0; k < N; k++)
	{
		double x = k*D2.dx;
		U[k] = exp(-pow(x- 2*PI,2.0)/5.0);
		Vext[k] =  -exp(-pow(x- 2*PI,4.0)/50.0);
		//U[k] = sin(x);

	}

	ofstream data1("1Dvid.bin", ios::out | ios::binary); 
	data1.close();

	dt = tf/(double) Nt*CX(0.0, -1.0);
	RK4_1D (U, N, Nt,dt, D2.dx, D2, Vext);

	dt = tf/(double) Nt*CX(1.0, 0.0);
	RK4_1D (U, N, Nt,dt, D2.dx, D2, Vext);

	
		
	ofstream data("1Ddata.bin", ios::out | ios::binary); 
	double temp;

	temp = N;
	data.write((char *) & temp, sizeof temp);

	temp = L;
	data.write((char *) & temp, sizeof temp);

	temp = D2.dx;
	data.write((char *) & temp, sizeof temp);

	for (int k = 0; k < N; k++)
	{
		temp = real(U[k]);
		data.write((char *) & temp, sizeof temp);
	}

	for (int k = 0; k < N; k++)
	{
		temp = real(Vext[k]);
		data.write((char *) & temp, sizeof temp);
	
	}

	
	data.close();


	boss.endTime();
}
