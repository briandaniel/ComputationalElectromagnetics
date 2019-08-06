// classFunctions.cpp
#include "SN_3D.h"


Domain::Domain()
{
	
	Config settings("settings.txt");
	
	// Import number of grid points
	settings.get("Domain/Nx", Nx);
	settings.get("Domain/Ny", Ny);
	settings.get("Domain/Nz", Nz);
	N3 = Nx*Ny*Nz;

	// Important length of each dimension
	settings.get("Domain/Lx", Lx); 
	settings.get("Domain/Ly", Ly);
	settings.get("Domain/Lz", Lz);
	dx = Lx/(Nx-1);
	dy = Ly/(Ny-1);
	dz = Lz/(Nz-1);

	// Import boundary condition
	settings.get("Domain/boundary", boundary);

	// Import time domain
	settings.get("Domain/tf_im", tf_im);
	settings.get("Domain/Nt_im", Nt_im);

	settings.get("Domain/tf_re", tf_re);
	settings.get("Domain/Nt_re", Nt_re);
	
	dt_im = tf_im/ (double) Nt_im;
	dt_re = tf_re/ (double) Nt_re;
	
}

Schrodinger::Schrodinger(Domain dom)
{
	Config settings("settings.txt");
	settings.get("Equations/alpha", alpha);
	settings.get("Equations/beta", beta);

	// cn	=   (			Real		,		Imag			)
	c1		= CX(			0.0			,		HBAR/(2*ME)		);
	c2		= CX(			0.0			,	-pow(E,2.0)/(HBAR)	);
	c3		= CX(			0.0			,	 pow(E,2.0)/(HBAR)	);

	U = new CX [dom.N3];
	rho = new double [dom.N3];

	for (int k = 0; k < dom.N3; k++)
		rho[k] = 0.0;

	Nx = dom.Nx;
	Ny = dom.Ny;
	Nz = dom.Nz;
	N3 = dom.N3;

	Lx = dom.Lx;
	Ly = dom.Ly;
	Lz = dom.Lz;

	xo = Lx/2.0;
	yo = Ly/2.0;
	zo = Lz/2.0;

	dz = dom.dz;
	dy = dom.dy;
	dx = dom.dx;


	// Inital U values
	for (int i = 0; i < dom.Nx; i++)
	{
		for (int j = 0; j < dom.Ny; j++)
		{
			for (int k = 0; k < dom.Nz; k++)
			{
				double z = dom.dz*k;
				double y = dom.dy*j;
				double x = dom.dx*i;

				U[ k + j*dom.Nz + i*dom.Ny*dom.Nz ] = exp ( ( -pow(x - xo, 2.0) -pow(y - yo, 2.0) -pow(z - zo, 2.0) )/5.0 );
				//U[ k + j*dom.Nz + i*dom.Ny*dom.Nz ] = sin(x) + sin(y) + sin(z) ;
				//U[ k + j*dom.Nz + i*dom.Ny*dom.Nz ] = 1.0;

			}
		}
	}

	density();

	double sum = trapz3D(rho, dom.Nx, dom.Ny, dom.Nz, dom.dx, dom.dy, dom.dz);
	cout << "Sum = " << sum;

	for (int i = 0; i < dom.N3; i++)
		U[i] = U[i]/sqrt(sum);
	
	density();
	sum = trapz3D(rho, dom.Nx, dom.Ny, dom.Nz, dom.dx, dom.dy, dom.dz);
	cout << "Sum = " << sum;

	Vext = new double [dom.N3];

	// Initial Vext values
	for (int i = 0; i < dom.Nx; i++)
	{
		for (int j = 0; j < dom.Ny; j++)
		{
			for (int k = 0; k < dom.Nz; k++)
			{
				double z = dom.dz*k;
				double y = dom.dy*j;
				double x = dom.dx*i;

				double r = sqrt( pow (z - zo,2.0) + pow (y - yo,2.0) + pow(x - xo,2.0) );
				Vext[ k + j*dom.Nz + i*dom.Ny*dom.Nz ] = -exp( -  pow(r ,4.0)/50);
			    //Vext[ k + j*dom.Nz + i*dom.Ny*dom.Nz ] = - exp ( ( -pow(x - xo, 2.0) -pow(y -yo, 2.0) -pow(z - zo, 2.0) )) ;

				//Vext[ k + j*dom.Nz + i*dom.Ny*dom.Nz ] = 10*( pow(x - xo, 2.0) + pow(y - yo, 2.0) + pow(z - zo, 2.0) )- 100.0;
			}
		}
	}


}

void Schrodinger::density (){
	for (int i = 0; i < N3; i++)
		rho[i] = pow( abs(U[i]), 2.0);
}