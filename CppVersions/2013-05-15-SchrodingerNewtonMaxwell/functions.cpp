#include "SNM.h"

// Constructs initial sech shape
void sech ( Domain dom, Schrodinger elc)
{
	for (int s = 0; s< dom.Npsi; s++)
	{
		for( int k = 0; k < dom.N; k++)
		{
			double zLoc = ( dom.zDomPsi[k] - dom.zMax/2 )/( dom.dz*10 );
			double zShift = 10*(double) s/ (double) dom.Npsi - 10/2.0;
			elc.PSI[s][k] = CX( 2*exp(-zLoc - zShift) / (1 + exp(-2*zLoc - zShift)), 0.0 );
		}
	}

}

/* Solves the NLS using a Crank-Nicolson method with a strang split to
%  compute the nonlinear part.
%  i.e. time dependent schrodinger equation with linear potential and 
%  nonlinear convolution term...
%-------------------------------------------------------------------------%
%                                                                         %  
%  dPsi      d^2Psi        dPsi       dA                                  % 
%  ---- = b1*------ + b2*A*---- + b3*----Psi + b4*A^2*Psi                 %
%   dt        dz^2          dz        dz                                  % 
%                                                                         % 
%                               + b5*V(z)*Psi + b6*F(rho)*Psi             %
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%      where rho(x) = sum |Psi_i(z)|^2                                    %
%                                                                         %
%      F(rho) = Int [ rho(z)*u(z-z') ] dz'                                %
%                i.e. F is the convolution of rho(z) and u(z)             %
%                                                                         %
%-------------------------------------------------------------------------%*/
void crank ( Domain dom, Schrodinger elc, CX dt, int S)
{ 
	
	CX g1 = elc.b1*dt / ( 2*pow(dom.dz,2.0) );
	CX g2 = elc.b2*dt / (4*dom.dz) ;
	CX g3 = elc.b3*dt / 2.0 ;
	CX g4 = elc.b4*dt;
	//CX g5 = elc.b5*dt / 2.0 ;


	g4 =( 0.0, 0.0 );

	CX* a = new CX [dom.N-1];
	CX* b = new CX [dom.N];
	CX* c = new CX [dom.N-1];
	CX* d = new CX [dom.N-1];
	CX* e = new CX [dom.N];
	CX* f = new CX [dom.N-1];


	for (int k = 0; k < dom.N-1; k++)
	{
		a[k] = -g1 - g2*elc.A[k];
		c[k] = -g1 + g2*elc.A[k+1];

		d[k] = g1 + g2*elc.A[k];
		f[k] = g1 - g2*elc.A[k+1];
	}

	for (int k = 0; k < dom.N; k++)
	{
		b[k] = 1.0 + g1*2.0 - g3*elc.Ax[k] - g4*pow(elc.A[k],2.0) - elc.pmlCoef[k];
		//b[k] = 1.0 + g1*2.0 - g3*elc.Ax[k] - g4*pow(elc.A[k],2.0);
		e[k] = 1.0 - g1*2.0 + g3*elc.Ax[k] + g4*pow(elc.A[k],2.0)  + elc.pmlCoef[k];
		//e[k] = 1.0 - g1*2.0 + g3*elc.Ax[k] + g4*pow(elc.A[k],2.0);
	}

	CXtridisolve (elc.PSI[S], a, b, c, d, e, f, dom.N);
	delete [] a;	delete [] b;	delete [] c; 	delete [] d;	delete [] e;	delete [] f;

}

// Determines electron wavefunction (not ion) ground states.
void ground (  Domain dom, Schrodinger elc, Newton ions, double imTol )
{
	ofstream datadel("new.bin", ios::out | ios::binary);
	datadel.close();

	ofstream datadel2("newIon.bin", ios::out | ios::binary);
	datadel2.close();


	ofstream data("new.bin", ios::out | ios::binary | ios::app);
	ofstream data2("newIon.bin", ios::out | ios::binary | ios::app);

	CX* oldPsi = new CX [dom.N];

	// Imaginary time step
	CX dt = CX(0.0, -dom.i_dt);

	for (int k = 0; k<dom.N; k++)
	{
		elc.A[k]  = 0.0;
		elc.Ax[k] = 0.0;
		elc.Vx[k] = 0.0;
	}

	int r = 0;
	double deltaPSI = 1e10;

	//while ( r < dom.i_tSteps && deltaPSI > imTol)

	while ( r < dom.i_tSteps && deltaPSI > imTol)
	{
	
		for (int k = 0; k<dom.N; k++)
		{
			oldPsi[k] = elc.PSI[elc.Npsi-1][k];
		}

		
		
		// Calculate linear part of solution
		for (int s = 0; s< dom.Npsi; s++)
			crank (dom, elc, dt/2.0, s);

		VxCalc( dom, elc, ions );
		// Calculate nonlinear part of solution
			strang ( dom, elc, dt);

			// Calculate linear part of solution
		for (int s = 0; s< dom.Npsi; s++)
			crank (dom, elc, dt/2.0, s);


		// Normalize during ground state evolution
		gramSchmidt( dom, elc );
		sqNormalize(elc.PSI, dom.N, dom.Npsi, dom.dz);
	

		deltaPSI = 0;
		for (int k = 0; k < (dom.N-1); k++)
		{
			deltaPSI = deltaPSI +  ( abs ( abs( elc.PSI[elc.Npsi-1][k] ) - abs( oldPsi[k]) )  + abs( abs(elc.PSI[elc.Npsi-1][k+1]) - abs( oldPsi[k+1])) )/2;
		}

		
		if (dom.gVid == 1 && mod(r, 10) == 0 )
		{

			for (int s=0; s< dom.Npsi; s++)
			{
				for (int k = 0; k < dom.N; k++)
				{
					double temp = real(elc.PSI[s][k]);
					data.write((char *) & temp, sizeof temp);
				}
			}
			

			for (int s=0; s < dom.Nion; s++)
			{
				data2.write((char *) & ions.Rion[s], sizeof ions.Rion[s]);
			}

		}


		// Notify user of current time step
		if ( mod(r, dom.i_tSteps/10) == 0 )
		{
			//if ( dom.groundIons != 1 )
				cout << "Imag. steps: " << r + 1 << " - " << r + dom.i_tSteps/10 << endl;	
		}
		r = r+1;
	}
	delete [] oldPsi;
	data.close();
	data2.close();
		
}

// Generates the density vector, i.e. rho(z) = sum_i |Psi_i|^2
void densitydz (CX** PSI, double* rho, int Npsi, int N, double dz)
{
	for (int k = 0; k< 2*N-1; k++)
		rho[k] = 0.0;

	// Parallel processing commands
	#ifdef USE_OPENMP	
		omp_set_num_threads(OMP_NUM_THREADS);
		//cout << "Max Threads: " << omp_get_max_threads()<< endl;
		#pragma omp parallel
		//cout << "Utilized threads: " << omp_get_num_threads()<< endl;
		#pragma omp for
	#endif

	for (int s = 0; s<Npsi; s++)
	{
		for (int k = 0; k<N; k++)
		{
			rho[k] = rho[k] + dz*2.0*pow( abs(PSI[s][k]),2.0);
		}
	}
}

// Generates the density vector, i.e. rho(z) = sum_i |Psi_i|^2
void density (CX** PSI, double* rho, int Npsi, int N)
{
	for (int k = 0; k< 2*N-1; k++)
		rho[k] = 0.0;

	// Parallel processing commands
	#ifdef USE_OPENMP	
		omp_set_num_threads(OMP_NUM_THREADS);
		//cout << "Max Threads: " << omp_get_max_threads()<< endl;
		#pragma omp parallel
		//cout << "Utilized threads: " << omp_get_num_threads()<< endl;
		#pragma omp for
	#endif

	for (int s = 0; s<Npsi; s++)
	{
		for (int k = 0; k<N; k++)
		{
			rho[k] = rho[k] + 2.0*pow( abs(PSI[s][k]),2.0);
		}
	}
}

// Computes the nonlinear part of the strang split exactly, this replaces "CX** PSI"
void strang ( Domain dom, Schrodinger elc, CX dt )		
{
	double* Vcoul = new double[dom.N];
	
	// Calculate total density of electrons
	densitydz (elc.PSI, elc.rho, dom.Npsi, dom.N, dom.dz);		

	// Calculate convolution with potential function
	fftConv (Vcoul, elc.rho, elc.coul, dom.N, elc.No);	

	// Calculate Exchange-Correlation Potential
	elc.VxcCalc();													

	double cof = - elc.eeCoef*pow(E,2.0)*elc.epsFactor/HBAR;

	// Coefficient for e-e potential
	CX coef = ( dt)*CX( 0.0, cof );

	// Coefficient for Vxc, exchange correlction potential
	CX coef2 = dt * CX( 0.0, elc.VxcConst/HBAR);

	CX g5 = elc.b5*dt ;
	// Apply strang split step
	for (int s = 0; s < dom.Npsi ; s++)
	{
		for (int k = 0; k< dom.N; k++)
		{
			elc.PSI[s][k] = exp( coef*Vcoul[k] + coef2*elc.Vxc[k] + g5*elc.Vx[k]) * elc.PSI[s][k] ;
		}
	}

	delete [] Vcoul;

}

// Computes ion position from velocities
void ionPosition (double* Rion, double* Vion, double dt, int Nion)
{
	for (int s = 0; s < Nion; s++)
		Rion[s] = Rion[s] + Vion[s]*dt;

}

// Calculates a single step of the nonlinear schrodinger equation in real time
void psiStep (  Domain dom, Schrodinger elc, Newton ions, Maxwell field )
{

	CX dt = CX(dom.r_dt, 0.0);

	// Calculate linear part of solution
	for (int s = 0; s< dom.Npsi; s++)
		crank (dom, elc, dt/2.0, s);

	// Calculate nonlinear part of solution
	VxCalc( dom, elc, ions );
		strang ( dom, elc, dt);

		// Calculate linear part of solution
	for (int s = 0; s< dom.Npsi; s++)
		crank (dom, elc, dt/2.0, s);
}

// Computes ion velocity from forces
void ionVelocity ( Domain dom, Schrodinger elc, Newton ions )
{
	double* Vcoul	= new double [dom.N];
	double* dVcoul	= new double [dom.N];
	double* Fpsi	= new double [dom.Nion];
	double* Fion	= new double [dom.Nion];

	// Calculate force generated by ions
	for (int s = 0; s< dom.Nion; s++)
	{
		Fion[s] = 0.0;
		for (int k = 0; k< dom.Nion; k++)
		{
			if (s != k)
			{
				double zPos = ions.Rion[s] - ions.Rion[k];
				double val  = dom.r_dt*ions.Zion[s]*ions.Zion[k]/ions.Mion[s];
				Fion[s] = Fion[s] + val*(zPos/ pow( (pow(elc.gamma,2.0) + pow(zPos,2.0) ), 1.5 ) );
			}
		}
	}

	// Calculate force generated by electrons
	densitydz (elc.PSI, elc.rho, dom.Npsi, dom.N, dom.dz);			// Calculate density of electrons
 	fftConv (Vcoul, elc.rho, elc.coul, dom.N, elc.No);				// Calculate convolution
	derivative( dVcoul, Vcoul, dom.dz, dom.N );						// Determine derivative of potential to generate force
	spline(dom.zDomPsi, dVcoul, ions.Rion, Fpsi, dom.N, dom.Nion);	// Interpolate force onto ion positions which are off-grid

	for (int k = 0; k < dom.Nion; k++)
	{
		Fpsi[k] = dom.r_dt*Fpsi[k]*ions.Zion[k]/ions.Mion[k];
	}


	for (int k = 0; k < dom.Nion; k++)
	{
		ions.Vion[k] = ions.Vion[k] + ions.c1*Fion[k] + ions.c2*Fpsi[k];
		//Vion[k] = 0.0;
	}

	delete [] dVcoul; delete [] Vcoul; delete [] Fpsi; delete [] Fion;

}

// Populates ion positions with "x" and calculates the energy and stores it in "y"
double evaluateGround ( Domain dom, Schrodinger elc, Newton ions, Energy eVals, double* x)
{
	int N = (int) floor( (dom.Nion + 1 )/2.0 ) ;
	
	for (int k = 0; k < ions.Nion; k++ )
		ions.Rion[k] = 0.0;

	ions.Rion[0] = dom.surfPos;


	double run = 0.0;

	if ( mod (dom.Nion,2.0) == 0 )
	{
		for (int i = 1; i < N; i++)
		{
			run = x[i-1] + run;
			ions.Rion[i] = dom.surfPos + run;
		}

		int j = N-1;
		for (int i = N; i < dom.Nion; i++)
		{
			run = x[j] + run;
			ions.Rion[i] = dom.surfPos + run;

			j--;
		}
	}
	else
	{
		for (int i = 1; i < N; i++)
		{
			run = x[i-1] + run;
			ions.Rion[i] = dom.surfPos + run;
		}

		int j = N-2;
		for (int i = N; i < dom.Nion; i++)
		{
			run = x[j] + run;
			ions.Rion[i] = dom.surfPos + run;
			j--;
		}
	}

	ground(dom, elc, ions, elc.imTol);
	eVals.matterEnergy( dom, elc, ions);
	return eVals.totalEnergy();

}

// Simplex search algorithm to find ground state of ions
void simplex ( double* xIn, double variation, double xTol, double yTol, Domain dom, Schrodinger elc, Newton ions, Energy eVals, int maxEvals )
{
	
	cout << "Starting positions: " << endl;
	for (int k = 0; k < dom.Nion; k++)
		cout << ions.Rion[k] << "  ";
	cout << endl;

	ofstream videoData("groundVideo.bin", ios::out | ios::binary);

	double alpha = 1.0 ;
	double rho	 = -0.5;
	double gamma = 2.0 ;
	double sigma = 0.5 ;

	int N = (int) floor( (dom.Nion + 1 )/2.0 ) ;
	int M = N+1;
	int iter = 0;

	// Initialize temporary arrays
	double ** xo = new double * [M];
	double ** x = new double * [M];

	for (int k = 0; k < M; k++ )
	{
		xo[k] = new double [N];
		x[k]  = new double [N];
	}

	double * y = new double [M];
	double * yo = new double [M];

	double * xBar = new double [N];
	double * xRef = new double [N];
	double * xe = new double [N];
	double * xc = new double [N];

	double fxr = 1e5;
	double fxc = 1e5;
	double fxe = 1e5;
	

	// Populate first set of x-values
	for ( int i = 0; i < N; i++ )
	{
		xo[0][i] = xIn[i]*0.5;
		x[0][i]  = xIn[i];
	}

	// Generate additional data points randomly with a maximum distance of "variation" along each dimension
	for ( int i = 1; i < M; i++ )
	{
		for (int j = 0; j < N; j++ )
		{
			x[i][j] = xIn[j] + variation*timeRand();
			xo[i][j] = x[i][j]*0.5;
		}
	}

	// Initial function evaluations
	for ( int i = 0; i < M; i++ )
	{
		y[i] = evaluateGround( dom, elc, ions, eVals, x[i]);
		iter = iter + 1; cout<< "Function evals: " << iter << endl;
		yo[i] = 0.5*y[i];
	}


	// Main simplex algorithm
	while ( ( abs(y[0] - y[1]) > yTol ) && (iter < maxEvals) )
	{

		// Store previous values
		for ( int i = 0; i < M; i++ )
		{
			for (int j = 0; j < N; j++ )
			{
				xo[i][j] = x[i][j];
			}
				yo[i] = y[i];
		}

		cout << "Current distances: ";
		for (int i = 0; i < dom.Nion -1; i++)
		{
			cout << ions.Rion[i+1] - ions.Rion[i] << ",  " ;
		}

		cout << endl <<  "Current energy: " << y[0] << endl;

		// Sort values
		xySort( x, y, M, N);
	
		// Calculate centroid
		for (int i = 0; i < N; i++ )
			xBar[i] = 0.0;

		for (int i = 0; i < N; i++ )
		{
			for (int j = 0; j < N; j++ )
				xBar[j] = abs( xBar[j] + ( 1.0/N )* x[i][j] );
		}

		// Calculate Reflection
		for (int i = 0; i < N; i++ )
			xRef[i] = abs( xBar[i] + alpha* ( xBar[i] - x[M-1][i] ) );

		// Evaluate reflection point
		fxr = evaluateGround( dom, elc, ions, eVals, xRef);
		iter = iter + 1; cout<< "Function evals: " << iter << endl;

		// Decision time
		if ( fxr >= y[0] && fxr < y[N] )
		{
			for (int i = 0; i < N; i++ )
				x[M-1][i] = xRef[i];

			y[M-1] = fxr;
		}

		else if ( fxr < y[0])
		{
			// If the reflection point is the new best point, try going further --
			
			// Calculate expansion
			for (int i = 0; i < N; i++ )
				xe[i] = abs( xBar[i] + gamma* ( xBar[i] - x[M-1][i] ) );

			double fxe = 0.0;
			fxe = evaluateGround( dom, elc, ions, eVals, xe);
			iter = iter + 1; cout<< "Function evals: " << iter << endl;

			if (fxe < fxr)
			{
				for (int i = 0; i < N; i++ )
					x[M-1][i] = xe[i];
				y[M-1] = fxe;
			}
			else
			{
				for (int i = 0; i < N; i++ )
					x[M-1][i] = xRef[i];
				y[M-1] = fxr;
			}
		}
		else
		{
			// Going outside the domain didnt work, try staying inside

			// Calculate contraction
			for (int i = 0; i < N; i++ )
				xc[i] = abs( xBar[i] + rho* ( xBar[i] - x[M-1][i] ) );
			fxc = evaluateGround( dom, elc, ions, eVals, xc); 
			iter = iter + 1; cout<< "Function evals: " << iter << endl;

			if ( fxc < y[M-1] )
			{
				for (int i = 0; i < N; i++ )
					x[M-1][i] = xc[i];
				y[M-1] = fxc;
			}

			// Nothing worked: reduce the domain keeping the best point, and moving all points closer to this best point
			else
			{
				for (int i = 1; i < M; i++)
				{
					for (int j = 0; j < N; j++)
						x[i][j] = abs( x[i][j] + sigma * ( x[i][j] - x[0][j] ) );
				
					y[i] = evaluateGround( dom, elc, ions, eVals, x[i]);
					iter = iter + 1; cout<< "Function evals: " << iter << endl;
				}
			}
		}
		
	}
	
	evaluateGround( dom, elc, ions, eVals, x[0]);

	cout << "Final positions: " << endl;
	for (int k = 0; k < dom.Nion; k++)
		cout << ions.Rion[k] << "  ";
	cout << endl;

	// Clean up
	delete [] y;
	delete [] yo;
	
	delete [] xBar;
	delete [] xRef;
	delete [] xe;
	delete [] xc;

	for (int k = 0; k < M; k++ )
	{
		delete [] x[k];
		delete [] xo[k];
	}
	delete [] x;
	delete [] xo;
}

void Maxwell::FDTDstep( Domain dom, Schrodinger elc, int r)
{
	double dt = dom.r_dt;

	source( r );

	// 3a. Efield Update
	for (int k = 1; k <NEx-1; k++)
		QHy[k] = bx[k]*QHy[k] - gx[k]*( Hy[k] - Hy[k-1] )/dz; 

	for (int k = 1; k < NEx-1 ; k++)
		Ex[k] = Ex[k] - a2*dt*( ( Hy[k] - Hy[k-1] )/dz + QHy[k] ); 


	// 3b. Hfield Update
	for (int k = 0; k < NHy; k++)
		QEx[k] = by[k]*QEx[k] - gy[k]*( Ex[k+1] - Ex[k] )/dz; 

	for (int k = 0; k < NHy; k++)
		Hy[k] = Hy[k] - a1*dt*( ( Ex[k+1] - Ex[k] )/dz + QEx[k] ); 

	// 3c. A [vector potential] Update
	for (int k = 1; k <NEx-1; k++)
		Aold[k] = A[k];

	for (int k = 1; k <NEx-1; k++)
		A[k] = A[k] - dt*Ex[k];


}

void Maxwell::source( int r )
{

	double sourceVal;
	double tempAmp = pAmp/3.0;
		sourceVal = 0.0;


		if (pulseSetting == "modulatedGauss")
			sourceVal = exp( -1.0*pow(maxT-(double)r*dt,2.0) / pow(width/5.0, 2.0) ) * tempAmp*sin(omega*double(r)*dt);
		 else if(pulseSetting == "sech")
		{
			double inputT = PI*(dt*r - width*2)/width;
			sourceVal = pAmp*sin(omega*double(r)*dt) * 2 / ( exp(inputT) + exp(-(inputT)) );
		}
		else if(pulseSetting == "gauss")
			sourceVal = pAmp*exp( -1.0*pow(maxT-r*dt,2.0) / pow(width/5.0, 2.0) );
		else if(pulseSetting == "CW")
			sourceVal = pAmp*sin(omega*double(r)*dt);
		else if(pulseSetting == "sin2")
			sourceVal = pAmp*sin(omega*double(r)*dt)*pow( sin( width - r*dt), 2.0 );
		else
			cout << "ERROR 666: Unrecognized source condition, please choose from list" << endl;

		Ex[sourcePos] = Ex[sourcePos] + sourceVal;

}

void interpolateA (Domain dom, Schrodinger elc, Maxwell field, int r)
{

	double* Atemp = new double [field.interpN];
	double* Axtemp = new double [field.interpN];

	for ( int k = 0; k < field.interpN; k++ )
	{
		int i = k + field.psiStart;
		Atemp[k] = ( field.A[ i ] + field.Aold[ i ] )/2.0;
	}

	derivative( Axtemp, Atemp, field.dz, field.interpN );

	spline( field.fieldDom, Atemp, elc.zDomPsi, elc.A, field.interpN, dom.N);
	spline( field.fieldDom, Axtemp, elc.zDomPsi, elc.Ax, field.interpN, dom.N);

	delete[] Atemp;
	delete[] Axtemp;
}

void Schrodinger::VxcCalc()
{

	if (VxcSetting == "rho"){
		for (int k = 0; k< N; k++)
		{
			Vxc[k] = 2.0*rho[k];
		}
	}else if (VxcSetting == "rho1/3"){
		for (int k = 0; k< N; k++)
		{
			Vxc[k] = 4.0/3.0*pow(rho[k], 1.0/3.0);
		}
	}else{
	}
}

void Schrodinger::ExcCalc()
{

	if (VxcSetting == "rho"){
		for (int k = 0; k< N; k++)
		{
			Exc[k] = rho[k];
		}
	}else if (VxcSetting == "rho1/3"){
		for (int k = 0; k< N; k++)
		{
			Exc[k] = pow(rho[k], 1.0/3.0);
		}
	}else{
	}
}

void VxCalc( Domain dom, Schrodinger elc, Newton ions)
{
	for (int k = 0; k<dom.N; k++)
	{
		elc.Vx[k] = 0.0;
	}

	double modify;
	if ( elc.pseudo == "exp"){
		for	 (int s = 0; s< dom.Nion; s++)
		{
			for (int k = 0; k< dom.N; k++)
			{
			double zLoc =  (dom.zDomPsi[k] - ions.Rion[s]);
			modify =  ( elc.pseudoCharge + ions.Zion[s] ) * 
				            (1 - exp( -pow (zLoc,2.0) /  pow (elc.coreRadius,2.0) )*elc.coreEffect );
			elc.Vx[k] = elc.Vx[k] - modify*abs( 1/sqrt(pow( elc.gamma,2.0) +
				                                  pow( (dom.zDomPsi[k] - ions.Rion[s]) ,2.0) ) );
			}
		}
	}else{
		for	 (int s = 0; s< dom.Nion; s++)
		{
			for (int k = 0; k< dom.N; k++)
			{
				elc.Vx[k] = elc.Vx[k] - ions.Zion[s]*abs( 1/sqrt(pow( elc.gamma,2.0) + pow( (dom.zDomPsi[k] - ions.Rion[s]) ,2.0) ) );
			}
		}
	}
}