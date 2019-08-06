// finiteDifference.cpp
#include "SN_3D.h"

Derivative::Derivative(double derivative, double accuracy, string boundary, int Nx, int Ny, int Nz, double Lx, double Ly, double Lz)
{

	int cdNum, ext, evals;
	double* o4 = NULL;	double* dax = NULL;	double* day = NULL;	double* daz = NULL;

	// Generate dx,dy,dz for use in sparse matrix generation
	dx = Lx/(Nx-1);
	dy = Ly/(Ny-1);
	dz = Lz/(Nz-1);
	N3 = Nx*Ny*Nz;
	
	// Order of matrix
	if (derivative == 1)
	{


	}else if (derivative == 2){
			
		if (accuracy >= 2.0 && accuracy < 3.0)
		{

			std::cout << "Generating 2nd deriv FD with 2nd order accuracy" ;
			cdNum = 3;
			ext = (cdNum-1) / 2;
			evals = cdNum*3*N3*2;
			o4 = new double [cdNum];

			o4[0] =	 1.0 ;
			o4[1] = -2.0 ;
			o4[2] =  1.0 ;

			dax = new double [cdNum];
			day = new double [cdNum];
			daz = new double [cdNum];
			for (int k = 0; k < cdNum; k++)
			{
				dax[k] = o4[k]/pow(dx,2.0);
				day[k] = o4[k]/pow(dy,2.0);
				daz[k] = o4[k]/pow(dz,2.0);
			}

		}
		else if (accuracy >= 3.0)
		{
			std::cout << "Generating 2nd deriv FD with 4th order accuracy" ;
			cdNum = 5;
			ext = (cdNum-1) / 2;
			evals = cdNum*3*N3*2;
			double* o4 = new double [cdNum];

			o4[0] =	-1.0 / 12.0;
			o4[1] =	 4.0 / 3.0;
			o4[2] = -5.0 / 2.0;
			o4[3] =  4.0 / 3.0;
			o4[4] = -1.0 / 12.0;
			dax = new double [cdNum];
			day = new double [cdNum];
			daz = new double [cdNum];
			for (int k = 0; k < cdNum; k++)
			{
				dax[k] = o4[k]/pow(dx,2.0);
				day[k] = o4[k]/pow(dy,2.0);
				daz[k] = o4[k]/pow(dz,2.0);
			}
		}

	}

	SparseDMat(evals, cdNum, ext, dax, day, daz, boundary, Nx, Ny, Nz, Lx, Ly, Lz);

	delete [] o4; delete [] dax; delete [] day; delete [] daz;
}

void Derivative::SparseDMat(int evals, int cdNum, int ext, double* dax, double* day,  double* daz,  string boundary, int Nx, int Ny, int Nz, double Lx, double Ly, double Lz)
{
	int ** index = new int * [2];
	index[0] = new int [evals];
	index[1] = new int [evals];

	double * D2temp = new double [evals];


	for (int k = 0; k < evals; k++)
	{
		index[0][k] = 0;
		index[1][k] = 0;
		 D2temp[k] = 0.0;
	}

	
	int k = 0;
	int m = 0;
	int n = 0;
	std::cout << " of size Nx = " << Nx << "   \t Ny = " << Ny << "   \t Nz = " << Nz << " " ;

	if ( boundary == "periodic" )
	{
		std::cout << " with periodic boundaries" << endl;

		// Cycle through the rows, i represents the i'th row of the full matrix.
		for (int i = 0; i < N3; i++)
		{
			m = i/Nz;
			n = i/(Nz*Ny);

			// Main diagonal values
			index[0][k] = i;
			index[1][k] = i;
			D2temp [k] = ( dax[ext] + day[ext] + daz[ext] );
			k++;

			// Place z values along main diagaonal
			for (int s = -ext; s <= ext; s++)
			{
				if (i + s - m*Nz >= 0 && i+s - m*Nz < Nz && s != 0)
				{
					index[0][k] = i;
					index[1][k] = i + s;
					D2temp[k] = daz[s+ext];

					k++;
				} else if ( i + s - m*Nz <= 0 && s != 0) {
					index[0][k] = i;
					index[1][k] = i + Nz + s;
					D2temp[k] = daz[s+ext];

					k++;
				} else if ( i + s - m*Nz  >= Nz && s != 0) {
					index[0][k] = i;
					index[1][k] = i - Nz + s;
					D2temp[k] = daz[s+ext];

					k++;
				}
			}

			// Place y values along middle band
			for (int s = -ext; s <= ext; s++)
			{
				if (i + s*Nz - n*Nz*Ny >= 0 && i+s*Nz - n*Nz*Ny < Nz*Ny && s != 0)
				{
					index[0][k] = i;
					index[1][k] = i + s*Nz;
					D2temp[k] = day[s+ext];

					k++;
				} else if ( i + s*Nz - n*Nz*Ny  <= 0 && s != 0) {
					index[0][k] = i;
					index[1][k] = i + Nz*s + Ny*Nz;
					D2temp[k] = day[s+ext];

					k++;
				} else if ( i + s*Nz - n*Nz*Ny >= Nz && s != 0) {
					index[0][k] = i;
					index[1][k] = i + Nz*s - Ny*Nz;
					D2temp[k] = day[s+ext];

					k++;
				}
			}

			// Place x values along outer band
			for (int s = -ext; s <= ext; s++)
			{
				if (i + s*Nz*Ny  >= 0 && i + s*Nz*Ny < N3 && s != 0)
				{
					index[0][k] = i;
					index[1][k] = i + s*Nz*Ny;
					D2temp[k] = dax[s+ext];

					k++;
				} else if ( i + s*Nz*Ny  <= 0 && s != 0) {
					index[0][k] = i;
					index[1][k] = i + s*Nz*Ny + N3;
					D2temp[k] = dax[s+ext];

					k++;
				} else if ( i + s*Nz*Ny >= Nz && s != 0) {
					index[0][k] = i;
					index[1][k] = i + s*Nz*Ny - N3 ;
					D2temp[k] = dax[s+ext];

					k++;
				}
			}
		}
	}else{
		std::cout << " with zero boundaries" << endl;

		// Cycle through the rows, i represents the i'th row of the full matrix.
		for (int i = 0; i < N3; i++)
		{
			m = i/Nz;
			n = i/(Nz*Ny);

			// Main diagonal values
			index[0][k] = i;
			index[1][k] = i;
			D2temp [k] = ( dax[ext] + day[ext] + daz[ext] );
			k++;

			// Place z values along main diagaonal
			for (int s = -ext; s <= ext; s++)
			{
				if (i + s - m*Nz >= 0 && i+s - m*Nz < Nz && s != 0)
				{
					index[0][k] = i;
					index[1][k] = i + s;
					D2temp[k] = daz[s+ext];

					k++;
				}
			}

			// Place y values along middle band
			for (int s = -ext; s <= ext; s++)
			{
				if (i + s*Nz - n*Nz*Ny >= 0 && i+s*Nz - n*Nz*Ny < Nz*Ny && s != 0)
				{
					index[0][k] = i;
					index[1][k] = i + s*Nz;
					D2temp[k] = day[s+ext];

					k++;
				}
			}

			// Place x values along outer band
			for (int s = -ext; s <= ext; s++)
			{
				if (i + s*Nz*Ny  >= 0 && i + s*Nz*Ny < N3 && s != 0)
				{
					index[0][k] = i;
					index[1][k] = i + s*Nz*Ny;
					D2temp[k] = dax[s+ext];

					k++;
				}
			}
		}
	}

	M = k;
	Di = new int * [2];
	Di[0] = new int [M];
	Di[1] = new int [M];

	D = new double [M];

	for (int k = 0; k < M; k++)
	{
		Di[0][k] = index[0][k];
		Di[1][k] = index[1][k];
		D[k] = D2temp[k];
	}

	delete [] D2temp;
	delete [] index[0];
	delete [] index[1];
	delete [] index;

}

// Takes n'th derivative defined in "Derivative", assumes dimension and derivative order as initialized into "Derivative"
void Derivative::deriv( double * U, double* dU ){

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	int k = 0;
	for (int i = 0; i < N3; i++)
	{
		//int k = idex[i];

		double sum = 0;
		while (Di[0][k] == i)
		{
			sum = sum + D[k]*U[Di[1][k]];
			k++;
		}

		dU[i] = sum;
	}

}

// Takes n'th derivative defined in "Derivative", assumes dimension and derivative order as initialized into "Derivative"
void Derivative::derivCX( CX *U, CX *dU ){

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	int k = 0;

	for (int i = 0; i < N3; i++)
	{
		CX sum = 0.0;
		while (Di[0][k] == i)
		{
			sum = sum + D[k]*U[Di[1][k]];
			k++;
		}

		dU[i] = sum;
	}

}

PoissonMG::PoissonMG( string boundary, int NxIn, int NyIn, int NzIn,  double LxIn, double LyIn, double LzIn)
{
	Jweight = 1.39;
	mu = 1;
	Ncyc = 2;
	nu1 = 3;
	nu2 = 1;
	cycle = "FMG";
	
	int minDim = 10000;

	if (NxIn < minDim)
		minDim = NxIn;

	if (NyIn < minDim)
		minDim = NyIn;

	if (NzIn < minDim)
		minDim = NzIn;

	int Ntemp = minDim;
	int i = 0;
	while  ( Ntemp >= 2 )
	{
		i++;
		Ntemp = (Ntemp + 1) / 2  - 1;
	}
	Ngrids = i;
	cout << endl << "Number of grids for Multigrid Poisson solver = " << Ngrids << endl;

	NxFinal = NxIn;
	NyFinal = NyIn;
	NzFinal = NzIn;
			
	indexR = new int[Ngrids];
	indexD = new int[Ngrids];

	Nx = new int [Ngrids];
	Ny = new int [Ngrids];
	Nz = new int [Ngrids];
	
	Lx = new double [Ngrids];
	Ly = new double [Ngrids];
	Lz = new double [Ngrids];

	dx = new double [Ngrids];
	dy = new double [Ngrids];
	dz = new double [Ngrids];

	for (int i = 0; i < Ngrids; i++)
	{
		Nx[i] = (int) ( ( NxFinal+1 )/ pow(2.0, i) ) - 1;
		Ny[i] = (int) ( ( NyFinal+1 )/ pow(2.0, i) ) - 1;
		Nz[i] = (int) ( ( NzFinal+1 )/ pow(2.0, i) ) - 1;

		if (i == 0)
		{
			Lx[i] = LxIn;
			Ly[i] = LyIn;
			Lz[i] = LzIn;

			dx[i] = Lx[i]/(Nx[i] - 1.0);
			dy[i] = Ly[i]/(Ny[i] - 1.0);
			dz[i] = Lz[i]/(Nz[i] - 1.0);

		}else{

			// Reduce size of domain for internal grid iterates
			Lx[i] = Lx[i-1] - 2*dx[i-1];
			Ly[i] = Ly[i-1] - 2*dy[i-1];
			Lz[i] = Lz[i-1] - 2*dz[i-1];

			dx[i] = Lx[i]/(Nx[i] - 1.0);
			dy[i] = Ly[i]/(Ny[i] - 1.0);
			dz[i] = Lz[i]/(Nz[i] - 1.0);

		}

		if (Nx[i] < 2)
			cout << "Warning, Nx < 2 on coarse grid, change Ngrid" << endl;
	}

	// Generate large storage arrays
	int Ntot = 0;
	pIndex = new int [Ngrids];
	for (int i = 0; i < Ngrids; i++)
	{
		pIndex[i] = Ntot;
		Ntot = Ntot + Nx[i]*Ny[i]*Nz[i];
	}

	pf = new double [Ntot];
	pU = new double [Ntot];

	for (int i = 0; i < Ntot; i ++)
	{
		pf[i] = 0.0;
		pU[i] = 0.0;
	}

	// Generate large derivative arrays to represent set of sparse matrices
	int nk = 0;
	double accuracy = 2.0;
	for (int i = 0; i < Ngrids; i++)
		nk = nk + Nx[i] * Ny[i] * Nz[i] * (int ) pow(accuracy,3.0);

	int k = 0;
	int j = 0;

	double* Dtemp = new double [nk];
	int** Ditemp = new int* [2];
	Ditemp[0] = new int [nk];
	Ditemp[1] = new int [nk];

	double* Rtemp= new double [nk];
	int** Ritemp= new int* [2];
	Ritemp[0] = new int [nk];
	Ritemp[1] = new int [nk];


	for (int s = 0; s < Ngrids; s++)
	{
		indexD[s] = j;
		indexR[s] = k;

		Derivative D2( 2.0, accuracy , boundary, Nx[s], Ny[s], Nz[s], Lx[s], Ly[s], Lz[s]);

		for (int m = 0; m < D2.M; m++)
		{
			if (D2.Di[0][m] == D2.Di[1][m])
			{
				Dtemp[j] = D2.D[m];
				Ditemp[0][j] = D2.Di[0][m];
				Ditemp[1][j] = D2.Di[1][m];

				j++;
				
			}else{
				Rtemp[k] = D2.D[m];
				Ritemp[0][k] = D2.Di[0][m];
				Ritemp[1][k] = D2.Di[1][m];
				k++;
			}
			
		}

	}
	
	int L1 = j;
	int L2 = k;

	D = new double [L1];
	Di = new int* [2];
	Di[0] = new int [L1];
	Di[1] = new int [L1];

	indexD[0] = 0;
	int count = 1;
	for (int i = 0; i < L1; i++)
	{
		D[i] = Dtemp[i];
		Di[0][i] = Ditemp[0][i];
		Di[1][i] = Ditemp[1][i];
	}

	R = new double [L2];
	Ri = new int* [2];
	Ri[0] = new int [L2];
	Ri[1] = new int [L2];
	

	indexR[0] = 0;
	count = 1;
	for (int i = 0; i < L2; i++)
	{
		R[i] = Rtemp[i];
		Ri[0][i] = Ritemp[0][i];
		Ri[1][i] = Ritemp[1][i];
	}

	delete [] Ditemp[0];
	delete [] Ditemp[1];
	delete [] Ditemp;
	delete [] Dtemp;

	delete [] Ritemp[0];
	delete [] Ritemp[1];
	delete [] Ritemp;
	delete [] Rtemp;

}

void PoissonMG::PoissonSolve( double* f, double* guess, double* U)
{
	
	int N3 = Nx[0]*Ny[0]*Nz[0];

	// Assign first set of data to general variable array
	for (int k = 0; k < N3; k ++)
	{
			pU[k] = guess[k];
			pf[k] = f[k];
	}

	if (cycle == "mu"){
	for (int i = 0; i < Ncyc; i++)
		muCycle( 0 );

	}else if ( cycle == "CG" ){
		CG(0,N3);
	}
	else{ 
		FMG(0);
	}

	for (int k = 0; k < N3; k ++)
		U[k] = pU[k];
}

void PoissonMG::muCycle( int i )
{
	int N3 = Nx[i]*Ny[i]*Nz[i];

	for (int k = 0; k < nu1; k++)
		gaussSeidel (i, N3);

	if (i != Ngrids - 1)
	{
		
		restrict3D (i);

		for (int k = pIndex[i+1]; k < pIndex[i+2]; k++)
			pU[k] = 0.0;

		for (int j = 0; j < mu; j++)
		{
			muCycle(i+1);
			
		}

		prolongate(i);
		
	}else{
		N3 = Nx[i]*Ny[i]*Nz[i];
		for (int k = 0; k < 5; k++)
			gaussSeidel (i, N3);

	}

	N3 = Nx[i]*Ny[i]*Nz[i];
	for (int k = 0; k < nu2; k++)
		gaussSeidel (i, N3);

}

// Full multigrid cycle, calls muCycle for intermediate stages.
void PoissonMG::FMG( int i )
{
	int N3 = Nx[i]*Ny[i]*Nz[i];

	if (i != Ngrids - 1)
	{

		restrictFMG(i);

		FMG(i + 1);

		prolongateFMG(i);
	
	} else {

		for (int k = pIndex[i]; k < pIndex[i+1]; k++)
			pU[k] = 0.0;

		N3 = Nx[i]*Ny[i]*Nz[i];
		for (int k = 0; k < 10; k++)
			gaussSeidel (i, N3);

	}

	for (int j = 0; j < Ncyc; j++)
	{
		muCycle( i );
	}

}

// Multiples sparse matrix A with vector v and puts the result in u, i.e. u = A*v;
void indexMult ( int N3, int startA, int startV, double* A, int** Ai, double* v, double* u )
{
	int k = startA;
	for (int i = 0; i < N3; i++)
	{
		double sum = 0;
		while (Ai[0][k] == i)
		{
			sum = sum + A[k]*v[Ai[1][k] + startV];
			k++;
		}

		u[i] = sum;
	}


}

// Iterates members of PoissonMG using standard jacobi iterate
void PoissonMG::jacobi (int i, int N3)
{			
	
	N3 = Nx[i]*Ny[i]*Nz[i];
	double* vtemp = new double [N3];

	indexMult (N3, indexR[i], pIndex[i], R, Ri, pU, vtemp );

	int m = indexD[i];
	for (int k = 0; k < N3; k++){
		pU[k + pIndex[i]] = (1-Jweight)*pU[k + pIndex[i]] + Jweight*( pf[k + pIndex[i]] - vtemp[k] ) / D[m];
		m++;
	}

	delete [] vtemp;
}

// Iterates members of PoissonMG using Gauss-seidel iterate [typically works 2x faster than jacobi iterate]
void PoissonMG::gaussSeidel (int i, int N3)
{			
	
	N3 = Nx[i]*Ny[i]*Nz[i];

	double* wtemp = new double [N3];
	for (int k = 0; k < N3; k ++)
		wtemp[k] = pU[k + pIndex[i]];

	int m = indexD[i];
	int startA = indexR[i];
	int startV =  pIndex[i];

	int k = startA;
	double coef1 = 1.0-Jweight;

	for (int n = 0; n < N3; n++)
	{
		double sum = 0.0;
		while (Ri[0][k] == n)
		{
			sum = sum + R[k]*pU[Ri[1][k] + startV];
			k++;
		}
		//Gauss seidel puts new values back into matrix so that in the remainder of the iterations the U_(n+1) values are used
		//whever they have been already updated
		pU[n + pIndex[i]] = coef1*wtemp[n] + Jweight*( pf[n + pIndex[i]] - sum ) / D[m];
		m++;
	}
	
	delete [] wtemp;
}

// Iterates members of PoissonMG using conjugate gradient (CG) method
void PoissonMG::CG (int i, int N3)
{			
	cout << i << endl;
	N3 = Nx[i]*Ny[i]*Nz[i];

	double* r = new double [N3];
	double* p = new double [N3];

	double* tempv = new double [N3];
	double* tempu = new double [N3];

	indexMult( N3, indexD[i], pIndex[i], D, Di, pU, p);
	indexMult( N3, indexR[i], pIndex[i], R, Ri, pU, r);
	
	for (int k = 0; k < N3; k++)
		r[k] = pf[pIndex[i] + k] - (p[k] + r[k] );
		
	for (int k = 0; k < N3; k++)
		p[k] = r[k];

	//-------------------------------------------------------------//

	double* Ap = new double [N3];

	double alpha;
	double rs_new = 1e10;
	double rs_old = 0.0;

	rs_old = 0.0;
	for (int k = 0; k < N3; k++)
		rs_old = rs_old + r[k] * r[k];
	int count = 0;

	int iters;
	

	iters = (int) pow(2,i+1);

	
	//-------------------------------------------------------------//
	while ( count < 1 )
	{
		indexMult( N3, indexD[i], 0, D, Di, p, tempv);
		indexMult( N3, indexR[i], 0, R, Ri, p, tempu);

		for (int k = 0; k < N3; k++)
			Ap[k] = tempv[k] + tempu[k];

		alpha = 0.0;
		for (int k = 0; k < N3; k++)
			alpha = alpha + Ap[k]*p[k];

		alpha = rs_old/alpha;

		for (int k = 0; k < N3; k++)
			pU[pIndex[i] + k] = pU[pIndex[i] + k] + alpha*p[k];

		for (int k = 0; k < N3; k++)
			r[k] = r[k] - alpha*Ap[k];

		rs_new = 0.0;
		for (int k = 0; k < N3; k++)
			rs_new = rs_new + r[k] * r[k];
		
		for (int k = 0; k < N3; k++)
			p[k] = r[k] + rs_new/rs_old * p[k];

		rs_old = rs_new;
		count ++;

	}
	

	cout << "iters : " << count << endl;
	delete [] p;
	delete [] r;
	delete [] Ap;
	delete [] tempv;
	delete [] tempu;
}


// Places the i'th fine grid onto the i+1'th coarser grid 
// (e.g., averages the 127x127x127 grid onto the 63x63x63 grid )
void PoissonMG::restrict3D (int s)
{

	int K_lg = pIndex[s];
	int K_sm = pIndex[s+1];

	int N3_lg = Nx[s]*Ny[s]*Nz[s];
	int N3_sm = Nx[s+1]*Ny[s+1]*Nz[s+1];

	// Something wrong here
	double* resid = new double [N3_lg];
	double* temp = new double[N3_lg];

	indexMult( N3_lg, indexD[s], pIndex[s], D, Di, pU, temp);
	indexMult( N3_lg, indexR[s], pIndex[s], R, Ri, pU, resid);

	for (int k = 0; k < N3_lg; k++)
		resid[k] = pf[K_lg + k] - ( temp[k] + resid[k] );
		
	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	for ( int i = 0; i < Nx[s+1]; i++ )
	{
		for (int j = 0; j < Ny[s+1]; j++)
		{
			for (int k = 0; k < Nz[s+1]; k++)
			{
				int ci = ( 2*i + 1) *Ny[s]*Nz[s] ;
				int cj = ( 2*j + 1) *Nz[s] ;
				int ck = ( 2*k + 1) ;

				pf[ K_sm + i*Ny[s+1]*Nz[s+1] + j*Nz[s+1] + k ] = 
		   1.0/64.0 * ( resid[ ci-1 + cj-1 + ck-1 ] + resid[ ci+1 + cj+1 + ck+1 ] + 
					    resid[ ci-1 + cj-1 + ck+1 ] + resid[ ci-1 + cj+1 + ck-1 ] + 
					    resid[ ci+1 + cj-1 + ck-1 ] + resid[ ci+1 + cj-1 + ck+1 ] + 
					    resid[ ci+1 + cj+1 + ck-1 ] + resid[ ci-1 + cj+1 + ck+1 ] + 

					2*( resid[ ci + cj-1 + ck-1 ] + resid[ ci-1 + cj + ck-1 ] + 
					    resid[ ci-1 + cj-1 + ck ] + resid[ ci + cj+1 + ck+1 ] + 
					    resid[ ci+1 + cj + ck+1 ] + resid[ ci+1 + cj+1 + ck ] + 
					    resid[ ci + cj-1 + ck+1 ] + resid[ ci-1 + cj + ck+1 ] + 
					    resid[ ci-1 + cj+1 + ck ] + resid[ ci + cj+1 + ck-1 ] + 
					    resid[ ci+1 + cj + ck-1 ] + resid[ ci+1 + cj-1 + ck ] ) + 

					4* (resid[ ci + cj + ck-1 ] + resid[ ci + cj + ck+1 ] + 
					    resid[ ci + cj-1 + ck ] + resid[ ci + cj+1 + ck ] + 
					    resid[ ci-1 + cj + ck ] + resid[ ci+1 + cj + ck ] ) + 

					8 * resid[ ci + cj + ck ]);

			}
		}
	}

	delete [] resid;
	delete [] temp;

}

// Places the i'th fine grid onto the i+1'th coarser grid 
// (e.g., averages the 127x127x127 grid onto the 63x63x63 grid )
void PoissonMG::restrictFMG (int s)
{

	int K_lg = pIndex[s];
	int K_sm = pIndex[s+1];

	int N3_lg = Nx[s]*Ny[s]*Nz[s];
	int N3_sm = Nx[s+1]*Ny[s+1]*Nz[s+1];

	// Something wrong here
	double* resid = new double [N3_lg];

	for (int k = 0; k < N3_lg; k++)
		resid[k] = pf[K_lg + k];
		
	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	for ( int i = 0; i < Nx[s+1]; i++ )
	{
		for (int j = 0; j < Ny[s+1]; j++)
		{
			for (int k = 0; k < Nz[s+1]; k++)
			{
				int ci = ( 2*i + 1) *Ny[s]*Nz[s] ;
				int cj = ( 2*j + 1) *Nz[s] ;
				int ck = ( 2*k + 1) ;

				pf[ K_sm + i*Ny[s+1]*Nz[s+1] + j*Nz[s+1] + k ] = 
		   1.0/64.0 * ( resid[ ci-1 + cj-1 + ck-1 ] + resid[ ci+1 + cj+1 + ck+1 ] + 
					    resid[ ci-1 + cj-1 + ck+1 ] + resid[ ci-1 + cj+1 + ck-1 ] + 
					    resid[ ci+1 + cj-1 + ck-1 ] + resid[ ci+1 + cj-1 + ck+1 ] + 
					    resid[ ci+1 + cj+1 + ck-1 ] + resid[ ci-1 + cj+1 + ck+1 ] + 

					2*( resid[ ci + cj-1 + ck-1 ] + resid[ ci-1 + cj + ck-1 ] + 
					    resid[ ci-1 + cj-1 + ck ] + resid[ ci + cj+1 + ck+1 ] + 
					    resid[ ci+1 + cj + ck+1 ] + resid[ ci+1 + cj+1 + ck ] + 
					    resid[ ci + cj-1 + ck+1 ] + resid[ ci-1 + cj + ck+1 ] + 
					    resid[ ci-1 + cj+1 + ck ] + resid[ ci + cj+1 + ck-1 ] + 
					    resid[ ci+1 + cj + ck-1 ] + resid[ ci+1 + cj-1 + ck ] ) + 

					4* (resid[ ci + cj + ck-1 ] + resid[ ci + cj + ck+1 ] + 
					    resid[ ci + cj-1 + ck ] + resid[ ci + cj+1 + ck ] + 
					    resid[ ci-1 + cj + ck ] + resid[ ci+1 + cj + ck ] ) + 

					8 * resid[ ci + cj + ck ]);

			}
		}
	}

	delete [] resid;

}

// Places the i+1'th coarser grid back onto the i'th fine grid
// (e.g., 63x63x63 grid --> 127x127x127 grid ) 
void PoissonMG::prolongate (int s)
{
	int K_sm = pIndex[s+1];
	int K_lg = pIndex[s];
	int N3_lg = Nx[s]*Ny[s]*Nz[s];
	double* v = new double [ N3_lg ];

	// Remove this later
	for (int k = 0; k < N3_lg; k++)
		v[k] = 0.0;

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Same point direct copying
	for ( int i = 1; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 1; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				
				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = pU[ K_sm +  ci + cj + ck];
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in z-direction
	for ( int i = 1; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 1; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 0; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.5*( pU[ K_sm +  ci + cj + ck ] + pU[ K_sm +  ci + cj + ck + 1] );
			}
		}
	}


	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in y-direction
	for ( int i = 1; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;
				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.5* ( pU[ K_sm + ci + cj + ck] + pU[ K_sm +  ci + cj + Nz[s+1] + ck] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in x-direction
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 1; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				
				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.5* ( pU[ K_sm +  ci + cj + ck] + pU[ K_sm + ci + Ny[s+1]*Nz[s+1]+ cj + ck] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in z-y face
	for ( int i = 1; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 0; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;
				
				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.25*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + ck + 1]  + pU[ K_sm +  ci + cj + Nz[s+1] + ck] + pU[ K_sm +  ci + cj + Nz[s+1] + ck + 1] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in x-z face
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 1; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 0; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.25*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + ck + 1]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck + 1] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in x-y face
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.25*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + Nz[s+1] + ck ]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + Nz[s+1] + ck ] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in center of cube
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.25*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + Nz[s+1] + ck ]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + Nz[s+1] + ck ] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in center of cube
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 0; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.125*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + Nz[s+1] + ck ]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + Nz[s+1] + ck ] + pU[ K_sm +  ci + cj + ck +1] 
				+ pU[ K_sm +  ci + cj + Nz[s+1] + ck +1 ]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck + 1] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + Nz[s+1] + ck  + 1] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	for (int k = 0; k < N3_lg; k++)
		pU[K_lg + k] = pU[K_lg + k] + v[k];

	delete [] v;
}

// Places the i+1'th coarser grid back onto the i'th fine grid
// (e.g., 63x63x63 grid --> 127x127x127 grid ) 
void PoissonMG::prolongateFMG (int s)
{
	int K_sm = pIndex[s+1];
	int K_lg = pIndex[s];
	int N3_lg = Nx[s]*Ny[s]*Nz[s];
	double* v = new double [ N3_lg ];

	// Remove this later
	for (int k = 0; k < N3_lg; k++)
		v[k] = 0.0;

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Same point direct copying
	for ( int i = 1; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 1; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				
				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = pU[ K_sm +  ci + cj + ck];
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in z-direction
	for ( int i = 1; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 1; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 0; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.5*( pU[ K_sm +  ci + cj + ck ] + pU[ K_sm +  ci + cj + ck + 1] );
			}
		}
	}


	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in y-direction
	for ( int i = 1; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;
				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.5* ( pU[ K_sm + ci + cj + ck] + pU[ K_sm +  ci + cj + Nz[s+1] + ck] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in x-direction
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 1; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				
				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.5* ( pU[ K_sm +  ci + cj + ck] + pU[ K_sm + ci + Ny[s+1]*Nz[s+1]+ cj + ck] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in z-y face
	for ( int i = 1; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 0; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;
				
				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.25*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + ck + 1]  + pU[ K_sm +  ci + cj + Nz[s+1] + ck] + pU[ K_sm +  ci + cj + Nz[s+1] + ck + 1] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in x-z face
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 1; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 0; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.25*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + ck + 1]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck + 1] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in x-y face
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.25*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + Nz[s+1] + ck ]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + Nz[s+1] + ck ] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in center of cube
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 1; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.25*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + Nz[s+1] + ck ]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + Nz[s+1] + ck ] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	// Averaging in center of cube
	for ( int i = 0; i < Nx[s] - 1; i = i + 2 )
	{
		for (int j = 0; j < Ny[s] - 1; j = j + 2)
		{
			for (int k = 0; k < Nz[s] - 1; k = k + 2)
			{
				int ci = ( (i+1)/2 - 1 ) *Ny[s+1]*Nz[s+1] ;
				int cj = ( (j+1)/2 - 1 ) *Nz[s+1] ;
				int ck = ( (k+1)/2 - 1 )  ;

				v[ i*Nz[s]*Ny[s] +  j*Nz[s] + k ] = 0.125*( pU[ K_sm +  ci + cj + ck ] 
				+ pU[ K_sm +  ci + cj + Nz[s+1] + ck ]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + Nz[s+1] + ck ] + pU[ K_sm +  ci + cj + ck +1] 
				+ pU[ K_sm +  ci + cj + Nz[s+1] + ck +1 ]  + pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + ck + 1] 
				+ pU[ K_sm +  ci + Ny[s+1]*Nz[s+1] + cj + Nz[s+1] + ck  + 1] );
			}
		}
	}

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	for (int k = 0; k < N3_lg; k++)
		pU[K_lg + k] = v[k];

	delete [] v;
}

Derivative1D::Derivative1D(double derivative, double accuracy, string boundary, int Nin, double Lin)
{

	int cdNum, ext, evals;
	double* o4 = NULL;	double* dax = NULL;	double* day = NULL;	double* daz = NULL;

	N = Nin;
	L = Lin;
	// Generate dx,dy,dz for use in sparse matrix generation
	dx = L/((double) N-1);

	// Order of matrix
	if (derivative == 1)
	{


	}else if (derivative == 2){
			
		if (accuracy >= 2.0 && accuracy < 3.0)
		{

			std::cout << "Generating 2nd deriv FD with 2nd order accuracy" ;
			cdNum = 3;
			ext = (cdNum-1) / 2;
			evals = cdNum*3*N*2;
			o4 = new double [cdNum];

			o4[0] =	 1.0 ;
			o4[1] = -2.0 ;
			o4[2] =  1.0 ;

			dax = new double [cdNum];
			day = new double [cdNum];
			daz = new double [cdNum];
			for (int k = 0; k < cdNum; k++)
			{
				dax[k] = o4[k]/pow(dx,2.0);
			
			}

		}
		else if (accuracy >= 3.0)
		{
			std::cout << "Generating 2nd deriv FD with 4th order accuracy" ;
			cdNum = 5;
			ext = (cdNum-1) / 2;
			evals = cdNum*3*N*2;
			double* o4 = new double [cdNum];

			o4[0] =	-1.0 / 12.0;
			o4[1] =	 4.0 / 3.0;
			o4[2] = -5.0 / 2.0;
			o4[3] =  4.0 / 3.0;
			o4[4] = -1.0 / 12.0;
			dax = new double [cdNum];
			day = new double [cdNum];
			daz = new double [cdNum];
			for (int k = 0; k < cdNum; k++)
			{
				dax[k] = o4[k]/pow(dx,2.0);
			}
		}

	}

	SparseDMat(evals, cdNum, ext, dax,boundary, N,  L);

	delete [] o4; delete [] dax; delete [] day; delete [] daz;
}

void Derivative1D::SparseDMat(int evals, int cdNum, int ext, double* dax,  string boundary, int N, double L)
{
	int ** index = new int * [2];
	index[0] = new int [evals];
	index[1] = new int [evals];

	double * D2temp = new double [evals];


	for (int k = 0; k < evals; k++)
	{
		index[0][k] = 0;
		index[1][k] = 0;
		 D2temp[k] = 0.0;
	}

	
	int k = 0;
	std::cout << " of size " << N ;

	if ( boundary == "periodic" )
	{
		std::cout << " with periodic boundaries" << endl;

		// Cycle through the rows, i represents the i'th row of the full matrix.
		for (int i = 0; i < N; i++)
		{
	
			// Main diagonal values
			index[0][k] = i;
			index[1][k] = i;
			D2temp [k] = ( dax[ext] );
			k++;

			// Place z values along main diagaonal
			for (int s = -ext; s <= ext; s++)
			{
				if (i + s >= 0 && i+s < N && s != 0)
				{
					index[0][k] = i;
					index[1][k] = i + s;
					D2temp[k] = dax[s+ext];

					k++;
				} else if ( i + s <= 0 && s != 0) {
					index[0][k] = i;
					index[1][k] = i + s + N;
					D2temp[k] = dax[s+ext];

					k++;
				} else if ( i + s >= N && s != 0) {
					index[0][k] = i ;
					index[1][k] = i - N + s;
					D2temp[k] = dax[s+ext];
					k++;
				}
			}

	
		}
	}else{
		std::cout << " with zero boundaries" << endl;

		// Cycle through the rows, i represents the i'th row of the full matrix.
		for (int i = 0; i < N; i++)
		{

			// Main diagonal values
			index[0][k] = i;
			index[1][k] = i;
			D2temp [k] = dax[ext];
			k++;

			// Place z values along main diagaonal
			for (int s = -ext; s <= ext; s++)
			{
				if (i + s >= 0 && i+s  < N && s != 0)
				{
					index[0][k] = i;
					index[1][k] = i + s;
					D2temp[k] = dax[s+ext];

					k++;
				}
			}
		}
	}

	M = k;
	Di = new int * [2];
	Di[0] = new int [M];
	Di[1] = new int [M];

	D = new double [M];

	for (int k = 0; k < M; k++)
	{
		Di[0][k] = index[0][k];
		Di[1][k] = index[1][k];
		D[k] = D2temp[k];
	}

	delete [] D2temp;
	delete [] index[0];
	delete [] index[1];
	delete [] index;

}

// Takes n'th derivative defined in "Derivative", assumes dimension and derivative order as initialized into "Derivative"
void Derivative1D::deriv( double * U, double* dU ){

	#ifdef USE_OPENMP
		#pragma omp parallel for
	#endif
	int k = 0;
	for (int i = 0; i < N; i++)
	{
		//int k = idex[i];

		double sum = 0;
		while (Di[0][k] == i)
		{
			sum = sum + D[k]*U[Di[1][k]];
			k++;
		}

		dU[i] = sum;
	}

}

// Takes n'th derivative defined in "Derivative", assumes dimension and derivative order as initialized into "Derivative"
void Derivative1D::derivCX( CX *U, CX *dU ){


	int k = 0;
	for (int i = 0; i < N; i++)
	{

		CX sum = 0.0;
		while (Di[0][k] == i)
		{
			sum = sum + D[k]*U[Di[1][k]];

			k++;
		}

		dU[i] = sum;
	}

}
