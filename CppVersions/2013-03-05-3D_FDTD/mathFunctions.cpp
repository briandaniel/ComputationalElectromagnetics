#include "3D_FDTD.h"

/* FUNCTION LIST 
1. round
2. trapz
3. linspace
4. dot
5. norm2
6. sqNormalize
7. product
8. CXproduct
9. gramSchmidt
10. tridisolve
*/

// "round": Rounds to nearest integer
double round(double d)
{
  return floor(d + 0.5);
}


// "trapz": Estimates integral on evenly spaced grid by trapezoidal rule
double trapz(double *F, int N, double dz)
{
	double sum = 0;
	for (int k = 0; k < (N-1); k++)
	{
		sum = sum + (F[k] + F[k+1])/2;
	}
	double summed = dz*sum;
	return summed;
}


// Generates array [x_1 ... x_N] with N evenly spaced points with the  x_1 = a, x_N = b
double* linspace( double a, double b, int N )
{
	double* x = new double [N];
	double dx = (b - a)/(N-1);
	for (int k = 0; k < N; k++)
	{
		x[k] = a + dx*k;
	}
	return x;
}

// Computes dot product of vectors U,V both of length N
double dot( double* U, double* V, int N)
{
	double sum = 0;
	for (int k = 0; k<N; k++)
	{
		sum = sum + U[k]*V[k];
	}
	return sum;
}

// Computes the 2-norm of the complex vector U of length N
double norm2 ( CX* U, int N)
{
	double sum = 0;
	for (int k = 0; k<N; k++)
	{
		sum = sum + pow( abs(U[k]), 2.0 );
	}
	return sqrt(sum);
}

// Computes the modulus
double mod (double x, double y)
{
	int n = (int) (x/y);
	double z = x - n*y;
	return z;
}

// Normalizes a function using the integral of |PSI|^2 dz as the norm; 
// i.e. after normalization integral |Psi|^2 dz = 1.
void sqNormalize( CX** PSI, int N, int Npsi, double dz )
{
	for (int s = 0; s< Npsi; s++)
	{
		// Calculate integral
		double PSIint = 0;
		for( int k = 0; k < N-1; k++)
		{
			PSIint = PSIint + dz*( pow( abs( PSI[s][k]), 2.0 ) +  pow( abs( PSI[s][k+1]), 2.0 )  )/2 ; 
		}
	
		// Divide by integral value
		for( int k = 0; k < N; k++)
		{
			PSI[s][k] = PSI[s][k]/sqrt(PSIint);
		}
	}
}

// Gram-Schmidt Orthogonalization
void gramSchmidt( CX** PSI, int N, int Npsi, double dz )
{
	
	CX* v = new CX [N];

	for (int s = 0; s<Npsi; s++)
	{

		double vNorm = norm2 (PSI[s],N);
		CX R = CX ( vNorm , 0.0);
		
		for (int n = 0; n < N; n++)
		{
			v[n] = PSI[s][n] / R;
		}

		for (int k = s+1; k < Npsi; k++)
		{
			R = CX(0.0,0.0);
			for (int n = 0; n < N; n++)
			{
				R = R + PSI[k][n]*v[n];
			}
			for (int n = 0; n < N; n++)
			{
				PSI[k][n] = PSI[k][n] - R*v[n];
			}
		}
	}

	delete [] v;
}

// Solves a complex tridiagonal system using gaussian elimination
// [a b c]*PSInew = [d e f]*PSIold; where the first (a, d) are subdiagonal, and third (c, f) are superdiagonal
void CXtridisolve (CX* PSI, CX* a, CX* b, CX* c, CX* d, CX* e, CX* f, int N)
{
	CX* RHS = new CX [N];

	// Compute RHS
	RHS[0]   = e[0]*PSI[0] + f[0]*PSI[1];
	RHS[N-1] = d[N-2]*PSI[N-2] + e[N-1]*PSI[N-1];
	for (int k = 1; k < (N-1); k++)
	{
		RHS[k] = d[k-1]*PSI[k-1] + e[k]*PSI[k] + f[k]*PSI[k+1];
	}

	for (int k = 0; k < N-1; k++ )
	{
		CX mu = a[k]/b[k];
		b[k+1] = b[k+1] - mu*c[k];
		RHS[k+1] = RHS[k+1] - mu*RHS[k];
	}

	RHS[N-1] = RHS[N-1]/b[N-1];

	for (int M = N-2; M > -1; M-- )
	{
		RHS[M] = ( RHS[M] - c[M]*RHS[M+1] )/ b[M];
	}

	for (int k = 0; k < N; k++ )
	{
		PSI[k] = RHS[k];
	}

	delete [] RHS;
}

// Solves a tridiagonal system using gaussian elimination
// [a b c]*x = d; where the first a is sub diagonal
void tridisolve (double* x, double* a, double* b, double* c, double* d, int N)
{
	for (int k = 0; k < N; k++ )
	{
		x[k] = d[k];
	}

	for (int k = 0; k < N-1; k++ )
	{
		double mu = a[k]/b[k];
		b[k+1] = b[k+1] - mu*c[k];
		x[k+1] = x[k+1] - mu*x[k];
	}

	x[N-1] = x[N-1]/b[N-1];

	for (int M = N-2; M > -1; M-- )
	{
		x[M] = (x[M] - c[M]*x[M+1] )/ b[M];
	}

}

// Moves second half of vector into first half, i.e. 
//    {a_1, a_2, .... a_n} ----> {a_(N/2), a_(N/2+1), ... a_1, a_2, ... a_(N/2-1)}
void fftShift (double* shiftedg, double* g, int N)
{
	for (int k = 0; k<N/2; k++)
		shiftedg[k] = g[N/2 + k];

	for (int k = N/2; k<N; k++)
		shiftedg[k] = g[k - N/2];
}

// Takes derivative of second vector "F" and records into first vector "dF" using periodic boundary
void derivative (double* dF, double* F, double dz, int N)
{
	dF[0]	= ( F[1] - F[N-1] )/ (2.0*dz);
	dF[N-1] = ( F[0] - F[N-2] )/ (2.0*dz);

	for (int k = 1; k < (N-1); k++ )
	{
		dF[k] = (F[k+1] - F[k-1] ) / ( 2.0*dz );
	}
}

// Takes second derivative of second vector "F" and records into first vector "d2F" using periodic boundary
void derivative2 (double* d2F, double* F, double dz, int N)
{
	d2F[0]	 = ( F[1] - 2.0*F[0] + F[N-1] ) / pow(dz,2.0);
	d2F[N-1] = ( F[0] - 2.0*F[N] + F[N-2] ) / pow(dz,2.0);

	for (int k = 1; k < (N-1); k++ )
	{
		d2F[k] = (F[k+1] - 2.0* F[k] + F[k-1] ) / pow(dz,2.0);
	}
}

// Takes second derivative of second vector "F" and records into first vector "d2F" using periodic boundary
void CXderivative2 (CX* d2F, CX* F, double dz, int N)
{
	d2F[0]	 = ( F[1] - 2.0*F[0] + F[N-1] ) / pow(dz,2.0);
	d2F[N-1] = ( F[0] - 2.0*F[N] + F[N-2] ) / pow(dz,2.0);

	for (int k = 1; k < (N-1); k++ )
	{
		d2F[k] = (F[k+1] - 2.0* F[k] + F[k-1] ) / pow(dz,2.0);
	}
}

// Interpolates the function f({x}) = {y} onto {xx} yielding {yy}; 
// where      length(x) = length(y) = N      and      length(xx) = M
void spline (double* x, double* y, double* xx, double* yy, int N, int M)
{
	double* d		= new double [N];
	double* delta	= new double [N-1];
	double* h		= new double [N-1];
	double* a		= new double [N-1];
	double* b		= new double [N];
	double* c		= new double [N-1];
	double* r		= new double [N];
	double* b2		= new double [N];
	double* c2		= new double [N];

	
	for ( int k = 0; k < (N-1) ; k++ )
	{
		h[k] = x[k+1] - x[k];
		delta[k] = ( y[k+1] - y[k] ) / h[k] ;
	}

	// Compute LHS
	for (int k = 1; k < (N-1); k++)
	{
		b[k] = 2*(h[k-1] + h[k]);
	}
	
	for (int k = 1; k < (N-2); k++)
	{
		a[k] = h[k-1];
		c[k] = h[k+1];	
	}	
	
	a[0] = h[1];
	b[0] = h[1];
	c[0] = h[0] + h[1];

	a[N-2] = h[N-2] + h[N-3];
	b[N-1] = h[N-3];
	c[N-2] = h[N-3];
	
	// Compute RHS of tridiagonal system
	for (int k = 1; k < (N-1) ; k++)
		r[k] = 3.0*( h[k]*delta[k-1] + h[k-1]*delta[k] );

	// Semi-erroneous boundary (dont spline at boundary without fixing)
	r[0] = ((h[0]+2*c[0])*h[1]*delta[0] + pow(h[0], 2.0)*delta[1])/c[0];
	r[N-1] = (pow(h[N-2],2.0)*delta[N-3] + (2*a[N-2] + h[N-2]) * h[N-3]*delta[N-2])/a[N-2];
	

	tridisolve(d, a, b, c, r, N);
	
	for ( int k = 0; k<(N-1) ; k++ )
	{
		b2[k] = ( d[k] - 2.0*delta[k] + d[k+1] )/pow(h[k],2.0);
		c2[k] = ( 3.0*delta[k] - 2*d[k] - d[k+1] )/h[k];
	}
	
	int* kVals = new int [M];
	for (int s = 0; s< N; s++)
	{
		for (int i = 0; i<M; i++)
		{
			if ( x[s] <= xx[i] )
			{
				kVals[i] = s;
			}
		}
	}

	
	for (int i = 0; i<M; i++)
	{
		double s = xx[i] - x[kVals[i]];
		yy[i] = y[kVals[i]] + s*( d[kVals[i]] + s*(c2[kVals[i]] + s*b2[kVals[i]] ) ) ;
	}
	
	delete [] d;
	delete [] h;
	delete [] delta;
	delete [] kVals;
	delete [] a;
	delete [] b;
	delete [] c;
	delete [] r;
	delete [] b2;
	delete [] c2;

}

//double nelderMead (