// finiteDifference.h
#include "SN_3D.h"

#ifndef FD_H
#define FD_H


class Derivative{

  public:
	int M;
	int Nx, Ny, Nz, N3;
	double dx, dy, dz;
	double Lx, Ly, Lz;
		
	double * D;
	int ** Di;
	//int * idex;
	
	string boundary;
	// Initializer, generates sparse derivative matrix [vector]
	Derivative::Derivative(double derivative, double accuracy, string boundary, int Nx, int Ny, int Nz, double Lx, double Ly, double Lz);
	void Derivative::SparseDMat(int evals, int cdNum, int ext, double* dax, double* day,  double* daz,  string boundary, int Nx, int Ny, int Nz, double Lx, double Ly, double Lz);

	// Takes derivative using stored derivative sparse matrix
	void deriv ( double * U , double * dU);
	void Derivative::derivCX( CX *U, CX *dU );

	// Destructor
	Derivative::~Derivative()
	{
		//delete [] D;
		//delete [] Di[0];
		//delete [] Di[1];
		//delete [] Di;
	}


};


class Derivative1D{

  public:
	int M;
	int N;
	double dx;
	double L;
		
	double * D;
	int ** Di;
	//int * idex;
	
	string boundary;
	// Initializer, generates sparse derivative matrix [vector]
	Derivative1D::Derivative1D(double derivative, double accuracy, string boundary, int N, double L);
	void Derivative1D::SparseDMat(int evals, int cdNum, int ext, double* dax,  string boundary, int N, double L);

	// Takes derivative using stored derivative sparse matrix
	void deriv ( double * U , double * dU);
	void Derivative1D::derivCX( CX *U, CX *dU );

	// Destructor
	Derivative1D::~Derivative1D()
	{
		//delete [] D;
		//delete [] Di[0];
		//delete [] Di[1];
		//delete [] Di;
	}


};


class PoissonMG{

  public:
	int Ngrids; // Governs entire algorithm by specifying the number of
				// grids to be used in the multigrid setup

	int NxFinal, NyFinal, NzFinal, N3Final;
	double* dx,* dy, *dz;
	string cycle;

	double* Lx,* Ly,* Lz;
	int *Nx, *Ny, *Nz;

	int mu;
	double Jweight;
	int Ncyc;
	
	int nu1;			// Number of jacobi iterations
	int nu2;
	

	double * D; // Diagonal vector containing all diagonals of matrices
	double * R; // Offdiagonal vector containing all off diagonals of matrices
	int ** Ri;
	int ** Di;
	int * indexD;
	int * indexR;

	double* pf;
	double* pU;
	int* pIndex;

	string boundary;

	// Initializer
	PoissonMG::PoissonMG( string boundary, int NxIn, int NyIn, int NzIn, double Lx, double Ly, double Lz);

	// solves laplac[U] = f; guess is starting guess vector.
	void PoissonMG::PoissonSolve( double* f, double* guess, double* U);	
 
private:
	// Solution vectors
	void PoissonMG::muCycle( int i );
	void PoissonMG::jacobi (int i, int N3);
	void PoissonMG::restrict3D (int s);
	void PoissonMG::prolongate (int s);
	void PoissonMG::FMG( int i );
	void PoissonMG::restrictFMG (int s);
	void PoissonMG::prolongateFMG (int s);
	void PoissonMG::gaussSeidel (int i, int N3);
	void PoissonMG::CG (int i, int N3);
};



#endif

