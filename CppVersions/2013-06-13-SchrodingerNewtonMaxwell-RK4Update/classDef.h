// finiteDifference.h
#include "SN_3D.h"

#ifndef CLASS_H
#define CLASS_H

class Domain{

  public:

	int Nx, Ny, Nz, N3; 
	double dx, dy, dz;
	double Lx, Ly, Lz;
	
	int Nt_im, Nt_re; 
	double dt_im, tf_im;
	double dt_re, tf_re;

	string boundary;

	// Initializer, generates sparse derivative matrix [vector]
	Domain();

};

class Schrodinger{

public:
	CX* U;
	CX c1, c2, c3;
	double alpha, beta;
	double *Vext, *rho;
	double dx, dy, dz;
	double xo, yo, zo;
	double Lx, Ly, Lz;
	int Nx, Ny, Nz, N3;


	Schrodinger::Schrodinger(Domain dom);
	void Schrodinger::density();
};


class Manager{

public:
	string videoBin;
	int frames, Nt_re, Nt_im, N3;
	double elapTime;
	clock_t beginT, endT;

	Manager(Domain dom);
	void Manager::vidWrite( int r, double* rho, CX* U, CX dt );
	void Manager::initalData( Domain dom, Schrodinger elc, Derivative D2 );
	void Manager::startTime();
	void Manager::midTime();
	void Manager::endTime();

};



#endif
