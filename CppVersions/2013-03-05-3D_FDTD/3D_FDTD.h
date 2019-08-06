/*%% Information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Copyright (C) Brian D Hong, 2012                                        %
%                                                                         %
% This program is free software: you can redistribute it and/or modify    %
% it under the terms of the GNU General Public License as published by    %
% the Free Software Foundation, either version 3 of the License, or       %
% (at your option) any later version.                                     %
%                                                                         %
% This program is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of	      %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           %
% GNU General Public License for more details.							  %
%																		  %
% You should have received a copy of the GNU General Public License       %
% along with this program.  If not, see <http://www.gnu.org/licenses/>.   %    
%																	      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef HEADER_FILE
#define HEADER_FILE

//#define USE_OPENMP

// Includes
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <complex>
#include "config.h"

#ifdef USE_OPENMP
	#include <omp.h>
	#define OMP_NUM_THREADS 3
#endif

using namespace std;

// Complex variable type definition "CX"
typedef std::complex<double> CX;

// Standard constant definitions
#define PI	   3.14159265
#define EPSNOT 8.854187817620e-12
#define MUNOT  1.25663706e-6
#define CNOT   2.99792458e8
#define ETANOT 3.767303132465403e02
#define PLANCK 6.62606957e-34
#define E      1.60217646e-19
#define ME     9.10938188e-31
#define MP     1.67262158e-27
#define HBAR   1.0545717253e-34  
#define JEV    6.24150974e18

// Make class definitions
#include "model.h"

/* mathFunctions.cpp declarations
OUTPUT                      FUNCTION    ( INPUT )																							*/
double						round		( double d )																						;
double						trapz		( double *F, int N, double dz )																		;
double*						linspace	( double a, double b, int N )																		;
double						dot			( double* U, double* V, int N )																		;
double						norm2		( CX* U, int N )																					;
double						mod			(double x, double y )																				;
void						sqNormalize	( CX** PSI, int N, int Npsi, double dz )															;
void 						gramSchmidt	( CX** PSI, int N, int Npsi, double dz )															;
void 						CXtridisolve( CX* PSI, CX* a, CX* b, CX* c, CX* d, CX* e, CX* f, int N )										;
void						tridisolve	( double* PSI, double* a, double* b, double* c, double* d, double* e, double* f, int N )			;
void						fftShift	( double* shiftedg, double* g, int N )																;
void						derivative	(double* dF, double* F, double dz, int N)															;
void						spline		(double* x, double* y, double* xx, double* yy, int N, int M)										;
void						CXderivative2 (CX* d2F, CX* F, double dz, int N);
void						derivative2 (double* d2F, double* F, double dz, int N);

/* functions.cpp declarations
OUTPUT                      FUNCTION		( INPUT )																*/
void						HxUpdate		( Domain dom, Hfield Hx, Efield Ey, Efield Ez )							;
void						HyUpdate		( Domain dom, Hfield Hy, Efield Ez, Efield Ex )							;
void						HzUpdate		( Domain dom, Hfield Hz, Efield Ex, Efield Ey )							;

void						ExUpdate		( Domain dom, Efield Ex, Hfield Hz, Hfield Hy )							;
void						EyUpdate		( Domain dom, Efield Ey, Hfield Hx, Hfield Hz )							;
void						EzUpdate		( Domain dom, Efield Ez, Hfield Hy, Hfield Hx )							;

void						sourceCondition	( Source aSource, Domain dom, Efield Ez, Efield Ex, Efield Ey, int r )	;
/* utilityFunctions.cpp declarations
OUTPUT                      FUNCTION    ( INPUT )																	*/
void						writeData		( Efield aEfield, string name )											;

#endif // HEADER_FILE