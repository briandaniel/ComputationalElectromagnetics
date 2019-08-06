/*%% Information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This program computes solutions to three models coupled into a single   %
% solver in a 1D toy setting. The models are:                             %
%                                                                         %
%  I. Classical Maxwells' Equations for light propagation                 %
%  II. Schrodinger Equation for quantum electron evolution                %
%  III. Classical Newton's equation for ion motion                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  I. FDTD solver for maxwells equations                                  %
%-------------------------------------------------------------------------%
%                                                                         %
%      dHy        dEx                                                     % 
%  1. ----- = a1*-----                                                    % 
%      dt         dz                                                      % 
%                                                                         %
%      dEx        dHy                                                     % 
%  2. ----- = a2*----- + a3*Je + a4*Jion                                  % 
%      dt         dz                                                      % 
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%      where                                                              %
%                          dPsi_j      dPsi*_j                            % 
%      Je =  Sum [ Psi*_j -------- +  --------- Psi_j ]                   % 
%                            dz          dz                               % 
%                                                                         %
%      Jion = Sum [ Zn * Vn * delta( z - Rn ) ]                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  II. Crank Nicholson scheme with Strang-Splitting to solve the          %
%      Nonlinear Schrodinger equation for quantum electrons               %
%-------------------------------------------------------------------------%
%                                                                         %  
%  dPsi      d^2Psi                                                       % 
%  ---- = b1*------ + b2*(z*Ex) + b3*F(rho(z)) + b4*G(Rn)                 % 
%   dt        dz^2                                                        % 
%                                                                         % 
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%      where                                                              %
%      rho(x) = Sum |Psi_i(x)|^2                                          %
%                                                                         %
%                    rho( z')                                             %
%      F(rho) = Int ---------- dz'                                        %
%                   | z - z' |                                            %
%                i.e. F is the convolution of rho(x) and u(x)             %
%                                                                         %
%                      Zn                                                 %
%      G(Rn) = Sum ---------- dz'                                         %
%                  | z - Rn |                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  III. Finite difference solution to Newton's equation                   %
%-------------------------------------------------------------------------%
%   The second order differential equation is split into a system of two  %
% first order differential equations.                                     %
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%                                                                         %
%          dRn                                                            %
%     1.  ----- = Vn                                                      %
%           dt                                                            %
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%     where                                                               %
%     Rn = position of nth ion (Rion)                                     %
%     Vn = velocity of nth ion (Vion)                                     %
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%                                                                         %
%          dVn                                                            %
%     2.  ----- = c1*F_Ion + c2*F_e + c3*F_Lorentz                        %
%          dt                                                             %
%                                                                         % 
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%     where                                                               %
%                            Rn - Rm                                      %
%     F_Ion = Zn*Zm* Sum  --------------                                  % 
%                    m~=n | Rn  - Rm |^3                                  %   
%                                                                         %
%                  d         rho(z)                                       %
%     F_Ion = Zn* ---- Int ---------- dz                                  % 
%                  dR      | z - Rn |                                     %   
%                                                                         %
%     F_Lorentz = Zn* ( Ex + Vn x By)                                     %
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

// Point arrangement for EH fields
//  |---PML---|  buffer  |-----EMpts=----|  buffer  |---PML---|    
//   x   o   x   o   x   o   x   o   x   o   x   o   x   o   x
//   E   H   E   H   E   H   E   H   E   H   E   H   E   H   E

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
	#define OMP_NUM_THREADS 8
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

// Include class definitions
#include "model.h"

/* mathFunctions.cpp declarations
OUTPUT                      FUNCTION		( INPUT )																						*/
double						round			( double d )																					;
double						trapz			( double *F, int N, double dz )																	;
CX							CXtrapz			( CX *F, int N, double dz )																		;
double*						linspace		( double a, double b, int N )																	;
double						dot				( double* U, double* V, int N )																	;
double						norm2			( CX* U, int N )																				;
double						mod				(double x, double y )																			;
void						sqNormalize		( CX** PSI, int N, int Npsi, double dz )														;
void 						gramSchmidt		( Domain dom, Schrodinger elc)																	;
void 						CXtridisolve	( CX* PSI, CX* a, CX* b, CX* c, CX* d, CX* e, CX* f, int N )									;
void						tridisolve		( double* PSI, double* a, double* b, double* c, double* d, double* e, double* f, int N )		;
void						fftShift		( double* shiftedg, double* g, int N )															;
void						derivative		( double* dF, double* F, double dz, int N )														;
void						spline			( double* x, double* y, double* xx, double* yy, int N, int M )									;
void						CXderivative2	( CX* d2F, CX* F, double dz, int N )															;
void						derivative2		( double* d2F, double* F, double dz, int N )													;
void						xySort			( double** x, double* y,  int M, int N )														;
double						timeRand		()																								;

/* functions.cpp declarations
OUTPUT                      FUNCTION		( INPUT )																						*/
void						sech			( Domain dom, Schrodinger elc)																	;
void						crank			( Domain dom, Schrodinger elc, CX dt, int s)													;
void						ground			( Domain dom, Schrodinger elc, Newton ions, double imTol )										;
void						densitydz		( CX** PSI, double* rho, int Npsi, int N, double dz )											;
void						density			(CX** PSI, double* rho, int Npsi, int N);
void						strang			( Domain dom, Schrodinger elc, CX dt)															;
void						ionPosition		( double* Rion, double* Vion, double dt, int Nion )												;
void						psiStep			( Domain dom, Schrodinger elc, Newton ions, Maxwell )											;
void						ionVelocity		( Domain dom, Schrodinger elc, Newton ions )													;
double						evaluateGround  ( Domain dom, Schrodinger elc, Newton ions, Energy eVals, double* x)							;
void						simplex			( double* xIn, double variation, double xTol, 
											  double yTol, Domain dom, Schrodinger elc, Newton ions, Energy eVals, int maxEvals)			;
void						interpolateA	( Domain dom, Schrodinger elc, Maxwell field, int r);
void						VxCalc			( Domain dom, Schrodinger elc, Newton ions);


/* utilityFunctions.cpp declarations
OUTPUT                      FUNCTION    ( INPUT )                                                              */
void						printData	( int r, int steps,  double Etotal, double* Rion, int Nion, Energy eVals, double* rho, int N  )		;
#endif // HEADER_FILE