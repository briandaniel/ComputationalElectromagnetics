
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Info:                                                           %
% Simulates wave propagation in 1D using the FDTD method which produces   %
% second order accuracy on the traditional YEE staggered grid.            %
%                                                                         %
% 1. Input file provides settings which are imported into the program     %
% 2. Utilizes Ex and Hy fields propagating in z direction                 %
% 3. CPML utilized to truncate spatial domain                             %
% 4. Source condition can be varied in the settings (e.g. SECH, Gaussian) %
% 5. Various Media can be implemented including                           %
%    I. Linear:                                                           %
%    	i.   Dielectric                                                   %
%	    ii.  Lorentz Poles                                                %
%		iii. Drude Poles                                                  %
%		iv.  Debye Poles                                                  %
%    II. Nonlinear                                                        %
%       i.   Kerr                                                         %
%       ii.  Raman                                                        %
%       iii. Multiple 2LVL Atom Quantum Model                             %
% 6. Frequency domain is recorded at intial and final points, and         %
%    reflectivity/transmission data is generated                          %
% 7. Data output is recorded in binary files which can be read by a       %
%    corresponding matlab code.                                           %
%                                                                         %
%                                                                         %
% Brian Hong, ACMS, 2012											      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#pragma region //%%%%%%%%%%%%%%%%%% LIBRARIES %%%%%%%%%%%%%%%%%%//
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "config.h"
using namespace std;
#pragma endregion //%%%%%%%%%%%%%%%% END LIBRARIES %%%%%%%%%%%%%%%%//
#pragma region //%%%%%%%%%%%%%%%% STD CONSTANTS %%%%%%%%%%%%%%%%//
#define PI	   3.14159265
#define EPSNOT 8.854187817620e-12
#define MUNOT  1.25663706e-6
#define CNOT   2.99792458e8
#define ETANOT 3.767303132465403e02
#define PLANCK 6.62606957e-34
#define E      1.60217646e-19
#define ME     9.10938188e-31
#define DIRAC  1.0545717253e-34  
#pragma endregion //%%%%%%%%%%%%%%%% END CONSTANTS %%%%%%%%%%%%%%%%//
#pragma region //%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%//
double round(double d)
{
  return floor(d + 0.5);
}
#pragma endregion //%%%%%%%%%%%%%%%% END FUNCTIONS %%%%%%%%%%%%%%%%//

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
int main () {

	#pragma region //%%%%%%%%%%%%%%% IMPORT SETTINGS %%%%%%%%%%%%%%%//
	int detInit1, detFin1, gridSize, maxIter, sourcePos, detInit, detFin, N, PML, FNUM, PARAM, medSect;
	string fileName, pulseSetting;
	double spaceSize, timeDom, f, pAmp, width, freqMin, freqMax, keepFrames, sigx, refractiveIndex;
	double pmlWidth, sigM, kappaMax, kappaM, aMax, aM, maxC;

	Config settings("FDTDsettings.txt");
	//settings.printSettings(std::cout);
	cout << "" << endl;
	cout << "" << endl;

    settings.get("Domain/gridSize", gridSize);
	N = gridSize;
	settings.get("Domain/spaceSize", spaceSize);
	settings.get("Domain/timeDomain", timeDom);

	settings.get("InputPulse/frequency", f);
	settings.get("InputPulse/amplitude", pAmp);
	settings.get("InputPulse/width", width);
	settings.get("InputPulse/position", sourcePos);
	settings.get("InputPulse/pulse", pulseSetting);

	settings.get("FrequencyDomain/frequencyNumber", FNUM);
	settings.get("FrequencyDomain/minimumFrequency", freqMin);
	settings.get("FrequencyDomain/maximumFrequency", freqMax);
	settings.get("FrequencyDomain/initialDetector", detInit1);
	detInit = sourcePos + detInit1;
	settings.get("FrequencyDomain/finalDetector", detFin1);
	detFin  = N + detFin1;

	settings.get("SimpleMedium/sigma", sigx);
	settings.get("SimpleMedium/mediumSections", medSect);
	settings.get("SimpleMedium/refractiveIndex", refractiveIndex);

	settings.get("PML/cells", pmlWidth);
	PML = (int) pmlWidth;
	settings.get("PML/sigmaExponent", sigM);
    settings.get("PML/sigmaMultiplier", maxC);
	settings.get("PML/kappaMax", kappaMax);
	settings.get("PML/aMax", aMax);
	settings.get("PML/kappaExponent", kappaM);
	settings.get("PML/aExponent", aM);

	settings.get("Output/fileName", fileName);
	settings.get("Output/frames", keepFrames);
	settings.get("Output/outputParameters", PARAM);

	int lorentzPoles, drudePoles, debyePoles;
	settings.get("LorentzMedium/poleNumber", lorentzPoles);
	settings.get("DrudeMedium/poleNumber", drudePoles);
	settings.get("DebyeMedium/poleNumber", debyePoles);

	string debyeOn, drudeOn, lorentzOn, kerrOn, ramanOn, blochOn, dielectricOn;
	settings.get("LorentzMedium/LorentzOn", lorentzOn);
	settings.get("DrudeMedium/DrudeOn", drudeOn);
	settings.get("DebyeMedium/DebyeOn", debyeOn);
	settings.get("2LVLAtom/2LVLAtomOn", blochOn);
	settings.get("KerrMedium/KerrOn", kerrOn);
	settings.get("RamanMedium/RamanOn", ramanOn);
	settings.get("DielectricMedium/DielectricOn", dielectricOn);

	#pragma endregion //%%%%%%%%%%%%% END IMPORT SETTINGS %%%%%%%%%%%%%//	
	#pragma region //%%%%%%%%%%%%%% INITIALIZE DOMAIN %%%%%%%%%%%%%%//
	// Initialize variables and constants
	double dx, dt, sigMax;											// Domain Variables
	double muy;											            // Coef. constants
	int n, k, m, T, on;   									        // Counters
	double source, maxT, lambda, omega, timeSwitch;                 // Source constants
	double curlE, curlH;											// Calculation variables
	double testFreq;    											// Frequency domain variables
	double xVec, sigVec, sigVecStag, kappaVec, aVec, aVecStag;		// PML calculation constants
	double kappaVecStag, tau, tauStag;							    // PML calculation constants
	int skipSteps;

	// MEDIUM ARRAYS	
	double* refIndx = new double [medSect];
	double* epsMed  = new double [medSect];						    
	int** medPos = new int* [medSect];
	for (k = 0; k < medSect; k++)
	    medPos[k] = new int [2];

	// PML ARRAYS
	double* bx = new double [PML];
	double* by = new double [PML];
	double* gx = new double [PML];
	double* gy = new double [PML];
	double* kx = new double [PML];
	double* ky = new double [PML];
	double* QEx = new double [PML*2];
	double* QHy = new double [PML*2];

	// CONSTANT ARRAYS
	double* c1 = new double [N];							
	double* c2 = new double [N];
	double* d1 = new double [N];

	// FREQUENCY DOMAIN ARRAYS
	double* realEref  = new double [FNUM];
	double* realEinc  = new double [FNUM];
	double* realEfin  = new double [FNUM];
	double* imagEref  = new double [FNUM];
	double* imagEinc  = new double [FNUM];
	double* imagEfin  = new double [FNUM];
	double* ampEref   = new double [FNUM];
	double* ampEinc   = new double [FNUM];
	double* ampEfin   = new double [FNUM];
	double* phaseEref = new double [FNUM];
	double* phaseEinc = new double [FNUM];
	double* phaseEfin = new double [FNUM];
	double* transmission = new double [FNUM];
	double* reflection   = new double [FNUM];
	for (k = 0; k<FNUM; k++)
	{
		realEref[k] = 0;
		realEinc[k] = 0;
		realEfin[k] = 0;
		imagEref[k] = 0;
		imagEinc[k] = 0;
		imagEfin[k] = 0;
	}
	for( k = 0; k < PML*2; k++)
	{
		QEx[k] = 0;
		QHy[k] = 0;
	}
	

	// SOLUTION ARRAYS
	double* Ex      = new double [N];
	double* Dx      = new double [N];
	double* Exlin   = new double [N];
	double* ExStore = new double [N];
	double* ExNL    = new double [N];
	double* PL      = new double [N]; 
	double* Jx      = new double [N];
 	double* epsx    = new double [N];
	for( k = 0; k < N; k++)
	{
		Ex[k]      = 0;
		Exlin[k]   = 0;
		ExNL[k]    = 0;
		ExStore[k] = 0;
		Dx[k]      = 0;
		PL[k]      = 0;
		Jx[k]      = 0;
		epsx[k]    = 1;
	}

	double* Hy      = new double [N-1];	
	for( k = 0; k < N-1; k++)
		Hy[k] = 0;

	double* parameterValues = new double[PARAM];
	for (k = 0; k<PARAM; k++)
		parameterValues[k] = 0.0;
	parameterValues[0] = PARAM;

	#pragma endregion //%%%%%%%%%%%% END INITIALIZE DOMAIN %%%%%%%%%%%%//
	#pragma region //%%%%%%%%%%%%%% INITIALIZE MEDIA %%%%%%%%%%%%%%%//

		#pragma region //************** DIELECTRIC **************//
		double dielStart, dielEnd;
		double epsRel;

		if (dielectricOn == "on")
		{
			settings.get("DielectricMedium/startPosition", dielStart);
			settings.get("DielectricMedium/endPosition", dielEnd);
			settings.get("DielectricMedium/epsRel", epsRel);
		}
		#pragma endregion
		#pragma region //*************** LORENTZ ****************//
		double lorzStart, lorzEnd;
		double epsInfLorentz;
		double *epsLorentz = NULL, *fLorentz = NULL, *delLorentz = NULL, *alphaLorz = NULL, *betaLorz = NULL, *gammaLorz = NULL;
		double **Plorz = NULL, **PlorzOld = NULL, **PlorzOld2 = NULL;

		if (lorentzOn == "on")
		{
			cout << "Lorentz Medium on" << endl;

			settings.get("LorentzMedium/startPosition", lorzStart);
			settings.get("LorentzMedium/endPosition", lorzEnd);
			settings.get("LorentzMedium/epsInf", epsInfLorentz);
			epsLorentz    = new double [lorentzPoles];
			fLorentz      = new double [lorentzPoles];
			delLorentz    = new double [lorentzPoles];
			for (k = 0; k<lorentzPoles; k++)
			{
				ostringstream integer;
				integer << k+1;
				string integer2 = integer.str();
				settings.get("LorentzMedium/eps" + integer2, epsLorentz[k]);
				settings.get("LorentzMedium/f" + integer2, fLorentz[k]);
				settings.get("LorentzMedium/del" + integer2, delLorentz[k]);	
			}

			alphaLorz = new double[lorentzPoles]; 
			betaLorz  = new double[lorentzPoles]; 
			gammaLorz = new double[lorentzPoles];
			for (int s = 0; s<lorentzPoles; s++)
			{
				alphaLorz[s] = 0;
				betaLorz[s] = 0;
				gammaLorz[s] = 0;
			}
			Plorz = new double* [N];
			for (k = 0; k< N; k++)
				Plorz[k] = new double [lorentzPoles];
			PlorzOld = new double* [N];
			for (k = 0; k< N; k++)
				PlorzOld[k] = new double [lorentzPoles];
			PlorzOld2= new double* [N];
			for (k = 0; k< N; k++)
				PlorzOld2[k] = new double [lorentzPoles];

			for( k = 0; k < N; k++)
			{
				for (int s = 0; s<lorentzPoles; s++)
				{
					Plorz[k][s] = 0;
					PlorzOld[k][s] = 0;
					PlorzOld2[k][s] = 0;
				}
			}
		}
		#pragma endregion
		#pragma region//**************** DEBYE *****************//
		double debyeStart, debyeEnd;
		double epsInfDebye;
		double *alphaDebye = NULL, *betaDebye = NULL, *tauDebye = NULL, *epsDebye = NULL;
		double **Pdebye = NULL, **PdebyeOld = NULL, **PdebyeOld2 = NULL;
	
		if (debyeOn == "on")
		{
			cout << "Debye Medium on" << endl;
			settings.get("DebyeMedium/startPosition", debyeStart);
			settings.get("DebyeMedium/endPosition", debyeEnd);
			settings.get("DebyeMedium/epsInf", epsInfDebye);

			tauDebye    = new double [debyePoles];
			epsDebye    = new double [debyePoles];

			for (k = 0; k<debyePoles; k++)
			{
				ostringstream integer;
				integer << k+1;
				string integer2 = integer.str();
				settings.get("DebyeMedium/tau" + integer2, tauDebye[k]);
				settings.get("DebyeMedium/eps" + integer2, epsDebye[k]);
		
			}

			alphaDebye= new double[debyePoles]; 
			betaDebye = new double[debyePoles]; 
			for (int s = 0; s<debyePoles; s++)
			{
				alphaDebye[s] = 0;
				betaDebye[s] = 0;
			}
			Pdebye = new double* [N];
			for (k = 0; k< N; k++)
				Pdebye[k] = new double [debyePoles];
			PdebyeOld = new double* [N];
			for (k = 0; k< N; k++)
				PdebyeOld[k] = new double [debyePoles];
			PdebyeOld2 = new double* [N];
			for (k = 0; k< N; k++)
				PdebyeOld2[k] = new double [debyePoles];

			for( k = 0; k < N; k++)
			{
				for (int s = 0; s<debyePoles; s++)
				{
					Pdebye[k][s] = 0;
					PdebyeOld[k][s] = 0;
					PdebyeOld2[k][s] = 0;
				}
			}

		}
		#pragma endregion
		#pragma region//**************** DRUDE *****************//
		double drudeStart, drudeEnd;
		double *fDrude = NULL, *alphaDrude = NULL, *betaDrude = NULL, *kappaDrude = NULL, *gammaDrude = NULL;
		double **Pdrude = NULL, **PdrudeOld = NULL, **PdrudeOld2 = NULL;
	
		if (drudeOn == "on")
		{
			cout << "Drude Medium on" << endl;
			settings.get("DrudeMedium/startPosition", drudeStart);
			settings.get("DrudeMedium/endPosition", drudeEnd);
			fDrude   = new double [drudePoles];
			gammaDrude = new double [drudePoles];

			for (k = 0; k<drudePoles; k++)
			{
				ostringstream integer;
				integer << k+1;
				string integer2 = integer.str();
				settings.get("DrudeMedium/gamma" + integer2, gammaDrude[k]);
				settings.get("DrudeMedium/f" + integer2, fDrude[k]);
			}

			alphaDrude= new double[drudePoles]; 
			betaDrude  = new double[drudePoles]; 
			kappaDrude = new double[drudePoles];
			for (int s = 0; s<drudePoles; s++)
			{
				alphaDrude[s] = 0;
				betaDrude[s] = 0;
				kappaDrude[s] = 0;
			}
			Pdrude = new double* [N];
			for (k = 0; k< N; k++)
				Pdrude[k] = new double [drudePoles];
			PdrudeOld = new double* [N];
			for (k = 0; k< N; k++)
				PdrudeOld[k] = new double [drudePoles];
			PdrudeOld2= new double* [N];
			for (k = 0; k< N; k++)
				PdrudeOld2[k] = new double [drudePoles];

			for( k = 0; k < N; k++)
			{
				for (int s = 0; s<drudePoles; s++)
				{
					Pdrude[k][s] = 0;
					PdrudeOld[k][s] = 0;
					PdrudeOld2[k][s] = 0;
				}
			}
		}
		#pragma endregion
		#pragma region//*********** 2LVL ATOM MODEL ************//
		double blochStart, blochEnd, blochTau1, blochTau2, TNatom, gammaBloch;
		double rho30, blochFreqMin, blochFreqMax, blochCm, blochCp, blochD;
		int odeNumber;
		double *omegaBloch = NULL, *Natom = NULL, *blochA = NULL, *blochB = NULL;
		double **u1 = NULL, **u2  = NULL, **u3 = NULL;
		double **u1store = NULL;

		if (blochOn == "on")
		{
			cout << "2LVL Atom Medium on" << endl;
			settings.get("2LVLAtom/startPosition", blochStart);
			settings.get("2LVLAtom/endPosition", blochEnd);
			settings.get("2LVLAtom/ODEnumber", odeNumber);
			settings.get("2LVLAtom/tau1", blochTau1);
			settings.get("2LVLAtom/tau2", blochTau2);
			settings.get("2LVLAtom/Natom", TNatom);
			settings.get("2LVLAtom/gammaN", gammaBloch);
			settings.get("2LVLAtom/minimumFrequency", blochFreqMin);
			settings.get("2LVLAtom/maximumFrequency", blochFreqMax);
			settings.get("2LVLAtom/rho30", rho30);

			Natom  = new double [odeNumber];
			blochA = new double [odeNumber];
			blochB = new double [odeNumber];
			omegaBloch = new double [odeNumber];
			for (int s = 0; s < odeNumber; s++)
			{
				omegaBloch[s] = 0;
				Natom[s] = 0;
				blochA[s] = 0;
				blochB[s] = 0;
			}

			u1      = new double* [N];
			u2      = new double* [N];
			u3      = new double* [N];
			u1store = new double* [N];

			for (k = 0; k< N; k++)
			{
				u1[k] = new double [odeNumber];
				u2[k] = new double [odeNumber];
				u3[k] = new double [odeNumber];
				u1store [k] = new double [odeNumber];

			}
			for( k = 0; k < N; k++)
			{
				for (int s = 0; s < odeNumber; s++)
				{
				u1[k][s] = 0;
				u2[k][s] = 0;
				u3[k][s] = 0;
				u1store [k][s] = 0;

				}
			}

		}
		#pragma endregion
		#pragma region//********** KERR NONLINEARITY ***********//
		double alphaN, chi3, NLiter;
		double kerrStart, kerrEnd;
		if (kerrOn == "on")
		{
			cout << "Kerr Nonlinearity on" << endl;
			settings.get("KerrMedium/alpha", alphaN);
			settings.get("KerrMedium/chi3", chi3);
			settings.get("KerrMedium/startPosition", kerrStart);
			settings.get("KerrMedium/endPosition", kerrEnd);
			settings.get("KerrMedium/iterations", NLiter);
		}
		#pragma endregion
		#pragma region//********** RAMAN NONLINEARITY **********//
		double ramanTau1, ramanTau2;
		double ramanStart, ramanEnd;
		double alphaRaman = NULL, betaRaman = NULL, kappaRaman = NULL;
		double *Praman = NULL, *PramanOld = NULL, *PramanOld2 = NULL;
		if (ramanOn == "on")
		{
			cout << "Raman Nonlinearity on" << endl;
			settings.get("RamanMedium/tau1", ramanTau1);
			settings.get("RamanMedium/tau2", ramanTau2);
			settings.get("KerrMedium/startPosition", ramanStart);
			settings.get("KerrMedium/endPosition", ramanEnd);
			settings.get("KerrMedium/alpha", alphaN);
			settings.get("KerrMedium/chi3", chi3);
			Praman = new double [N];
			PramanOld = new double [N];
			PramanOld2= new double [N];

			for( k = 0; k < N; k++)
			{
					Praman[k] = 0;
					PramanOld[k] = 0;
					PramanOld2[k] = 0;
			}
		}
		#pragma endregion
		cout<< "" << endl;

	#pragma endregion //%%%%%%%%%%%% END INITIALIZE MEDIA %%%%%%%%%%%%%//
	#pragma region //%%%%%%%% DOMAIN PRE-LOOP COMPUTATIONS %%%%%%%%%//
	
	sigx = 0;
	muy = 1;	

	// Start file and input pertinent data
	parameterValues[1] = N;
	ofstream data(fileName.c_str(), ios::out | ios::binary | ios::trunc);    // Open file to store data
	for(k = 0; k<PARAM; k++)
		data.write((char *) &parameterValues[k], sizeof parameterValues[k]);
	data.close();

	// Domain computations
	gridSize = N;
	dx = spaceSize/gridSize;
	dt = .5*dx/CNOT;
	maxIter = (int) round(timeDom/dt);
	skipSteps = (int) round(maxIter/keepFrames);
	T = 0;

	
	// Coefficents
	m = 0;
	on = 0;
	for (k = 0; k < N; k++)
	{
		double epsTemp = 1;
		c1[k] = (1-(sigx*dt)/(2*EPSNOT*epsTemp))/(1+(sigx*dt)/(2*EPSNOT*epsTemp));
		c2[k] = (dt/(EPSNOT*epsTemp))/(1+(sigx*dt)/(2*EPSNOT*epsTemp));
		d1[k] = dt/(MUNOT*muy);
	}
	// Source
	maxT = width*2;
	timeSwitch = maxT*(2);
	lambda = CNOT/f;
	omega = 2*PI*f;

	// PML computations
	sigMax = maxC*(sigM+1)*.8/(ETANOT*dx);   // Max value of sigma in PML

	m = PML - 1;
	for (k = 0; k < PML; k++)
	{
		xVec = m + .5;
		m = m-1;

		sigVec		= sigMax * pow(xVec + .5, sigM)/(pow(pmlWidth, sigM));
		sigVecStag  = sigMax * pow(xVec, sigM)/(pow(pmlWidth, sigM));

		kappaVec	= 1 + (kappaMax-1)* pow(xVec + .5, kappaM)/(pow(pmlWidth, kappaM));
		kappaVecStag= 1 + (kappaMax-1)* pow(xVec , kappaM)/(pow(pmlWidth, kappaM));

		aVec		= aMax * pow(xVec + .5, aM)/(pow(pmlWidth, aM));
		aVecStag	= aMax * pow(xVec, aM)/(pow(pmlWidth, aM));

		tau			= (kappaVec*EPSNOT)/(kappaVec*aVec + sigVec);
		tauStag		= (kappaVecStag*EPSNOT)/(kappaVecStag*aVecStag + sigVecStag);
		 
		bx[k] = exp(-dt/tau);
		by[k] = exp(-dt/tauStag);

		gx[k] = (sigVec/(kappaVec*(kappaVec*aVec+sigVec)))*(1-bx[k]);
		gy[k] = (sigVecStag/(kappaVecStag*(kappaVecStag*aVecStag+sigVecStag)))*(1-by[k]);

		kx[k] = kappaVec;
		ky[k] = kappaVecStag;
	}


	#pragma endregion //%%%%%% END DOMAIN PRE-LOOP COMPUTATIONS %%%%%%%//
	#pragma region //%%%%%%%%% MEDIA PRE-LOOP COMPUTATIONS %%%%%%%%%//

		#pragma region //************** DIELECTRIC **************//
		if (dielectricOn == "on")
		{
			dielStart = (int) round(dielStart/dx);
			dielEnd   = (int) round(dielEnd/dx);
			for ( k = (int) dielStart-1 ; k < (int) dielEnd; k++ )
			{ 
				epsx[k] = epsRel;
			}
		}
		#pragma endregion
		#pragma region //*************** LORENTZ ****************//
		if (lorentzOn == "on")
		{
			for (k = 0; k<lorentzPoles; k++)
			{
				double delEps   = epsLorentz[k] - epsInfLorentz;
				double omega    = 2*PI*fLorentz[k];
				double omegaSQ  = pow(omega,2.0);
				double denom    = delLorentz[k]*dt + 1;
				alphaLorz[k] = (EPSNOT*delEps*omegaSQ*pow(dt,2))/denom;
				betaLorz[k]  = (2 - omegaSQ*pow(dt,2))/denom;
				gammaLorz[k] = (delLorentz[k]*dt - 1)/denom;
			}
			lorzStart = (int) round(lorzStart/dx);
			lorzEnd   = (int) round(lorzEnd/dx);

			for ( k = (int) lorzStart-1 ; k < (int) lorzEnd; k++ )
			{ 
				epsx[k] = epsInfLorentz;
			}
		}
		#pragma endregion
		#pragma region//**************** DEBYE *****************//
		if (debyeOn == "on")
		{
			for (k = 0; k<debyePoles; k++)
			{
				double delEps   = epsDebye[k] - epsInfDebye;
				double denom    = tauDebye[k]/(2*dt);
				alphaDebye[k] = (EPSNOT*delEps)/denom;
				betaDebye[k]  = -1/denom;
			}
			debyeStart = (int) round(debyeStart/dx);
			debyeEnd   = (int) round(debyeEnd/dx);

			for ( k = (int) debyeStart-1 ; k < (int) debyeEnd; k++ )
			{ 
				epsx[k] = epsInfDebye;
			}
		}
		#pragma endregion
		#pragma region//**************** DRUDE *****************//
		if (drudeOn == "on")
		{
			for (k = 0; k<drudePoles; k++)
			{
				double omega    = 2*PI*fDrude[k];
				double omegaSQ  = pow(omega,2.0);
				double denom    = 1 + ( gammaDrude[k]*dt/2 );
				alphaDrude[k] = ( EPSNOT*pow(dt,2.0)*omegaSQ )/denom;
				betaDrude[k]  = 2/denom;
				kappaDrude[k] = ( ( gammaDrude[k]*dt/2 ) - 1 )/denom;
			}
			drudeStart = (int) round(drudeStart/dx);
			drudeEnd   = (int) round(drudeEnd/dx);
		}
		#pragma endregion
		#pragma region//*********** 2LVL ATOM MODEL ************//
		if (blochOn == "on")
		{
			double cs = (odeNumber)/( 4*sqrt( 2*log(2.0) ) );
			double as = 3*TNatom/(odeNumber);
			double tempSum = 0;
			double delFreq = (blochFreqMax - blochFreqMin)/(odeNumber-1);

			for (k = 0; k<odeNumber; k++)
			{
				double u = k - odeNumber/2;
				Natom[k] = as*exp( - pow(u,2.0) / pow(cs,2.0));
				tempSum = tempSum + Natom[k];
				omegaBloch[k] = 2*PI*(blochFreqMin + k*delFreq);
			}

			blochStart = (int) round(blochStart/dx);
			blochEnd   = (int) round(blochEnd/dx);
		}
		#pragma endregion
		#pragma region//********** RAMAN NONLINEARITY **********//
		if (ramanOn == "on")
		{
			double omegaRam = sqrt( ( pow(ramanTau1,2.0) + pow(ramanTau2,2.0) ) / ( pow(ramanTau1,2.0) * pow(ramanTau2,2.0) ) );
			double deltaRam = 1/ramanTau2;

			alphaRaman = ((1-alphaN)*chi3*pow(omegaRam,2.0)*pow(dt,2.0))/(deltaRam*dt+1);
			betaRaman  = (2-pow(dt,2.0)*pow(omegaRam,2.0))/(deltaRam*dt+1);
			kappaRaman = (deltaRam*dt-1)/(deltaRam*dt+1);

			ramanStart = (int) round(ramanStart/dx);
			ramanEnd   = (int) round(ramanEnd/dx);
		}
		#pragma endregion
		cout << "Beginning Calculations ..." << endl;
		double start = clock();

	#pragma endregion //%%%%%%% END MEDIA PRE-LOOP COMPUTATIONS %%%%%%%//

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	for ( n = 1; n <= maxIter; n++)
	{

		#pragma region //%%%%%%%%%%%%%%%%%%% SOURCE %%%%%%%%%%%%%%%%%%%%//
		T = T+1; // Loop counter

		if (pulseSetting == "modulatedGauss")
		source = exp( -1.0*pow(maxT-(double)T*dt,2.0) / pow(width/5.0, 2.0) ) * pAmp*sin(omega*double(T)*dt);
		if(pulseSetting == "sech")
		{
			double inputT = PI*(dt*T - width*2)/width;
			source = pAmp*sin(omega*double(T)*dt) * 2 / ( exp(inputT) + exp(-(inputT)) );
		}
		if(pulseSetting == "gauss")
			source = pAmp*exp( -1.0*pow(maxT-T*dt,2.0) / pow(width/5.0, 2.0) );
		if(pulseSetting == "CW")
			source = pAmp*sin(omega*double(T)*dt);


		if (blochOn == "on"	)
			Ex[sourcePos] = source+ Ex[sourcePos];
		else
			Dx[sourcePos] = source*EPSNOT + Dx[sourcePos];
		
		# pragma endregion //%%%%%%%%%%%%%%%%% END SOURCE %%%%%%%%%%%%%%%%%%//
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEDIUM UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%%%% LORENTZ MEDIUM UPDATE %%%%%%%%%%%%//
		if (lorentzOn == "on")
		{
			for ( k = (int) lorzStart-1 ; k < (int) lorzEnd; k++ )
			{   
				for (int s = 0; s<lorentzPoles; s++ )
				{
				Plorz[k][s] = alphaLorz[s]*Ex[k] + betaLorz[s]*Plorz[k][s] + gammaLorz[s]*PlorzOld2[k][s];
				PL[k] = PL[k] + Plorz[k][s];
				PlorzOld2[k][s] = PlorzOld[k][s];
				PlorzOld[k][s] = Plorz[k][s];
				}
			}
		}
		
		#pragma endregion //%%%%%%%%%% END LORENTZ MEDIUM UPDATE %%%%%%%%%%//
		#pragma region //%%%%%%%%%%%%% DEBYE MEDIUM UPDATE %%%%%%%%%%%%%//
		if (debyeOn == "on")
		{
			for ( k = (int) debyeStart-1 ; k < (int) debyeEnd; k++ )
			{   
				for (int s = 0; s<debyePoles; s++ )
				{
				Pdebye[k][s] = alphaDebye[s]*Ex[k] + betaDebye[s]*Pdebye[k][s] + PdebyeOld2[k][s];
				PL[k] = PL[k] + Pdebye[k][s];
				PdebyeOld2[k][s] = PdebyeOld[k][s];
				PdebyeOld[k][s] = Pdebye[k][s];
				}
			}
		}
		#pragma endregion //%%%%%%%%%%% END DEBYE MEDIUM UPDATE %%%%%%%%%%%//
		#pragma region //%%%%%%%%%%%%% DRUDE MEDIUM UPDATE %%%%%%%%%%%%%//
		if (drudeOn == "on")
		{
			for ( k = (int) drudeStart-1 ; k < (int) drudeEnd; k++ )
			{   
				for (int s = 0; s<drudePoles; s++ )
				{
				Pdrude[k][s] = alphaDrude[s]*Ex[k] + betaDrude[s]*Pdrude[k][s] + kappaDrude[s]*PdrudeOld2[k][s];
				PL[k] = PL[k] + Pdrude[k][s];
				PdrudeOld2[k][s] = PdrudeOld[k][s];
				PdrudeOld[k][s] = Pdrude[k][s];
				}
			}
		}
		#pragma endregion //%%%%%%%%%%% END DRUDE MEDIUM UPDATE %%%%%%%%%%%//
		#pragma region //%%%%%%%%%%% 2LVL ATOM MEDIUM UPDATE %%%%%%%%%%%//
		if (blochOn == "on")
		{
			double tHalf = dt*(T-0.5);
			double t = dt*(T);

			blochCm = ( (2*gammaBloch)/DIRAC) * exp( -(tHalf)*(1/blochTau2-1/blochTau1) );
			blochCp = ( (2*gammaBloch)/DIRAC) * exp( -(t)*(1/blochTau1-1/blochTau2) );
			blochD = ( (2*gammaBloch*rho30)/DIRAC ) * exp (t/blochTau2);

			for ( int s = 0; s < odeNumber; s++ )
			{
				blochA[s] = (  (Natom[s]*gammaBloch)/(EPSNOT*blochTau2)  ) * exp(-tHalf/blochTau2);
				blochB[s] = ( (Natom[s]*gammaBloch*omegaBloch[s])/EPSNOT ) * exp(-tHalf/blochTau2);
				double denom = ( 1 + (pow(omegaBloch[s],2.0)*pow(dt,2.0))/4 );
				double cof1 = dt*(blochCm/2);
				double coef1 = 1 - ( pow(omegaBloch[s],2.0)*pow(dt,2.0) )/4 ;
				double coef2 = omegaBloch[s]*dt;
				double coef3 = omegaBloch[s]*pow(dt,2.0)/2;
				double coef4 = (omegaBloch[s]*dt/2);

				for (k = (int) blochStart-1; k< blochEnd; k++ )
				{

					u3[k][s] = u3[k][s] - (dt*blochCm/2)*(Ex[k] + ExStore[k])*u2[k][s];
					
					double Q;
					Q = blochCp*u3[k][s]*Ex[k]  +  blochD*Ex[k];
					u1[k][s] = ( coef1*u1[k][s] + coef2* u2[k][s]  + coef3*Q ) / denom;
					u2[k][s] = u2[k][s] - coef4*(u1[k][s] + u1store[k][s]) + dt*Q;

					PL[k] = PL[k] + dt*(-blochA[s]*u1[k][s] + blochB[s]*u2[k][s]);
					
					u1store[k][s] = u1[k][s];

				}
				
			}

			for (k = (int) blochStart-1; k< blochEnd; k++ )
				ExStore[k] = Ex[k];
		}


		#pragma endregion //%%%%%%%%% END 2LVL ATOM MEDIUM UPDATE %%%%%%%%%//
		#pragma region //%%%%%%%%%%%%% RAMAN MEDIUM UPDATE %%%%%%%%%%%%%//
		if (ramanOn == "on")
		{
			for ( k = (int) ramanStart-1 ; k < (int) ramanEnd; k++ )
			{   
			
				Praman[k] = alphaRaman*pow(Ex[k],2.0) + betaRaman*Praman[k] + kappaRaman*PramanOld2[k];
				PramanOld2[k] = PramanOld[k];
				PramanOld[k] = Praman[k];
			}
		}
		#pragma endregion //%%%%%%%%%%% END RAMAN MEDIUM UPDATE %%%%%%%%%%%//
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MEDIUM UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIELD UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%%%%%%% H-FIELD UPDATE %%%%%%%%%%%%%%%%//
			for ( k = PML; k < N-1-PML; k++ )
			{
				curlE = -(Ex[k+1] - Ex[k])/dx;
				Hy[k] = Hy[k] + d1[k]*curlE;
			}
		#pragma endregion //%%%%%%%%%%%%% END H-FIELD UPDATE %%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%%%% LINEAR E-FIELD UPDATE %%%%%%%%%%%%//
		if (ramanOn != "on"	&& kerrOn != "on"  && blochOn != "on" )
		{
			for ( k = PML ; k < N-PML; k++ )
			{   
				curlH = -(Hy[k] - Hy[k-1])/dx;
				Dx[k] = Dx[k] + dt*curlH;
				Ex[k] = (Dx[k] - PL[k])/(EPSNOT*epsx[k]);
				PL[k] = 0;
			}
		}
		#pragma endregion //%%%%%%%%%% END LINEAR E-FIELD UPDATE %%%%%%%%%%//
		#pragma region //%%%%%%%%%% NONLINEAR E-FIELD UPDATE %%%%%%%%%%%//
		if (ramanOn == "on"	|| kerrOn == "on" )
		{
			for ( k = PML ; k < N-PML; k++ )
			{   
				curlH = -(Hy[k] - Hy[k-1])/dx;
				Dx[k] = Dx[k] + dt*curlH;
				Exlin[k] = (Dx[k] - PL[k]);
				PL[k] = 0;
				ExNL[k] = Ex[k];
			}

			double counter = 0;
			
			while (counter < NLiter)
			{
				for ( k = PML ; k < N-PML; k++ )
				{   
					Ex[k] = Exlin[k]/( EPSNOT*(epsx[k] + Praman[k] + alphaN*chi3*pow(ExNL[k],2.0) ) );
					ExNL[k] = Ex[k];
				}
				counter = counter + 1;
			}
			
		}
		if (blochOn == "on"	)
		{
			for ( k = PML ; k < N-PML; k++ )
			{  
				curlH = -(Hy[k] - Hy[k-1])/dx;
				Ex[k] = c1[k]*Ex[k] + c2[k]*curlH + PL[k];
				PL[k] = 0;
			}
		}
		#pragma endregion //%%%%%%%% NONEND LINEAR E-FIELD UPDATE %%%%%%%%%//
		#pragma region //%%%%%%%%%%%%%% PML REGION UPDATE %%%%%%%%%%%%%%//
		//Hy Update
		for ( k = 0; k < PML; k++ )
		{
			curlE = -(Ex[k+1] - Ex[k])/dx;
			QEx[k] = by[k]*QEx[k] - gy[k]*curlE;
			Hy[k] = Hy[k] + d1[k]*ky[k]*( curlE + QEx[k] );
		}
		m = PML-1;
		for ( k = N-PML-1; k < N-1; k++ )
		{
			curlE = -(Ex[k+1] - Ex[k])/dx;
			QEx[m+PML] = by[m]*QEx[m+PML] - gy[m]*curlE;
			Hy[k] = Hy[k] + d1[k]*ky[m]*( curlE + QEx[m+PML] );
			m--;
		}
		//Ex Update
		for ( k = 1; k < PML; k++ )
		{
			curlH = -(Hy[k] - Hy[k-1])/dx;
			QHy[k] = bx[k]*QHy[k] - gx[k]*curlH;
			Ex[k] = c1[k]*Ex[k] + c2[k]*kx[k]*(curlH + QHy[k]);
		} 
		m = PML-1;
		for ( k = N-PML; k < N-1; k++ )
		{
			curlH = -(Hy[k] - Hy[k-1])/dx;
			QHy[m+PML] = bx[m]*QHy[m+PML] - gx[m]*curlH;
			Ex[k] = c1[k]*Ex[k] + c2[k]*kx[m]*(curlH + QHy[m+PML]);
			m--;
		}

		// Ex Boundary Update //
		Ex[1] = 0;
		Ex[N-1] = 0;
		#pragma endregion //%%%%%%%%%%%% END PML REGION UPDATE %%%%%%%%%%%%//
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END FIELD UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

		#pragma region //%%%%%%%%%%% FREQUENCY DOMAIN UPDATE %%%%%%%%%%%//
		for ( k = 0; k<FNUM; k++)
		{
			testFreq = freqMin + k*(freqMax - freqMin)/FNUM;

			if (n*dt < timeSwitch)
			{
			realEinc[k] = realEinc[k] + Ex[detInit]*cos( 2.0*PI * (double)T * dt*testFreq);
			imagEinc[k] = imagEinc[k] + Ex[detInit]*sin( 2.0*PI * (double)T * dt*testFreq);
			}
			else
			{
			realEref[k] = realEref[k] + Ex[detInit]*cos( 2.0*PI * (double)T * dt*testFreq);
			imagEref[k] = imagEref[k] + Ex[detInit]*sin( 2.0*PI * (double)T * dt*testFreq);
			}

			realEfin[k] = realEfin[k] + Ex[detFin ]*cos( 2.0*PI * (double)T * dt*testFreq);
			imagEfin[k] = imagEfin[k] + Ex[detFin ]*sin( 2.0*PI * (double)T * dt*testFreq);
		}
		#pragma endregion //%%%%%%%%% END FREQUENCY DOMAIN UPDATE %%%%%%%%%//	            
		#pragma region //%%%%%%%%%%%%%% WRITE VIDEO DATA %%%%%%%%%%%%%%%//
		if	(T%skipSteps == 0)
		{
			ofstream data(fileName.c_str(), ios::out | ios::binary | ios::app);    // Open file to store data
		    //os.write(reinterpret_cast<const char*> (s), sizeof(double) * ns);
			for (k=0; k<N; k++)
				data.write( (char*) &Ex[k], sizeof Ex[k]);
			data.close();

			/*
			ofstream rho3data("rho3data.bin", ios::out | ios::binary | ios::app);    // Open file to store data
			for (int s = 0; s< odeNumber; s++)
			{
				for (k=0; k<N; k++)
				{
					double rho3[1];
					//rho3[0] = rho30 + exp(-dt*T/blochTau1)*u3[k][s];
					rho3[0] = exp(-dt*T/blochTau2)*u1[k][s];
					rho3data.write( (char*) &rho3[0], sizeof rho3[0]);
				}
			}
			rho3data.close();
			*/
		}
		#pragma endregion //%%%%%%%%%%%% END WRITE VIDEO DATA %%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%%% PRINT SIMULATION INFO %%%%%%%%%%%%%//
		if (T == 1 && maxIter > 999)
		cout<< "Calculating step " << (double) T << " to " << T+998 << " out of " << maxIter << endl;
		if (T == 1 && maxIter < 999)
		cout<< "Calculating step " << (double) T << " to " << maxIter << " out of " << maxIter << endl;

		if ( T%1000 == 0 && T < (maxIter - 1000) )
		// Display info
		{
			cout<< "Calculating step " << (double) T << " to " << T+999 << " out of " << maxIter << endl;
			double end = clock();
			double time = end-start;
			cout << "                Current simulation time = " <<  time/CLOCKS_PER_SEC <<  " sec" <<endl;
			cout << "" <<endl;
		}

		if ( T%1000 == 0 && T < (maxIter)  && T > (maxIter - 1000) )
		// Display info
		{
			cout<< "Calculating step " << (double) T << " to " << maxIter << " out of " << maxIter << endl;
			double end = clock();
			double time = end-start;
			cout << "                Current simulation time = " <<  time/CLOCKS_PER_SEC <<  " sec" <<endl;
			cout << "" <<endl;
		}
		#pragma endregion //%%%%%%%%% END PRINT SIMULATION INFO %%%%%%%%%%%//

	}//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//					

	#pragma region //%%%%%%% PROCESS FREQUENCY DOMAIN OUTPUT %%%%%%%//
	for ( k = 0; k<FNUM; k++)
	{ 
		ampEinc[k] = pow(sqrt( pow(realEinc[k],2.0) + pow(imagEinc[k],2.0) ),2.0);
		ampEref[k] = pow(sqrt( pow(realEref[k],2.0) + pow(imagEref[k],2.0) ),2.0);
		ampEfin[k] = pow(sqrt( pow(realEfin[k],2.0) + pow(imagEfin[k],2.0) ),2.0);
		
		phaseEinc[k] = atan2(imagEinc[k], realEinc[k]);
		phaseEref[k] = atan2(imagEref[k], realEref[k]);
		phaseEfin[k] = atan2(imagEfin[k], realEfin[k]);

		transmission[k] = ampEfin[k]/ampEinc[k];
		reflection[k] = ampEref[k]/ampEinc[k];	
	}
	#pragma endregion //%%%%% END PROCESS FREQUENCY DOMAIN OUTPUT %%%%%//	
	#pragma region //%%%%%%%% WRITE FREQUENCY DOMAIN OUTPUT %%%%%%%%//
	// Write binary data to file
	cout << "Writing Data..." << endl;

	double fnum[3];
	fnum[0] = FNUM;
	fnum[1] = freqMin;
	fnum[2] = freqMax;


	ofstream data2("FDTDdataFREQ.bin", ios::out | ios::binary );// Open file to store data
	data2.write((char *) &fnum, sizeof fnum);
	data2.close();

	ofstream data3("FDTDdataFREQ.bin", ios::out | ios::binary | ios::app);// Open file to store data
	for (k=0; k<FNUM; k++)
		data3.write((char *) &ampEinc[k], sizeof ampEinc[k]);
	for (k=0; k<FNUM; k++)
	    data3.write((char *) &ampEref[k], sizeof ampEref[k]);
	for (k=0; k<FNUM; k++)
		data3.write((char *) &ampEfin[k], sizeof ampEfin[k]);
	for (k=0; k<FNUM; k++)
		data3.write((char *) &reflection[k], sizeof reflection[k]);
	for (k=0; k<FNUM; k++)
		data3.write((char *) &transmission[k], sizeof transmission[k]);
	data3.close();

	#pragma endregion //%%%%%% END WRITE FREQUENCY DOMAIN OUTPUT %%%%%%//
	#pragma region //%%%%%%%%%%%%%%% DISPLAY OUTPUT %%%%%%%%%%%%%%%%//
	// Display run time
	double end = clock();
	double time = end-start;
	cout << "Simulation complete!" << endl;
	cout << "Total Simulation Time = " <<  time/CLOCKS_PER_SEC <<  " sec" <<endl;

	// Display simulation parameters
	cout << "" <<endl;
	cout << "SIMULATION PARAMETERS:" << endl;
	cout << "x-Domain size = " << dx*N*1e6        << "microns           " << "dx = " << dx <<endl;
	cout << "t-Domain size = " << dt*maxIter*1e15 << "fs           " << "dt = " << dt <<endl;
	cout << "Freq domain interval = " << freqMin/1e12 << " - " << freqMax/1e12 << "Thz, with " << FNUM << " steps" <<endl;
	cout << "" <<endl;
	cout << "Press enter to exit cmd window" << endl;

	cin.get();
	return(0);

	#pragma endregion //%%%%%%%%%%%%% END DISPLAY OUTPUT %%%%%%%%%%%%%%//

}//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


