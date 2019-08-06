
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Info:                                                           %
% Simulates wave propagation in axisymmetric 2D using the FDTD method     %
% which yields second order accuracy on the traditional YEE staggered     %
% grid.																	  %	
%																		  %
% 1. Computational domain:												  %
%																		  %
%      ^Z																  %
%      |																  %
%      |																  %
%      |																  %
%	   |_______________________________________________________           %
%	   |				PML REGION (1) 				      |    |          %
%	   |__________________________________________________|    |          %
%	   |											      |    |          %
%	  A|-----------Transmission Caculated Here------------|    |          %
%	  X|											      |    |          %
%	  I|											      |    |          %
%	  S|---------------------------------			      |    |          %
%	   |								|			      |    |          %
%	  O|								|			      |    |          %
%	  F|								|			      |    |          %
%	   |								|			      |    |          %
%	  S|		Scattering				|			      | P  |          %
%	  Y|		  Object				|			      | M  |          %
%	  M|								|			      | L  |          %
%	  M|								|			      |(2) |          %
%	  E|								|			      |    |          %
%	  T|								|			      |    |          %
%	  R|---------------------------------   		      |    |          %
%	  Y|			^								      |    |           %
%	   |			| Principle						      |    |          %
%	   |			| Propagation					      |    |          %
%	   |			|								      |    |          %
%	   |			|								      |    |          %
%	   |			|								      |    |          %
%	   |-------Source Condition-------    			      |    |          %
%	   |											      |    |          %
%	   |------------Reflection Caculated Here-------------|    |          %
%	   |__________________________________________________|    |          %
%	   |				PML REGION (3) 				      |    |          %
%	   |__________________________________________________|____|----->R   %
%																		  %
%																		  %
% 2. Input file provides settings which are imported into the program     %
% 3. Utilizes full vector representations for E(r,z), H(r,z), i.e.        %
%    componenets Er, Ez, Ephi, Hr, Hz, Hphi.                              %
% 4. CPML utilized to truncate spatial domain                             %
% 5. Source condition can be varied in the settings (e.g. SECH, Gaussian) %
% 6. Various Media can be implemented in the scattering object            %
%    I. Linear:                                                           %
%    	i.   Dielectric                                                   %
%	    ii.  Lorentz Poles                                                %
%		iii. Drude Poles                                                  %
%		iv.  Debye Poles                                                  %
%    II. Nonlinear                                                        %
%       i.   Kerr                                                         %
%       ii.  Raman                                                        %
%       iii. Multiple 2LVL Atom Quantum Model                             %
% 7. Frequency domain is recorded at intial and final points, and         %
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
	double count = 0;
	#pragma region //%%%%%%%%%%%%%%% IMPORT SETTINGS %%%%%%%%%%%%%%%//
	int detInit1, detFin1, gridR, gridZ, maxIter, sourcePos, detInit, detFin, N, PML, FNUM, PARAM, medSect, M;
	string fileName, pulseSetting;
	double rRange, zRange, timeDom, f, pAmp, width, freqMin, freqMax, keepFrames, sigx, refractiveIndex;
	double pmlWidth, sigM, kappaMax, kappaM, aMax, aM, maxC;
	double curlEphi, curlEr, curlEz, curlHphi, curlHz, curlHr;
	double zNot, focusWidth;
	
	Config settings("FDTDsettings.txt");
	//settings.printSettings(std::cout);
	std::cout << "" << endl;
	std::cout << "" << endl;

    settings.get("Domain/radialGrid", gridR);
	N = gridR;
	settings.get("Domain/zGrid", gridZ);
	settings.get("Domain/radialRange", rRange);
	settings.get("Domain/zRange", zRange);
	settings.get("Domain/timeDomain", timeDom);
	settings.get("Domain/harmonicIndex", M);


	double dt, dr, dz;
	double dx = 0;
	dr = rRange/gridR;
	dz = zRange/gridZ;
	if (M == 0)
		dt = .5/( (1)*(CNOT)*sqrt( 1/pow(dr,2.0)+1/pow(dz,2.0) ) ); 
	else
		dt = .5/( (M)*(CNOT)*sqrt( 1/pow(dr,2.0)+1/pow(dz,2.0) ) ); 
	maxIter = (int) round(timeDom/dt);

	std::cout << "Domain (R x Z) " << rRange*1e6       << " microns x " << zRange*1e6       << " microns   ---    " << 
	"Time Domain = " << timeDom*1e15 << "fs " <<endl;
	cout << "GRID (R x Z x T) : " << gridR << " x " << gridZ << " x " << maxIter << endl;
	cout << " " <<endl;
	settings.get("InputPulse/frequency", f);
	settings.get("InputPulse/amplitude", pAmp);
	settings.get("InputPulse/temporalWidth", width);
	settings.get("InputPulse/position", sourcePos);
	settings.get("InputPulse/pulse", pulseSetting);
	settings.get("InputPulse/zNot", zNot);
	settings.get("InputPulse/focusWidth", focusWidth);
	
	if (pulseSetting == "modulatedGauss")
	{
		std::cout <<  "SOURCE-- Z: modulated gaussian, R: plane wave" << endl;
		std::cout <<  "Frequency: " << f/1e12 << " THz                Pulse Temporal Width: " << width*1e15 << " fs" << endl;
	}
	if(pulseSetting == "sech")
	{
		std::cout <<  "SOURCE-- Z: modulated sech, R: plane wave" << endl;
		std::cout <<  "Frequency: " << f/1e12 << " THz                Pulse Temporal Width: " << width*1e15 << " fs" << endl;
	}
	if(pulseSetting == "gauss")
	{
		std::cout <<  "SOURCE-- Z: gaussian, R: plane wave" << endl;
		std::cout <<  "Frequency: " << "N/A" << " THz                Pulse Temporal Width: " << width*1e15 << " fs" << endl;
	}
	if(pulseSetting == "CW")
	{
		std::cout <<  "SOURCE-- Z:CW, R: plane wave " << endl;
		std::cout <<  "Frequency: " << f/1e12  << " THz                Pulse Temporal Width: " << "N/A" << " fs" << endl;
	}
	if( pulseSetting == "gaussBeamCW" )
	{
		std::cout <<  "SOURCE-- Z,R: CW gaussian beam" << endl;
		std::cout << "focal width: " <<  focusWidth*1e6 << " microns          focal distance: " << zNot*1e6  << " microns" << endl;
		std::cout <<  "Frequency: " << f/1e12  << " THz,                Pulse Temporal Width: " << "N/A" << " fs" << endl;
	}
	if( pulseSetting == "gaussBeamPulse" )
	{
		std::cout <<  "SOURCE-- Z,R: Pulsed gaussian beam, focal width:" << endl;
		std::cout << "focal width: " <<  focusWidth*1e6 << " microns          focal distance: " << zNot*1e6  << " microns" << endl;
		std::cout <<  "Frequency: " << f/1e12 << " THz                Pulse Temporal Width: " <<  width*1e15 << " fs" << endl;
	}

	std::cout << " " << endl;

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

	string debyeOn, drudeOn, lorentzOn, nonlinOn, blochOn, dielectricOn;
	settings.get("LorentzMedium/LorentzOn", lorentzOn);
	settings.get("DrudeMedium/DrudeOn", drudeOn);
	settings.get("DebyeMedium/DebyeOn", debyeOn);
	settings.get("2LVLAtom/2LVLAtomOn", blochOn);
	settings.get("NonlinearMedium/NonlinearOn", nonlinOn);
	settings.get("DielectricMedium/DielectricOn", dielectricOn);

	#pragma endregion //%%%%%%%%%%%%% END IMPORT SETTINGS %%%%%%%%%%%%%//	
	#pragma region //%%%%%%%%%%%%%% INITIALIZE DOMAIN %%%%%%%%%%%%%%//
	// Initialize variables and constants
	double sigMax;									                // Domain Variables
	double muy;											            // Coef. constants
	int n, k, m, T, s, on;   							            // Counters
	double source, maxT, lambda, omega, timeSwitch;                 // Source constants
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
	double* bn = new double [PML];
	double* bs = new double [PML];
	double* gn = new double [PML];
	double* gs = new double [PML];
	double* kn = new double [PML];
	double* ks = new double [PML];
	double* bn2 = new double [PML];
	double* bs2 = new double [PML];
	double* gn2 = new double [PML];
	double* gs2 = new double [PML];
	double* kn2 = new double [PML];
	double* ks2 = new double [PML];

	double** QEzphi = new double* [gridR];
	double** QEzr   = new double* [gridR];
	double** QErz   = new double* [gridR];
	double** QErphi = new double* [gridR];
	double** QHzphi = new double* [gridR];
	double** QHrphi = new double* [gridR];
	double** QHzr   = new double* [gridR];
	double** QHrz   = new double* [gridR];

	for (k = 0; k< gridR; k++)
	{
		QEzphi[k] = new double [gridZ-1];
		QEzr[k]   = new double [gridZ-1];
		QErz[k]   = new double [gridZ-1];
		QErphi[k] = new double [gridZ];
		QHzphi[k] = new double [gridZ];
		QHrphi[k] = new double [gridZ-1];
		QHzr[k]   = new double [gridZ];
		QHrz [k]  = new double [gridZ];
	}

	for (k = 0; k< gridR; k++)
	{
		for (s = 0; s<gridZ; s++)
		{
			
			QErphi[k][s] = 0;
			QHzphi[k][s] = 0;
			QHzr[k][s]   = 0;
			QHrz [k][s]  = 0;
		}

		for (s = 0; s<gridZ-1; s++)
		{
		QEzphi[k][s] = 0;
		QEzr[k][s]   = 0;
		QErz[k][s]   = 0;
		QHrphi[k][s] = 0;
		}
	}

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

	// SOLUTION ARRAYS
	double** Hr   = new double* [gridR];
	double** Hz   = new double* [gridR];
	double** Hphi = new double* [gridR];

	double** Er   = new double* [gridR];
	double** Ez   = new double* [gridR];
	double** Ephi = new double* [gridR];

	double** ErNL   = new double* [gridR];
	double** EzNL   = new double* [gridR];
	double** EphiNL = new double* [gridR];
	
	double** Erlin   = new double* [gridR];
	double** Ezlin   = new double* [gridR];
	double** Ephilin = new double* [gridR];

	double** Dr   = new double* [gridR];
	double** Dz   = new double* [gridR];
	double** Dphi = new double* [gridR];

	double** PLr   = new double* [gridR];
	double** PLz   = new double* [gridR];
	double** PLphi = new double* [gridR];

	double** Jr   = new double* [gridR];
	double** Jz   = new double* [gridR];
	double** Jphi = new double* [gridR];

	double** epsr   = new double* [gridR];
	double** epsz   = new double* [gridR];
	double** epsphi = new double* [gridR];

	double** sigr   = new double* [gridR];
	double** sigz   = new double* [gridR];
	double** sigphi = new double* [gridR];

	double** mur   = new double* [gridR];
	double** muz   = new double* [gridR];
	double** muphi = new double* [gridR];

	for( k = 0; k < gridR; k++)
	{
		Hz[k]     = new double [gridZ];
		muz[k]    = new double [gridZ];
		Er[k]     = new double [gridZ];
		Ephi[k]   = new double [gridZ];
		ErNL[k]     = new double [gridZ];
		EphiNL[k]   = new double [gridZ];
		Erlin[k]     = new double [gridZ];
		Ephilin[k]   = new double [gridZ];
		Dr[k]     = new double [gridZ];
		Dphi[k]   = new double [gridZ];
		PLr[k]    = new double [gridZ];
		PLphi[k]  = new double [gridZ];
		Jr[k]     = new double [gridZ];
		Jphi[k]   = new double [gridZ];
		epsr[k]   = new double [gridZ];
		epsphi[k] = new double [gridZ];
		sigr[k]   = new double [gridZ];
		sigphi[k] = new double [gridZ];

	}
	for( k = 0; k < gridR; k++)
	{
		Hr[k]    = new double [gridZ-1];
		Hphi[k]  = new double [gridZ-1];
		mur[k]   = new double [gridZ-1];
		muphi[k] = new double [gridZ-1];
		Ez[k]    = new double [gridZ-1];
		EzNL[k]  = new double [gridZ-1];
		Ezlin[k]  = new double [gridZ-1];
		Dz[k]    = new double [gridZ-1];
		PLz[k]   = new double [gridZ-1];
		Jz[k]    = new double [gridZ-1];
		epsz[k]  = new double [gridZ-1];
		sigz[k]  = new double [gridZ-1];
	}


	for( k = 0; k < gridR; k++)
	{
		for( s = 0; s < gridZ; s++)
		{
		Hz[k][s]     = 0;
		muz[k][s]    = 1;
		Er[k][s]     = 0;
		Ephi[k][s]   = 0;
		ErNL[k][s]     = 0;
		EphiNL[k][s]   = 0;
		Erlin[k][s]     = 0;
		Ephilin[k][s]   = 0;
		Dr[k][s]     = 0;
		Dphi[k][s]   = 0;
		PLr[k][s]    = 0;
		PLphi[k][s]  = 0;
		Jr[k][s]     = 0;
		Jphi[k][s]   = 0;
		epsr[k][s]   = 1;
		epsphi[k][s] = 1;
		sigr[k][s]   = 0;
		sigphi[k][s] = 0;
		}

		for( s = 0; s < gridZ-1; s++)
		{
			Hr[k][s]    = 0;
			Hphi[k][s]  = 0;
			mur[k][s]   = 1;
			muphi[k][s] = 1;
			Ez[k][s]    = 0;
			EzNL[k][s]    = 0;
			Ezlin[k][s]    = 0;
			Dz[k][s]    = 0;
			PLz[k][s]   = 0;
			Jz[k][s]    = 0;
			epsz[k][s]  = 1;
			sigz[k][s]  = 0;
		}
	}

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
		double dielStartR, dielEndR, dielStartZ, dielEndZ;
		double epsRel;

		if (dielectricOn == "on")
		{
			settings.get("DielectricMedium/radialStart", dielStartR);
			settings.get("DielectricMedium/radialEnd", dielEndR);
			settings.get("DielectricMedium/zStart", dielStartZ);
			settings.get("DielectricMedium/zEnd", dielEndZ);
			settings.get("DielectricMedium/epsRel", epsRel);
		}
		#pragma endregion
		#pragma region //*************** LORENTZ ****************//
		double lorzStartR, lorzEndR, lorzStartZ, lorzEndZ;
		double epsInfLorentz;
		double *epsLorentz = NULL, *fLorentz = NULL, *delLorentz = NULL, *alphaLorz = NULL, *betaLorz = NULL, *gammaLorz = NULL;
		double ***PlorzR = NULL, ***PlorzROld = NULL, ***PlorzROld2 = NULL;
		double ***PlorzZ = NULL, ***PlorzZOld = NULL, ***PlorzZOld2 = NULL;
		double ***PlorzPHI = NULL, ***PlorzPHIOld = NULL, ***PlorzPHIOld2 = NULL;

		if (lorentzOn == "on")
		{
			cout << "Lorentz Medium on" << endl;


			settings.get("LorentzMedium/radialStart", lorzStartR);
			settings.get("LorentzMedium/radialEnd", lorzEndR);
			settings.get("LorentzMedium/zStart", lorzStartZ);
			settings.get("LorentzMedium/zEnd", lorzEndZ);

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
			PlorzR = new double** [lorentzPoles];
			PlorzROld = new double** [lorentzPoles];
			PlorzROld2= new double** [lorentzPoles];

			PlorzZ = new double** [lorentzPoles];
			PlorzZOld = new double** [lorentzPoles];
			PlorzZOld2= new double** [lorentzPoles];

			PlorzPHI = new double** [lorentzPoles];
			PlorzPHIOld = new double** [lorentzPoles];
			PlorzPHIOld2= new double** [lorentzPoles];

			for (k = 0; k< lorentzPoles; k++)
			{
				PlorzR[k] = new double* [gridR];
				PlorzROld[k] = new double* [gridR];
				PlorzROld2[k] = new double* [gridR];

				PlorzZ[k] = new double* [gridR];
				PlorzZOld[k] = new double* [gridR];
				PlorzZOld2[k] = new double* [gridR];

				PlorzPHI[k] = new double* [gridR];
				PlorzPHIOld[k] = new double* [gridR];
				PlorzPHIOld2[k] = new double* [gridR];
			}


			for( k = 0; k < lorentzPoles; k++)
			{
				for (int s = 0; s<gridR; s++)
				{
					PlorzR[k][s]     = new double [gridZ];
					PlorzROld[k][s]  = new double [gridZ];
					PlorzROld2[k][s] = new double [gridZ];

					PlorzPHI[k][s]     = new double [gridZ];
					PlorzPHIOld[k][s]  = new double [gridZ];
					PlorzPHIOld2[k][s] = new double [gridZ];

					PlorzZ[k][s]     = new double [gridZ];
					PlorzZOld[k][s]  = new double [gridZ];
					PlorzZOld2[k][s] = new double [gridZ];

				}
			}
			for( k = 0; k < lorentzPoles; k++)
			{
				for (int s = 0; s<gridR; s++)
				{
					for (m = 0; m < gridZ; m++)
					{
						PlorzR[k][s][m]     = 0;
						PlorzROld[k][s][m]  = 0;
						PlorzROld2[k][s][m] = 0;

						PlorzZ[k][s][m]     = 0;
						PlorzZOld[k][s][m]  = 0;
						PlorzZOld2[k][s][m] = 0;

						PlorzPHI[k][s][m]     = 0;
						PlorzPHIOld[k][s][m]  = 0;
						PlorzPHIOld2[k][s][m] = 0;
					}
				}
			}
		}
		#pragma endregion
		#pragma region//**************** DEBYE *****************//
		double debyeStartR, debyeEndR, debyeStartZ, debyeEndZ;
		double epsInfDebye;
		double *alphaDebye = NULL, *betaDebye = NULL, *tauDebye = NULL, *epsDebye = NULL;
		double ***PdebyeR = NULL, ***PdebyeROld = NULL, ***PdebyeROld2 = NULL;
		double ***PdebyeZ = NULL, ***PdebyeZOld = NULL, ***PdebyeZOld2 = NULL;
		double ***PdebyePHI = NULL, ***PdebyePHIOld = NULL, ***PdebyePHIOld2 = NULL;
		if (debyeOn == "on")
		{
			cout << "Debye Medium on" << endl;
			

			settings.get("DebyeMedium/radialStart", debyeStartR);
			settings.get("DebyeMedium/radialEnd", debyeEndR);
			settings.get("DebyeMedium/zStart", debyeStartZ);
			settings.get("DebyeMedium/zEnd", debyeEndZ);
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

			PdebyeR = new double** [debyePoles];
			PdebyeROld = new double** [debyePoles];
			PdebyeROld2= new double** [debyePoles];

			PdebyeZ = new double** [debyePoles];
			PdebyeZOld = new double** [debyePoles];
			PdebyeZOld2= new double** [debyePoles];

			PdebyePHI = new double** [debyePoles];
			PdebyePHIOld = new double** [debyePoles];
			PdebyePHIOld2= new double** [debyePoles];

			for (k = 0; k< debyePoles; k++)
			{
				PdebyeR[k] = new double* [gridR];
				PdebyeROld[k] = new double* [gridR];
				PdebyeROld2[k] = new double* [gridR];

				PdebyeZ[k] = new double* [gridR];
				PdebyeZOld[k] = new double* [gridR];
				PdebyeZOld2[k] = new double* [gridR];

				PdebyePHI[k] = new double* [gridR];
				PdebyePHIOld[k] = new double* [gridR];
				PdebyePHIOld2[k] = new double* [gridR];
			}


			for( k = 0; k < debyePoles; k++)
			{
				for (int s = 0; s<gridR; s++)
				{
					PdebyeR[k][s]     = new double [gridZ];
					PdebyeROld[k][s]  = new double [gridZ];
					PdebyeROld2[k][s] = new double [gridZ];

					PdebyePHI[k][s]     = new double [gridZ];
					PdebyePHIOld[k][s]  = new double [gridZ];
					PdebyePHIOld2[k][s] = new double [gridZ];

					PdebyeZ[k][s]     = new double [gridZ];
					PdebyeZOld[k][s]  = new double [gridZ];
					PdebyeZOld2[k][s] = new double [gridZ];

				}
			}
			for( k = 0; k < debyePoles; k++)
			{
				for (int s = 0; s<gridR; s++)
				{
					for (m = 0; m < gridZ; m++)
					{
						PdebyeR[k][s][m]     = 0;
						PdebyeROld[k][s][m]  = 0;
						PdebyeROld2[k][s][m] = 0;

						PdebyeZ[k][s][m]     = 0;
						PdebyeZOld[k][s][m]  = 0;
						PdebyeZOld2[k][s][m] = 0;

						PdebyePHI[k][s][m]     = 0;
						PdebyePHIOld[k][s][m]  = 0;
						PdebyePHIOld2[k][s][m] = 0;
					}
				}
			}



		}
		#pragma endregion
		#pragma region//**************** DRUDE *****************//
		double drudeStartR, drudeEndR, drudeStartZ, drudeEndZ;
		double *fDrude = NULL, *alphaDrude = NULL, *betaDrude = NULL, *kappaDrude = NULL, *gammaDrude = NULL;
		double ***PdrudeR = NULL, ***PdrudeROld = NULL, ***PdrudeROld2 = NULL;
		double ***PdrudeZ = NULL, ***PdrudeZOld = NULL, ***PdrudeZOld2 = NULL;
		double ***PdrudePHI = NULL, ***PdrudePHIOld = NULL, ***PdrudePHIOld2 = NULL;
	
		if (drudeOn == "on")
		{
			cout << "Drude Medium on" << endl;
			settings.get("DrudeMedium/radialStart", drudeStartR);
			settings.get("DrudeMedium/radialEnd", drudeEndR);
			settings.get("DrudeMedium/zStart", drudeStartZ);
			settings.get("DrudeMedium/zEnd", drudeEndZ);
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


			PdrudeR = new double** [drudePoles];
			PdrudeROld = new double** [drudePoles];
			PdrudeROld2= new double** [drudePoles];

			PdrudeZ = new double** [drudePoles];
			PdrudeZOld = new double** [drudePoles];
			PdrudeZOld2= new double** [drudePoles];

			PdrudePHI = new double** [drudePoles];
			PdrudePHIOld = new double** [drudePoles];
			PdrudePHIOld2= new double** [drudePoles];

			for (k = 0; k< drudePoles; k++)
			{
				PdrudeR[k] = new double* [gridR];
				PdrudeROld[k] = new double* [gridR];
				PdrudeROld2[k] = new double* [gridR];

				PdrudeZ[k] = new double* [gridR];
				PdrudeZOld[k] = new double* [gridR];
				PdrudeZOld2[k] = new double* [gridR];

				PdrudePHI[k] = new double* [gridR];
				PdrudePHIOld[k] = new double* [gridR];
				PdrudePHIOld2[k] = new double* [gridR];
			}


			for( k = 0; k < drudePoles; k++)
			{
				for (int s = 0; s<gridR; s++)
				{
					PdrudeR[k][s]     = new double [gridZ];
					PdrudeROld[k][s]  = new double [gridZ];
					PdrudeROld2[k][s] = new double [gridZ];

					PdrudePHI[k][s]     = new double [gridZ];
					PdrudePHIOld[k][s]  = new double [gridZ];
					PdrudePHIOld2[k][s] = new double [gridZ];

					PdrudeZ[k][s]     = new double [gridZ];
					PdrudeZOld[k][s]  = new double [gridZ];
					PdrudeZOld2[k][s] = new double [gridZ];

				}
			}
			for( k = 0; k < drudePoles; k++)
			{
				for (int s = 0; s<gridR; s++)
				{
					for (m = 0; m < gridZ; m++)
					{
						PdrudeR[k][s][m]     = 0;
						PdrudeROld[k][s][m]  = 0;
						PdrudeROld2[k][s][m] = 0;

						PdrudeZ[k][s][m]     = 0;
						PdrudeZOld[k][s][m]  = 0;
						PdrudeZOld2[k][s][m] = 0;

						PdrudePHI[k][s][m]     = 0;
						PdrudePHIOld[k][s][m]  = 0;
						PdrudePHIOld2[k][s][m] = 0;
					}
				}
			}
		}
		#pragma endregion
		#pragma region//*********** 2LVL ATOM MODEL ************//
		double blochStartR, blochEndR, blochStartZ, blochEndZ;
		double blochTau1, blochTau2, TNatom, gammaBloch;
		double rho30, blochFreqMin, blochFreqMax, blochCm, blochCp, blochD;
		int odeNumber;
		double *omegaBloch = NULL, *Natom = NULL, *blochA = NULL, *blochB = NULL;
		double ***u1r = NULL, ***u2r  = NULL, ***u3r = NULL, ***u1rstore = NULL;
		double ***u1z = NULL, ***u2z  = NULL, ***u3z = NULL, ***u1zstore = NULL;
		double ***u1phi = NULL, ***u2phi  = NULL, ***u3phi = NULL, ***u1phistore = NULL;
		double  **EzStore = NULL, **ErStore = NULL, **EphiStore = NULL;
		double **c1bloch = NULL, **c2bloch = NULL;
		double blochEps, blochSig;

		if (blochOn == "on")
		{
			std::cout << "2LVL Atom Medium on" << endl;
		
			settings.get("LorentzMedium/radialStart", blochStartR);
			settings.get("LorentzMedium/radialEnd", blochEndR);
			settings.get("LorentzMedium/zStart", blochStartZ);
			settings.get("LorentzMedium/zEnd", blochEndZ);

			settings.get("2LVLAtom/ODEnumber", odeNumber);
			settings.get("2LVLAtom/tau1", blochTau1);
			settings.get("2LVLAtom/tau2", blochTau2);
			settings.get("2LVLAtom/Natom", TNatom);
			settings.get("2LVLAtom/gammaN", gammaBloch);
			settings.get("2LVLAtom/minimumFrequency", blochFreqMin);
			settings.get("2LVLAtom/maximumFrequency", blochFreqMax);
			settings.get("2LVLAtom/rho30", rho30);

			settings.get("2LVLAtom/eps", blochEps);
			settings.get("2LVLAtom/sigma", blochSig);

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

			ErStore   = new double* [gridR];
			EzStore   = new double* [gridR];
			EphiStore = new double* [gridR];
			c1bloch = new double* [gridR];
			c2bloch = new double* [gridR];

			for (int s = 0; s<gridR; s++)
			{
			ErStore[s]    = new double [gridZ];
			EzStore[s] = new double [gridZ];
			EphiStore[s] = new double [gridZ];
			c1bloch[s] = new double [gridZ];
			c2bloch[s] = new double [gridZ];
			}

			for (int s = 0; s<gridR; s++)
			{
				for (m = 0; m < gridZ; m++)
				{
				double sigTemp = 0;
				double epsTemp = 1;

				ErStore[s][m]   = 0;
				EzStore[s][m]   = 0;
				EphiStore[s][m] = 0;
				c1bloch[s][m] = (1-(sigTemp*dt)/(2*EPSNOT*epsTemp))/(1+(sigTemp*dt)/(2*EPSNOT*epsTemp));
				c2bloch[s][m] = (dt/(EPSNOT*epsTemp))/(1+(sigTemp*dt)/(2*EPSNOT*epsTemp));
				}
			}

			u1r = new double** [odeNumber];
			u2r = new double** [odeNumber];
			u3r = new double** [odeNumber];
			u1rstore = new double** [odeNumber];
			u1z = new double** [odeNumber];
			u2z = new double** [odeNumber];
			u3z = new double** [odeNumber];
			u1zstore = new double** [odeNumber];
			u1phi = new double** [odeNumber];
			u2phi = new double** [odeNumber];
			u3phi = new double** [odeNumber];
			u1phistore = new double** [odeNumber];
			for (k = 0; k< odeNumber; k++)
			{
				u1r[k] = new double* [gridR];
				u2r[k] = new double* [gridR];
				u3r[k] = new double* [gridR];
				u1rstore[k] = new double* [gridR];
				u1z[k] = new double* [gridR];
				u2z[k] = new double* [gridR];
				u3z[k] = new double* [gridR];
				u1zstore[k] = new double* [gridR];
				u1phi[k] = new double* [gridR];
				u2phi[k] = new double* [gridR];
				u3phi[k] = new double* [gridR];
				u1phistore[k] = new double* [gridR];
			}
			for( k = 0; k < odeNumber; k++)
			{
				for (int s = 0; s<gridR; s++)
				{
					u1r[k][s]        = new double [gridZ];
					u2r[k][s]        = new double [gridZ];
					u3r[k][s]        = new double [gridZ];
					u1rstore[k][s]   = new double [gridZ];
					u1z[k][s]        = new double [gridZ];
					u2z[k][s]        = new double [gridZ];
					u3z[k][s]        = new double [gridZ];
					u1zstore[k][s]   = new double [gridZ];
					u1phi[k][s]      = new double [gridZ];
					u2phi[k][s]      = new double [gridZ];
					u3phi[k][s]      = new double [gridZ];
					u1phistore[k][s] = new double [gridZ];
				}
			}
			for( k = 0; k < odeNumber; k++)
			{
				for (int s = 0; s<gridR; s++)
				{
					for (m = 0; m < gridZ; m++)
					{
						u1r[k][s][m] = 0;
						u2r[k][s][m] = 0;
						u3r[k][s][m] = 0;
						u1rstore[k][s][m] = 0;
						u1z[k][s][m] = 0;
						u2z[k][s][m] = 0;
						u3z[k][s][m] = 0;
						u1zstore[k][s][m] = 0;
						u1phi[k][s][m] = 0;
						u2phi[k][s][m] = 0;
						u3phi[k][s][m] = 0;
						u1phistore[k][s][m] = 0;
					}
				}
			}  

		}
		#pragma endregion
		#pragma region//************* NONLINEARITY *************//
		double alphaN, chi3, NLiter;
		double nonlinStartR, nonlinEndR, nonlinStartZ, nonlinEndZ;
		double ramanTau1, ramanTau2;
		double alphaRaman = NULL, betaRaman = NULL, kappaRaman = NULL;
		
		double **PramanR = NULL, **PramanROld = NULL, **PramanROld2 = NULL;
		double **PramanZ = NULL, **PramanZOld = NULL, **PramanZOld2 = NULL;
		double **PramanPHI = NULL, **PramanPHIOld = NULL, **PramanPHIOld2 = NULL;


		if (nonlinOn == "on")
		{
			cout << "Kerr/Raman Nonlinearity on" << endl;

			settings.get("NonlinearMedium/radialStart", nonlinStartR);
			settings.get("NonlinearMedium/radialEnd", nonlinEndR);
			settings.get("NonlinearMedium/zStart", nonlinStartZ);
			settings.get("NonlinearMedium/zEnd", nonlinEndZ);

			settings.get("NonlinearMedium/alpha", alphaN);
			settings.get("NonlinearMedium/chi3", chi3);
			settings.get("NonlinearMedium/iterations", NLiter);
			settings.get("NonlinearMedium/tau1", ramanTau1);
			settings.get("NonlinearMedium/tau2", ramanTau2);

			PramanR = new double* [gridR];
			PramanROld = new double* [gridR];
			PramanROld2= new double* [gridR];

			PramanZ = new double* [gridR];
			PramanZOld = new double* [gridR];
			PramanZOld2= new double* [gridR];

			PramanPHI = new double* [gridR];
			PramanPHIOld = new double* [gridR];
			PramanPHIOld2= new double* [gridR];

			for (k = 0; k< gridR; k++)
			{
				PramanR[k]		= new double [gridZ];
				PramanROld[k]	= new double [gridZ];
				PramanROld2[k]	= new double [gridZ];

				PramanZ[k]		= new double [gridZ];
				PramanZOld[k]	= new double [gridZ];
				PramanZOld2[k]	= new double [gridZ];

				PramanPHI[k]	= new double [gridZ];
				PramanPHIOld[k]	= new double [gridZ];
				PramanPHIOld2[k] = new double [gridZ];
			}

			for (int s = 0; s<gridR; s++)
			{
				for (m = 0; m < gridZ; m++)
				{
					PramanR[s][m]     = 0;
					PramanROld[s][m]  = 0;
					PramanROld2[s][m] = 0;

					PramanZ[s][m]     = 0;
					PramanZOld[s][m]  = 0;
					PramanZOld2[s][m] = 0;

					PramanPHI[s][m]     = 0;
					PramanPHIOld[s][m]  = 0;
					PramanPHIOld2[s][m] = 0;
				}
			}

		}
		#pragma endregion
		std::cout<< "" << endl;

	#pragma endregion //%%%%%%%%%%%% END INITIALIZE MEDIA %%%%%%%%%%%%%//
	#pragma region //%%%%%%%% DOMAIN PRE-LOOP COMPUTATIONS %%%%%%%%%//
	
	// Domain computations
		
    dx = 0;
	dr = rRange/gridR;
	dz = zRange/gridZ;
	if (M == 0)
		dt = .5/( (1)*(CNOT)*sqrt( 1/pow(dr,2.0)+1/pow(dz,2.0) ) ); 
	else
		dt = .5/( (M)*(CNOT)*sqrt( 1/pow(dr,2.0)+1/pow(dz,2.0) ) ); 
	maxIter = (int) round(timeDom/dt);

	skipSteps = (int) round(maxIter/keepFrames);
	T = 0;

	
	sigx = 0;
	muy = 1;	

	// Start file and input pertinent data
	parameterValues[1] = gridR;
	parameterValues[2] = gridZ;
	parameterValues[3] = keepFrames;
	ofstream data(fileName.c_str(), ios::out | ios::binary | ios::trunc);    // Open file to store data
	for(k = 0; k<PARAM; k++)
		data.write((char *) &parameterValues[k], sizeof parameterValues[k]);
	data.close();

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

	double* ri     = new double [gridR];
	double* riHalf = new double [gridR];
	for (k = 0; k<gridR; k++)
	{
		ri[k] = (k + 0.5)*dr;
		riHalf[k] = k*dr;
	}


	double** dz1 = new double* [gridR];
	double** dz2 = new double* [gridR];
	double*  dz3  = new double [gridZ];

	for (k = 0; k<gridR; k++)
	{
		dz1[k] = new double [gridZ];
		dz2[k] = new double [gridZ];
	}

	double** dr1 = new double* [gridR];
	double** dr2 = new double* [gridR];
	double** dphi1 = new double* [gridR];

	for (k = 0; k<gridR; k++)
	{
		dr1[k] = new double [gridZ-1];
		dr2[k] = new double [gridZ-1];
		dphi1[k] = new double [gridZ-1];
	}

	for (k = 0; k<gridR; k++)
	{
		for( s = 0; s<gridZ; s++ )
		{
			dz1[k][s] = (M*dt)/(MUNOT*muz[k][s]);
			dz2[k][s] = dt/(MUNOT*muz[k][s]);
		}
		for( s = 0; s<gridZ-1; s++ )
		{
			dr1[k][s] = (M*dt)/(MUNOT*mur[k][s]);
			dr2[k][s] = dt/(MUNOT*mur[k][s]);

			dphi1[k][s] = (dt)/(MUNOT*muphi[k][s]);
		}
	}

	k = 0;
	for( s = 0; s < gridZ; s++)
		dz3[s] = (4*dt)/(MUNOT*muz[k][s]);
	


	// Source
	maxT = width*3;
	timeSwitch = maxT*(2);
	lambda = CNOT/f;
	omega = 2*PI*f;

	// Gaussian Beam source calculations
	int gridFocusWidth = (int) round(focusWidth/dr);
	double totalW = gridR-2;
	double zR = PI*pow(focusWidth, 2.0)/lambda;
	double rNot = zNot*(1+pow(zR/zNot, 2.0) );

	double* rVec = new double [gridR];
	double* zVec = new double [gridR];
	double widthBeam = (double) gridFocusWidth*(sqrt(1+pow(zNot/zR,2.0)));
	double rPos = 0;
	for ( k = 0; k<gridR; k++)
	{	
		rVec[k] = rPos;
		rPos = rPos + dr;
	}

	for ( k = 0; k<gridR; k++)
	{	
		double radicand = abs(pow(rNot,2.0) - pow(rVec[k],2.0));
		zVec[k] = - sqrt( radicand ) - rNot;
	}

	int **positionsB = new int *[2];
			 positionsB[0] = new int [gridR];
	         positionsB[1] = new int [gridR];

	for (k = 0; k < gridR; k++)
	{
		positionsB[0][k] = (int) round(rVec[k]/dr);
		positionsB[1][k] = (int) round( zVec[k]/dr + sourcePos );
	}
	// end Gaussian Beam source calculations

	// PML computations
	sigMax = maxC*(sigM+1)*.8/(ETANOT*dr);   // Max value of sigma in PML

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
		 
		bn[k] = exp(-dt/tau);
		bs[k] = exp(-dt/tauStag);

		gn[k] = (sigVec/(kappaVec*(kappaVec*aVec+sigVec)))*(1-bn[k]);
		gs[k] = (sigVecStag/(kappaVecStag*(kappaVecStag*aVecStag+sigVecStag)))*(1-bs[k]);

		kn[k] = kappaVec;
		ks[k] = kappaVecStag;

		/////////////////////////////////////////////////////////////////////////////////////////////
		/*
		double rMax = ri[gridR-1];
		double rVal = ri[gridR-PML-1+k];
		double rVal2 = riHalf[gridR-PML-1+k];

		double sigVec2		= sigMax * pow((rVal/rMax)*(xVec + .5), sigM)/(pow(pmlWidth, sigM));
		double sigVecStag2  = sigMax * pow((rVal2/rMax)*xVec, sigM)/(pow(pmlWidth, sigM));

		double kappaVec2	= 1 + (kappaMax-1)* pow(xVec + .5, kappaM)/(pow(pmlWidth, kappaM));
		double kappaVecStag2= 1 + (kappaMax-1)* pow(xVec , kappaM)/(pow(pmlWidth, kappaM));

		double aVec2		= aMax * pow(xVec + .5, aM)/(pow(pmlWidth, aM));
		double aVecStag2	= aMax * pow(xVec, aM)/(pow(pmlWidth, aM));

		double tau2			= (kappaVec*EPSNOT)/(kappaVec*aVec + sigVec);
		double tauStag2	= (kappaVecStag*EPSNOT)/(kappaVecStag*aVecStag + sigVecStag);

		bn2[k] = exp(-dt/tau2);
		bs2[k] = exp(-dt/tauStag2);

		gn2[k] = (sigVec2/(kappaVec2*(kappaVec2*aVec2+sigVec2)))*(1-bn2[k]);
		gs2[k] = (sigVecStag2/(kappaVecStag2*(kappaVecStag2*aVecStag2+sigVecStag2)))*(1-bs2[k]);

		kn2[k] = kappaVec2;
		ks2[k] = kappaVecStag2;
		*/

	}


	#pragma endregion //%%%%%% END DOMAIN PRE-LOOP COMPUTATIONS %%%%%%%//
	#pragma region //%%%%%%%%% MEDIA PRE-LOOP COMPUTATIONS %%%%%%%%%//

		#pragma region //************** DIELECTRIC **************//
		if (dielectricOn == "on")
		{
			dielStartZ = (int) round(dielStartZ /dz);
			dielStartR = (int) round(dielStartR /dr);
			dielEndZ   = (int) round(dielEndZ/dz);
			dielEndR   = (int) round(dielEndR/dr);
			for ( k = (int) dielStartR ; k < (int) dielEndR; k++ )
			{	
				for ( s = (int) dielStartZ ; s < (int) dielEndZ; s++ )
				{ 
				epsr[k][s] = epsRel;
				epsz[k][s] = epsRel;
				epsphi[k][s] = epsRel;
				}
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

			lorzStartZ = (int) round(lorzStartZ/dz);
			lorzEndZ   = (int) round(lorzEndZ/dz);
			lorzStartR = (int) round(lorzStartR/dr);
			lorzEndR   = (int) round(lorzEndR/dr);

			for ( k = (int) lorzStartR ; k < (int) lorzEndR; k++ )
			{	
				for ( s = (int) lorzStartZ ; s < (int) lorzEndZ; s++ )
				{ 
				epsr[k][s] = epsInfLorentz;
				epsz[k][s] = epsInfLorentz;
				epsphi[k][s] = epsInfLorentz;
				}
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
			debyeStartZ = (int) round(debyeStartZ/dz);
			debyeEndZ   = (int) round(debyeEndZ/dz);
			debyeStartR = (int) round(debyeStartR/dr);
			debyeEndR   = (int) round(debyeEndR/dr);

		for ( k = (int) debyeStartR ; k < (int) debyeEndR; k++ )
			{	
				for ( s = (int) debyeStartZ ; s < (int) debyeEndZ; s++ )
				{ 
				epsr[k][s] = epsInfDebye;
				epsz[k][s] = epsInfDebye;
				epsphi[k][s] = epsInfDebye;
				}
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
			drudeStartZ = (int) round(drudeStartZ/dz);
			drudeEndZ   = (int) round(drudeEndZ/dz);
			drudeStartR = (int) round(drudeStartR/dr);
			drudeEndR   = (int) round(drudeEndR/dr);
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
			blochStartZ = (int) round(blochStartZ/dz);
			blochEndZ   = (int) round(blochEndZ/dz);
			blochStartR = (int) round(blochStartR/dr);
			blochEndR   = (int) round(blochEndR/dr);

			for (int s = (int) blochStartR; s < (int) blochEndR; s++)
			{
				for (m = (int) blochStartZ; m < (int) blochEndZ; m++)
				{
				c1bloch[s][m] = (1-(blochSig*dt)/(2*EPSNOT*blochEps))/(1+(blochSig*dt)/(2*EPSNOT*blochEps));
				c2bloch[s][m] = (dt/(EPSNOT*blochEps))/(1+(blochSig*dt)/(2*EPSNOT*blochEps));
				}
			}

		}
		#pragma endregion
		#pragma region//********** RAMAN NONLINEARITY **********//
		if (nonlinOn == "on")
		{
			double omegaRam = sqrt( ( pow(ramanTau1,2.0) + pow(ramanTau2,2.0) ) / ( pow(ramanTau1,2.0) * pow(ramanTau2,2.0) ) );
			double deltaRam = 1/ramanTau2;

			alphaRaman = ((1-alphaN)*chi3*pow(omegaRam,2.0)*pow(dt,2.0))/(deltaRam*dt+1);
			betaRaman  = (2-pow(dt,2.0)*pow(omegaRam,2.0))/(deltaRam*dt+1);
			kappaRaman = (deltaRam*dt-1)/(deltaRam*dt+1);

			nonlinStartZ = (int) round(nonlinStartZ /dz);
			nonlinStartR = (int) round(nonlinStartR /dr);
			nonlinEndZ   = (int) round(nonlinEndZ/dz);
			nonlinEndR   = (int) round(nonlinEndR/dr);
		}
		#pragma endregion
		std::cout << "Beginning Calculations ..." << endl;
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

		if( pulseSetting == "gaussBeamCW" )
		{
			for (k = 0; k<gridR-1; k++)
			{
				source = pAmp* (gridFocusWidth/widthBeam) * exp( - pow(k,2.0) / pow(widthBeam, 2.0) ) * sin(omega*double(T)*dt);
				int pos1 = positionsB[0][k];
				int pos2 = positionsB[1][k];
				Hphi[pos1][pos2] = source*2.105e-3 + Hphi[pos1][pos2];
				Dphi[pos1][pos2] = source*EPSNOT + Dphi[pos1][pos2];
			}
		}

		if( pulseSetting == "gaussBeamPulse" )
		{
			for (k = 1; k<gridR-1; k++)
			{
				double gauss = exp( -1.0*pow(maxT-(double)T*dt,2.0) / pow(width/5.0, 2.0) ) ;
				source = pAmp* (gridFocusWidth/widthBeam) * exp( - pow(k,2.0) / pow(widthBeam, 2.0) ) *gauss* sin(omega*double(T)*dt);
				int pos1 = positionsB[0][k];
				int pos2 = positionsB[1][k];
				Hphi[pos1][pos2] = source*2.105e-3 + Hphi[pos1][pos2];
				Dphi[pos1][pos2] = source*EPSNOT + Dphi[pos1][pos2];

			}
		}

		# pragma endregion //%%%%%%%%%%%%%%%%% END SOURCE %%%%%%%%%%%%%%%%%%//
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEDIUM UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%%%% LORENTZ MEDIUM UPDATE %%%%%%%%%%%%//
		if (lorentzOn == "on")
		{

			for( k = 0; k < lorentzPoles; k++)
			{
				for (int s = (int) lorzStartR; s < (int) lorzEndR; s++)
				{
					for (m = (int) lorzStartZ; m < (int) lorzEndZ; m++)
					{
						PlorzR[k][s][m]  = alphaLorz[k]*Er[s][m] + betaLorz[k]*PlorzR[k][s][m]  + gammaLorz[k]*PlorzROld2[k][s][m] ;
						PLr[s][m] = PLr[s][m] + PlorzR[k][s][m]   ;
						PlorzROld2[k][s][m]  = PlorzROld[k][s][m] ;
						PlorzROld[k][s][m]  = PlorzR[k][s][m]     ;
			
						PlorzZ[k][s][m]  = alphaLorz[k]*Ez[s][m] + betaLorz[k]*PlorzZ[k][s][m]  + gammaLorz[k]*PlorzZOld2[k][s][m] ;
						PLz[s][m] = PLz[s][m] + PlorzZ[k][s][m]   ;
						PlorzZOld2[k][s][m]  = PlorzZOld[k][s][m] ;
						PlorzZOld[k][s][m]  = PlorzZ[k][s][m]     ;

						PlorzPHI[k][s][m]  = alphaLorz[k]*Ephi[s][m] + betaLorz[k]*PlorzPHI[k][s][m]  + gammaLorz[k]*PlorzPHIOld2[k][s][m] ;
						PLphi[s][m] = PLphi[s][m] + PlorzPHI[k][s][m]   ;
						PlorzPHIOld2[k][s][m]  = PlorzPHIOld[k][s][m] ;
						PlorzPHIOld[k][s][m]  = PlorzPHI[k][s][m]     ;
					}
				}
			}
		}
		
		#pragma endregion //%%%%%%%%%% END LORENTZ MEDIUM UPDATE %%%%%%%%%%//
		#pragma region //%%%%%%%%%%%%% DEBYE MEDIUM UPDATE %%%%%%%%%%%%%//
		if (debyeOn == "on")
		{
			for( k = 0; k < debyePoles; k++)
			{
				for (int s = (int) debyeStartR; s < (int) debyeEndR; s++)
				{
					for (m = (int) debyeStartZ; m < (int) debyeEndZ; m++)
					{
						PdebyeR[k][s][m]  = alphaDebye[k]*Er[s][m] + betaDebye[k]*PdebyeR[k][s][m]  + PdebyeROld2[k][s][m] ;
						PLr[s][m] = PLr[s][m] + PdebyeR[k][s][m]   ;
						PdebyeROld2[k][s][m]  = PdebyeROld[k][s][m] ;
						PdebyeROld[k][s][m]  = PdebyeR[k][s][m]     ;
			
						PdebyeZ[k][s][m]  = alphaDebye[k]*Ez[s][m] + betaDebye[k]*PdebyeZ[k][s][m]  + PdebyeZOld2[k][s][m] ;
						PLz[s][m] = PLz[s][m] + PdebyeZ[k][s][m]   ;
						PdebyeZOld2[k][s][m]  = PdebyeZOld[k][s][m] ;
						PdebyeZOld[k][s][m]  = PdebyeZ[k][s][m]     ;

						PdebyePHI[k][s][m]  = alphaDebye[k]*Ephi[s][m] + betaDebye[k]*PdebyePHI[k][s][m]  + PdebyePHIOld2[k][s][m] ;
						PLphi[s][m] = PLphi[s][m] + PdebyePHI[k][s][m]   ;
						PdebyePHIOld2[k][s][m]  = PdebyePHIOld[k][s][m] ;
						PdebyePHIOld[k][s][m]  = PdebyePHI[k][s][m]     ;
					}
				}
			}
		}
		#pragma endregion //%%%%%%%%%%% END DEBYE MEDIUM UPDATE %%%%%%%%%%%//
		#pragma region //%%%%%%%%%%%%% DRUDE MEDIUM UPDATE %%%%%%%%%%%%%//
		if (drudeOn == "on")
		{			
			for( k = 0; k < lorentzPoles; k++)
			{
				for (int s = (int) drudeStartR; s < (int) drudeEndR; s++)
				{
					for (m = (int) drudeStartZ; m < (int) drudeEndZ; m++)
					{
						PdrudeR[k][s][m]  = alphaDrude[k]*Er[s][m] + betaDrude[k]*PdrudeR[k][s][m]  + kappaDrude[k]*PdrudeROld2[k][s][m] ;
						PLr[s][m] = PLr[s][m] + PdrudeR[k][s][m]   ;
						PdrudeROld2[k][s][m]  = PdrudeROld[k][s][m] ;
						PdrudeROld[k][s][m]  = PdrudeR[k][s][m]     ;
			
						PdrudeZ[k][s][m]  = alphaDrude[k]*Ez[s][m] + betaDrude[k]*PdrudeZ[k][s][m]  + kappaDrude[k]*PdrudeZOld2[k][s][m] ;
						PLz[s][m] = PLz[s][m] + PdrudeZ[k][s][m]   ;
						PdrudeZOld2[k][s][m]  = PdrudeZOld[k][s][m] ;
						PdrudeZOld[k][s][m]  = PdrudeZ[k][s][m]     ;

						PdrudePHI[k][s][m]  = alphaDrude[k]*Ephi[s][m] + betaDrude[k]*PdrudePHI[k][s][m]  + kappaDrude[k]*PdrudePHIOld2[k][s][m] ;
						PLphi[s][m] = PLphi[s][m] + PdrudePHI[k][s][m]   ;
						PdrudePHIOld2[k][s][m]  = PdrudePHIOld[k][s][m] ;
						PdrudePHIOld[k][s][m]  = PdrudePHI[k][s][m]     ;
					}
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


			for( k = 0; k < odeNumber; k++)
			{

				blochA[k] = (  (Natom[k]*gammaBloch)/(EPSNOT*blochTau2)  ) * exp(-tHalf/blochTau2);
				blochB[k] = ( (Natom[k]*gammaBloch*omegaBloch[k])/EPSNOT ) * exp(-tHalf/blochTau2);
				double denom = ( 1 + (pow(omegaBloch[k],2.0)*pow(dt,2.0))/4 );
				double cof1 = dt*(blochCm/2);
				double coef1 = 1 - ( pow(omegaBloch[k],2.0)*pow(dt,2.0) )/4 ;
				double coef2 = omegaBloch[k]*dt;
				double coef3 = omegaBloch[k]*pow(dt,2.0)/2;
				double coef4 = omegaBloch[k]*dt/2;
	
				for (int s = (int) blochStartR; s < (int) blochEndR; s++)
				{
					for (m = (int) blochStartZ; m < (int) blochEndZ; m++)
					{
					u3r[k][s][m] = u3r[k][s][m]  - (dt*blochCm/2)*(Er[s][m]+ ErStore[s][m])*u2r[k][s][m] ;
					u3z[k][s][m] = u3z[k][s][m]  - (dt*blochCm/2)*(Ez[s][m]+ EzStore[s][m])*u2z[k][s][m] ;
					u3phi[k][s][m] = u3phi[k][s][m]  - (dt*blochCm/2)*(Ephi[s][m]+ EphiStore[s][m])*u2phi[k][s][m] ;

					double Qr = blochCp*u3r[k][s][m] *Er[s][m]  +  blochD*Er[s][m];
					double Qz = blochCp*u3z[k][s][m] *Ez[s][m]  +  blochD*Ez[s][m];
					double Qphi = blochCp*u3phi[k][s][m] *Ephi[s][m]  +  blochD*Ephi[s][m];

					u1r[k][s][m] = ( coef1*u1r[k][s][m]  + coef2* u2r[k][s][m]   + coef3*Qr ) / denom;
					u2r[k][s][m] = u2r[k][s][m] - coef4*( u1r[k][s][m]  + u1rstore[k][s][m] ) + dt*Qr;

					u1z[k][s][m] = ( coef1*u1z[k][s][m]  + coef2* u2z[k][s][m]   + coef3*Qz ) / denom;
					u2z[k][s][m] = u2z[k][s][m] - coef4*( u1z[k][s][m]  + u1zstore[k][s][m] ) + dt*Qz;

					u1phi[k][s][m] = ( coef1*u1phi[k][s][m]  + coef2* u2phi[k][s][m]   + coef3*Qphi ) / denom;
					u2phi[k][s][m] = u2phi[k][s][m] - coef4*( u1phi[k][s][m]  + u1phistore[k][s][m] ) + dt*Qphi;

					PLr[s][m] = PLr[s][m] + dt*(-blochA[k]*u1r[k][s][m]  + blochB[k]*u2r[k][s][m] );
					PLz[s][m] = PLz[s][m] + dt*(-blochA[k]*u1z[k][s][m]  + blochB[k]*u2z[k][s][m] );
					PLphi[s][m] = PLphi[s][m] + dt*(-blochA[k]*u1phi[k][s][m]  + blochB[k]*u2phi[k][s][m] );

					u1rstore[k][s][m]   = u1r[k][s][m] ;
					u1zstore[k][s][m]   = u1z[k][s][m] ;
					u1phistore[k][s][m] = u1phi[k][s][m] ;
					}
				}
			}

			for (int s = (int) blochStartR; s < (int) blochEndR; s++)
				{
					for (m = (int) blochStartZ; m < (int) blochEndZ; m++)
					{
						ErStore[s][m] = Er[s][m];
						EzStore[s][m] = Ez[s][m];
						EphiStore[s][m] = Ephi[s][m];
					}
			}
		}
		#pragma endregion //%%%%%%%%% END 2LVL ATOM MEDIUM UPDATE %%%%%%%%%//
		#pragma region //%%%%%%%%%%%%% RAMAN MEDIUM UPDATE %%%%%%%%%%%%%//
		if (nonlinOn == "on")
		{
			for (int s = (int) nonlinStartR; s < (int) nonlinEndR; s++)
			{
				for (m = (int) nonlinStartZ; m < (int) nonlinEndZ; m++)
				{
				PramanR[s][m]   = alphaRaman*pow(Er[s][m],2.0) + betaRaman*PramanR[s][m] + kappaRaman*PramanROld2[s][m];
				PramanZ[s][m]   = alphaRaman*pow(Ez[s][m],2.0) + betaRaman*PramanZ[s][m] + kappaRaman*PramanZOld2[s][m];
				PramanPHI[s][m] = alphaRaman*pow(Ephi[s][m],2.0) + betaRaman*PramanPHI[s][m] + kappaRaman*PramanPHIOld2[s][m];

				PramanROld2[s][m]   = PramanROld[s][m];
				PramanZOld2[s][m]   = PramanZOld[s][m];
				PramanPHIOld2[s][m] = PramanPHIOld[s][m];

				PramanROld[s][m]   = PramanR[s][m];
				PramanZOld[s][m]   = PramanZ[s][m];
				PramanPHIOld[s][m] = PramanPHI[s][m];
				}
			}
		}
		#pragma endregion //%%%%%%%%%%% END RAMAN MEDIUM UPDATE %%%%%%%%%%%//
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MEDIUM UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIELD UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%%%%%%% H-FIELD UPDATE %%%%%%%%%%%%%%%%//

		//%%%%%%%%%% Hr FIELD %%%%%%%%%%
		for ( k = 0; k < gridR-PML; k++ )
		{
			for ( s = PML; s < gridZ - PML - 1; s++)
			{
				curlEphi = (Ephi[k][s+1] - Ephi[k][s])/dz;
				Hr[k][s] = Hr[k][s] - dr1[k][s]*Ez[k][s]/ri[k] + dr2[k][s]*curlEphi;
			}
		}

		//%%%%%%%%%% Hz FIELD %%%%%%%%%%
		for ( k = 1; k < gridR-PML; k++ )
		{
			for ( s = PML; s < gridZ - PML; s++)
			{
				curlEphi = (ri[k]*Ephi[k][s] - ri[k-1]*Ephi[k-1][s])/dr;
				Hz[k][s] = Hz[k][s] + dz1[k][s]*Er[k][s]/riHalf[k] - dz2[k][s]*curlEphi/riHalf[k];
			}

		}
		if (M == 0)
		{
			for ( s = PML; s < gridZ - PML; s++)
			{
				Hz[0][s] = Hz[0][s] - dz3[s]*Ephi[0][s]/dr;
			}
		}
		else
		{
			for ( s = PML; s < gridZ - PML; s++)
				Hz[0][s] = 0;
		}

		//%%%%%%%%%% Hphi FIELD %%%%%%%%%%
		for ( k = 1; k < gridR-PML; k++ )
		{
			for ( s = PML; s < gridZ - PML - 1; s++)
			{
				curlEr = (Er[k][s+1] - Er[k][s])/dz;
				curlEz = (Ez[k][s] - Ez[k-1][s])/dr;
				Hphi[k][s] = Hphi[k][s] - dphi1[k][s]*( curlEr - curlEz);
			}
		}

		for ( s = PML; s < gridZ - PML - 1; s++)
				Hphi[0][s] = 0;
				

		
		#pragma endregion //%%%%%%%%%%%%% END H-FIELD UPDATE %%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%%%% LINEAR E-FIELD UPDATE %%%%%%%%%%%%//
		if (nonlinOn != "on"  && blochOn != "on" )
		{
			
			//%%%%%%%%%% Er FIELD %%%%%%%%%%
			for ( k = 1 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHphi  = (Hphi[k][s] - Hphi[k][s-1])/dz;
					Dr[k][s]  = Dr[k][s] - dt*(curlHphi + M*Hz[k][s]/riHalf[k]);
					Er[k][s]  = ( Dr[k][s] - PLr[k][s] )/( EPSNOT*epsr[k][s] );
					PLr[k][s] = 0;
				}
			}

			//%%%%%%%%%% Ez FIELD %%%%%%%%%%
			for ( k = 0 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHphi  = (riHalf[k+1]*Hphi[k+1][s] - riHalf[k]*Hphi[k][s])/dr;
					Dz[k][s]  = Dz[k][s] + dt*(curlHphi + M*Hr[k][s])/ri[k];
					Ez[k][s]  = ( Dz[k][s] - PLz[k][s] )/( EPSNOT*epsz[k][s] );
					PLz[k][s] = 0;
				}
			}
			
			//%%%%%%%%%% Ephi FIELD %%%%%%%%%%
			for ( k = 0 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHr      = ( Hr[k][s] - Hr[k][s-1] )/dz;
					curlHz      = ( Hz[k+1][s] - Hz[k][s] )/dr;
					Dphi[k][s]  = Dphi[k][s] + dt*(curlHr - curlHz);
					Ephi[k][s]  = ( Dphi[k][s] - PLphi[k][s] )/( EPSNOT*epsphi[k][s] );
					PLphi[k][s] = 0;
				}
			}


		}
		#pragma endregion //%%%%%%%%%% END LINEAR E-FIELD UPDATE %%%%%%%%%%//
		#pragma region //%%%%%%%%%% NONLINEAR E-FIELD UPDATE %%%%%%%%%%%//
		if (nonlinOn == "on" )
		{

			//%%%%%%%%%% Er FIELD %%%%%%%%%%
			for ( k = 1 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHphi    = (Hphi[k][s] - Hphi[k][s-1])/dz;
					Dr[k][s]    = Dr[k][s] - dt*(curlHphi + M*Hz[k][s]/riHalf[k]);
					Erlin[k][s] = ( Dr[k][s] - PLr[k][s] );
					ErNL[k][s]  = Er[k][s];
					PLr[k][s]   = 0;
				}
			}

			//%%%%%%%%%% Ez FIELD %%%%%%%%%%
			for ( k = 0 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHphi    = (riHalf[k+1]*Hphi[k+1][s] - riHalf[k]*Hphi[k][s])/dr;
					Dz[k][s]    = Dz[k][s] + dt*(curlHphi + M*Hr[k][s])/ri[k];
					Ezlin[k][s] = ( Dz[k][s] - PLz[k][s] );
					EzNL[k][s]  = Ez[k][s];
					PLz[k][s]   = 0;
				}
			}
			
			//%%%%%%%%%% Ephi FIELD %%%%%%%%%%
			for ( k = 0 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHr        = ( Hr[k][s] - Hr[k][s-1] )/dz;
					curlHz        = ( Hz[k+1][s] - Hz[k][s] )/dr;
					Dphi[k][s]    = Dphi[k][s] + dt*(curlHr - curlHz);
					Ephilin[k][s] = ( Dphi[k][s] - PLphi[k][s] );
					EphiNL[k][s]  = Ephi[k][s];
					PLphi[k][s]   = 0;
				}
			}


			double counter = 0;
			
			while (counter < NLiter)
			{
				for ( k = 0 ; k < gridR-PML; k++ )
				{
					for ( s = PML; s < gridZ - PML; s++ )
					{
						Er[k][s] = Erlin[k][s] /( EPSNOT*(epsr[k][s] + PramanR[k][s] + alphaN*chi3*pow(ErNL[k][s],2.0) ) );
						Ez[k][s] = Ezlin[k][s] /( EPSNOT*(epsz[k][s] + PramanZ[k][s] + alphaN*chi3*pow(EzNL[k][s],2.0) ) );
						Ephi[k][s] = Ephilin[k][s] /( EPSNOT*(epsphi[k][s] + PramanPHI[k][s] + alphaN*chi3*pow(EphiNL[k][s],2.0) ) );
						ErNL[k][s] = Er[k][s];
						EzNL[k][s] = Ez[k][s];
						EphiNL[k][s] = Ephi[k][s];
					}
				}
				counter = counter + 1;
			}
			
			
		}

		if (blochOn == "on"	)
		{
					
			//%%%%%%%%%% Er FIELD %%%%%%%%%%
			for ( k = 1 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHphi  = (Hphi[k][s] - Hphi[k][s-1])/dz;
					Er[k][s]  = c1bloch[k][s]*Er[k][s] - c2bloch[k][s] *( curlHphi + M*Hz[k][s]/riHalf[k]) + PLr[k][s];
					PLr[k][s] = 0;
				}
			}

			//%%%%%%%%%% Ez FIELD %%%%%%%%%%
			for ( k = 0 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHphi  = (riHalf[k+1]*Hphi[k+1][s] - riHalf[k]*Hphi[k][s])/dr;
					Ez[k][s]  = c1bloch[k][s]*Ez[k][s] + c2bloch[k][s] *(curlHphi + M*Hr[k][s])/ri[k] + PLz[k][s];
					PLz[k][s] = 0;
				}
			}
			
			//%%%%%%%%%% Ephi FIELD %%%%%%%%%%
			for ( k = 0 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHr      = ( Hr[k][s] - Hr[k][s-1] )/dz;
					curlHz      = ( Hz[k+1][s] - Hz[k][s] )/dr;
					Ephi[k][s]  = c1bloch[k][s]*Ephi[k][s] + c2bloch[k][s] *(curlHr - curlHz) + PLphi[k][s];
					PLphi[k][s] = 0;
				}
			}
		}
		#pragma endregion //%%%%%%%% NONEND LINEAR E-FIELD UPDATE %%%%%%%%%//
		#pragma region //%%%%%%%%%%%%%% PML REGION UPDATE %%%%%%%%%%%%%%//

		#pragma region //%%%%%%%%%% Hr FIELD %%%%%%%%%%//
		for ( k = 0; k < gridR; k++ )
		{
			for ( s = 0; s < PML ; s++)
			{
				curlEphi = (Ephi[k][s+1] - Ephi[k][s])/dz;
				QEzphi[k][s] = bs[s]*QEzphi[k][s] - gs[s]*curlEphi;
				Hr[k][s] = Hr[k][s] - dr1[k][s]*Ez[k][s]/ri[k] + dr2[k][s]*ks[s]*( curlEphi + QEzphi[k][s] );
			}
		}
		for ( k = 0; k < gridR; k++ )
		{
			int counter = PML-1;
			for ( s = gridZ-PML-1; s < gridZ-1 ; s++)
			{
				curlEphi = (Ephi[k][s+1] - Ephi[k][s])/dz;
				QEzphi[k][s] = bs[counter]*QEzphi[k][s]- gs[counter]*curlEphi;
				Hr[k][s] = Hr[k][s] - dr1[k][s]*Ez[k][s]/ri[k] + dr2[k][s]*ks[counter]*( curlEphi + QEzphi[k][s] );
				counter--;
			}
		}
		for ( k = gridR-PML; k < gridR; k++ )
		{
			for ( s = PML; s < gridZ - PML - 1; s++)
			{
				curlEphi = (Ephi[k][s+1] - Ephi[k][s])/dz;
				Hr[k][s] = Hr[k][s] - dr1[k][s]*Ez[k][s]/ri[k] + dr2[k][s]*curlEphi;
			}
		}  
		# pragma endregion //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%% Hz FIELD %%%%%%%%%%//
		for ( k = 1; k < gridR-PML; k++ )
		{
			for ( s = 0; s < PML ; s++)
			{
				curlEphi = (ri[k]*Ephi[k][s] - ri[k-1]*Ephi[k-1][s])/dr;
				Hz[k][s] = Hz[k][s] + dz1[k][s]*Er[k][s]/riHalf[k] - dz2[k][s]*curlEphi/riHalf[k];
			}
		}
		for ( k = 1; k < gridR-PML; k++ )
		{
			for ( s = gridZ-PML; s < gridZ ; s++)
			{
				curlEphi = (ri[k]*Ephi[k][s] - ri[k-1]*Ephi[k-1][s])/dr;
				Hz[k][s] = Hz[k][s] + dz1[k][s]*Er[k][s]/riHalf[k] - dz2[k][s]*curlEphi/riHalf[k];
			}
		}
		m = PML-1;
		for ( k = gridR-PML; k < gridR; k++ )
		{
			
			for ( s = 0; s < gridZ ; s++)
			{
					
				curlEphi = (ri[k]*Ephi[k][s] - ri[k-1]*Ephi[k-1][s])/dr;
				QErphi[k][s] = bs[m]*QErphi[k][s] - gs[m]*curlEphi;
				Hz[k][s] = Hz[k][s] + dz1[k][s]*Er[k][s]/riHalf[k] - dz2[k][s]*ks[m]*( curlEphi + QErphi[k][s]) /riHalf[k];
			}
			m--;
		} 
		# pragma endregion //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%% Hphi FIELD %%%%%%%%%//
		
		for ( k = 1; k < gridR; k++ )
		{
			for ( s = 0; s < PML ; s++)
			{
				curlEr = (Er[k][s+1] - Er[k][s])/dz;
				QEzr[k][s] = bs[s]*QEzr[k][s] - gs[s]*curlEr;
			}
		}

		for ( k = 1; k < gridR; k++ )
		{
			int counter = PML-1;
			for ( s = gridZ-PML-1; s < gridZ-1 ; s++)
			{
				curlEr = (Er[k][s+1] - Er[k][s])/dz;
				QEzr[k][s] = bs[counter]*QEzr[k][s]- gs[counter]*curlEr;
				counter = counter - 1;
			}
		}

		m = PML-1;
		for ( k = gridR-PML; k < gridR; k++ )
		{
			for ( s = 0; s < gridZ; s++)
			{
				curlEz = (Ez[k][s] - Ez[k-1][s])/dr;
				QErz[k][s] = bs[m]*QErz[k][s] - gs[m]*curlEz;
			}
			m--;
		}

		///////////////////////////
		for ( k = 1; k < gridR; k++ )
		{
			for ( s = 0; s < PML ; s++)
			{
				curlEr = (Er[k][s+1] - Er[k][s])/dz;
				curlEz = (Ez[k][s] - Ez[k-1][s])/dr;
				Hphi[k][s] = Hphi[k][s] - dphi1[k][s]*( ( curlEr + QEzr[k][s] ) - ( curlEz + QErz[k][s] ) );
			}
		}
		
		for ( k = 1; k < gridR; k++ )
		{
			for ( s = gridZ-PML-1; s < gridZ-1 ; s++)
			{
				curlEr = (Er[k][s+1] - Er[k][s])/dz;
				curlEz = (Ez[k][s] - Ez[k-1][s])/dr;
				Hphi[k][s] = Hphi[k][s] - dphi1[k][s]*( ( curlEr + QEzr[k][s] ) - ( curlEz + QErz[k][s] ) );
			}
		}

		for ( k = gridR-PML; k < gridR; k++ )
		{
			for (s = PML; s < gridZ-PML-1 ; s++)
			{
				curlEr = (Er[k][s+1] - Er[k][s])/dz;
				curlEz = (Ez[k][s] - Ez[k-1][s])/dr;
				Hphi[k][s] = Hphi[k][s] - dphi1[k][s]*( ( curlEr + QEzr[k][s] ) - ( curlEz + QErz[k][s] ) );
			}
		}
		# pragma endregion //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		
		#pragma region //%%%%%%%%%% Er FIELD %%%%%%%%%%//
		for ( k = 1; k < gridR; k++ )
		{
			for ( s = 0; s < PML ; s++)
			{
				curlHphi  = (Hphi[k][s] - Hphi[k][s-1])/dz;
				QHzphi[k][s] = bn[s]*QHzphi[k][s] - gn[s]*curlHphi;
				double coef = dt/(EPSNOT*epsr[k][s]);
				Er[k][s]  = Er[k][s] - coef*( ( curlHphi + QHzphi[k][s] + M*Hz[k][s]/riHalf[k]) );
			}
		}
		for ( k = 1; k < gridR; k++ )
		{
			int counter = PML-1;
			for ( s = gridZ-PML; s < gridZ ; s++)
			{
				curlHphi  = (Hphi[k][s] - Hphi[k][s-1])/dz;
				QHzphi[k][s] = bn[counter]*QHzphi[k][s]- gn[counter]*curlHphi;
				double coef = dt/(EPSNOT*epsr[k][s]);
				Er[k][s]  = Er[k][s] - coef*( ( curlHphi + QHzphi[k][s] + M*Hz[k][s]/riHalf[k]) );
				counter--;
			}
		}
		for ( k = gridR-PML; k < gridR; k++ )
		{
			for ( s = PML; s < gridZ - PML ; s++)
			{
				curlHphi  = (Hphi[k][s] - Hphi[k][s-1])/dz;
				double coef = dt/(EPSNOT*epsr[k][s]);
				Er[k][s]  = Er[k][s] - coef*( ( curlHphi + QHzphi[k][s] + M*Hz[k][s]/riHalf[k]) );
			}
		}
		# pragma endregion //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%%% Ez FIELD %%%%%%%%%%//
		for ( k = 0; k < gridR-PML; k++ )
		{
			for ( s = 0; s < PML ; s++)
			{
				curlHphi  = (riHalf[k+1]*Hphi[k+1][s] - riHalf[k]*Hphi[k][s])/dr;
				double coef = dt/(EPSNOT*epsz[k][s]);
				Ez[k][s] = Ez[k][s] + coef*( (curlHphi + Hr[k][s])/ri[k] );
			}
		}

		for ( k = 0; k < gridR-PML; k++ )
		{
			for ( s = gridZ-PML; s < gridZ ; s++)
			{
				curlHphi  = (riHalf[k+1]*Hphi[k+1][s] - riHalf[k]*Hphi[k][s])/dr;
				double coef = dt/(EPSNOT*epsz[k][s]);
				Ez[k][s] = Ez[k][s] + coef*( (curlHphi + M*Hr[k][s])/ri[k] );
			}
		}

		m = PML-1;
		for ( k = gridR-PML; k < gridR-1; k++ )
		{
			
			for ( s = 0; s < gridZ ; s++)
			{
				curlHphi  = (riHalf[k+1]*Hphi[k+1][s] - riHalf[k]*Hphi[k][s])/dr;
				double coef = dt/(EPSNOT*epsz[k][s]);
				QHrphi[k][s] = bn[m]*QHrphi[k][s] - gn[m]*curlHphi;
				Ez[k][s] = Ez[k][s] + coef*( (curlHphi + QHrphi[k][s] + M*Hr[k][s])/ri[k] );
			}
			m--;
		}
		# pragma endregion //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		#pragma region //%%%%%%%%% Ephi FIELD %%%%%%%%%//
		
		for ( k = 0; k < gridR; k++ )
		{
			for ( s = 0; s < PML ; s++)
			{
				curlHr = ( Hr[k][s] - Hr[k][s-1] )/dz;
				QHzr[k][s] = bn[s]*QHzr[k][s] - gn[s]*curlHr;
			}
		}

		for ( k = 0; k < gridR; k++ )
		{
			int counter = PML-1;
			for ( s = gridZ-PML; s < gridZ-1 ; s++)
			{
				curlHr = ( Hr[k][s] - Hr[k][s-1] )/dz;
				QHzr[k][s] = bn[counter]*QHzr[k][s]- gn[counter]*curlHr;
				counter = counter - 1;
			}
		}

		m = PML-1;
		for ( k = gridR-PML; k < gridR-1; k++ )
		{
			for ( s = 0; s < gridZ; s++)
			{
				curlHz = ( Hz[k+1][s] - Hz[k][s] )/dr;
				QHrz[k][s] = bn[m]*QHrz[k][s] - gn[m]*curlHz;
			}
			m--;
		}

		
		///////////////////////////
		for ( k = 0; k < gridR-1; k++ )
		{
			for ( s = 0; s < PML ; s++)
			{
				curlHr      = ( Hr[k][s] - Hr[k][s-1] )/dz;
				curlHz      = ( Hz[k+1][s] - Hz[k][s] )/dr;
				double coef = dt/(EPSNOT*epsphi[k][s]);
				Ephi[k][s]  = Ephi[k][s] + coef*( ( curlHr + QHzr[k][s] ) - ( curlHz  + QHrz[k][s] ) );
			}
		}
		
		for ( k = 0; k < gridR-1; k++ )
		{
			for ( s = gridZ-PML; s < gridZ-1 ; s++)
			{
				curlHr      = ( Hr[k][s] - Hr[k][s-1] )/dz;
				curlHz      = ( Hz[k+1][s] - Hz[k][s] )/dr;
				double coef = dt/(EPSNOT*epsphi[k][s]);
				Ephi[k][s]  = Ephi[k][s] + coef*( ( curlHr + QHzr[k][s] ) - ( curlHz  + QHrz[k][s] ) );
			}
		}

		for ( k = gridR-PML; k < gridR-1; k++ )
		{
			for (s = PML; s < gridZ-PML ; s++)
			{
				curlHr      = ( Hr[k][s] - Hr[k][s-1] )/dz;
				curlHz      = ( Hz[k+1][s] - Hz[k][s] )/dr;
				double coef = dt/(EPSNOT*epsphi[k][s]);
				Ephi[k][s]  = Ephi[k][s] + coef*( ( curlHr + QHzr[k][s] ) - ( curlHz  + QHrz[k][s] ) );
			}
		}
		/*

			for ( k = 0 ; k < gridR-PML; k++ )
			{
				for ( s = PML; s < gridZ - PML; s++ )
				{
					curlHr      = ( Hr[k][s] - Hr[k][s-1] )/dz;
					curlHz      = ( Hz[k+1][s] - Hz[k][s] )/dr;
					Dphi[k][s]  = Dphi[k][s] + dt*(curlHr - curlHz);
					Ephi[k][s]  = ( Dphi[k][s] - PLphi[k][s] )/( EPSNOT*epsphi[k][s] );
					PLphi[k][s] = 0;
				}
			}

			*/
		# pragma endregion //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
		#pragma region //%%%%%%%%% E-FIELD Boundary Update %%%%%%%%%//
		for (k = 0 ; k < gridR; k++)
		{
			Er[k][0]         = 0;
			Er[k][gridZ-1]   = 0;
			Ephi[k][0]       = 0;
			Ephi[k][gridZ-1] = 0;
		}
		for (k = 0; k < gridZ; k++)
			Ephi[gridR-1][k] = 0;
			
		for (k = 0; k < gridZ-1; k++)
			Ez[gridR-1][k] = 0;
		#pragma endregion //%%%%%%% END E-FIELD Boundary Update %%%%%%%//

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
			count = count + 1;
			double** EzWrite = new double* [gridR];
			for (k = 0; k<gridR; k++)
				EzWrite[k] = new double [gridZ];

			for (k = 0; k<gridR; k++)
			{
				for (s = 0; s< gridZ-1; s++)
					EzWrite[k][s] = Ez[k][s];
			}


			ofstream data(fileName.c_str(), ios::out | ios::binary | ios::app);    // Open file to store data
		    //os.write(reinterpret_cast<const char*> (s), sizeof(double) * ns);
			for (k=0; k<gridR; k++)
			{
				for (s = 0; s<gridZ; s++)
				data.write( (char*) &Er[k][s], sizeof Er[k][s]);
			}
			for (k=0; k<gridR; k++)
			{
				for (s = 0; s<gridZ; s++)
				data.write( (char*) &EzWrite[k][s], sizeof EzWrite[k][s]);
			}
				for (k=0; k<gridR; k++)
			{
				for (s = 0; s<gridZ; s++)
				data.write( (char*) &Ephi[k][s], sizeof Ephi[k][s]);
			}

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
		if (T == 1 && maxIter > 99)
		cout<< "Calculating step " << (double) T << " to " << T+98 << " out of " << maxIter << endl;
		if (T == 1 && maxIter < 99)
		cout<< "Calculating step " << (double) T << " to " << maxIter << " out of " << maxIter << endl;

		if ( T%100 == 0 && T < (maxIter - 100) )
		// Display info
		{
			double end = clock();
			double time = end-start;
			cout << "                Current simulation time = " <<  time/CLOCKS_PER_SEC <<  " sec" <<endl;
			cout << "" <<endl;
			cout<< "Calculating step " << (double) T << " to " << T+99 << " out of " << maxIter << endl;
		}

		if ( T%100 == 0 && T < (maxIter)  && T > (maxIter - 100) )
		// Display info
		{
			double end = clock();
			double time = end-start;
			cout << "                Current simulation time = " <<  time/CLOCKS_PER_SEC <<  " sec" <<endl;
			cout<< "Calculating step " << (double) T << " to " << maxIter << " out of " << maxIter << endl;
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
	cin.get();
	return(0);

	#pragma endregion //%%%%%%%%%%%%% END DISPLAY OUTPUT %%%%%%%%%%%%%%//

}//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


