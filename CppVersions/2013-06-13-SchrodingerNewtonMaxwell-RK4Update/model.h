// model.h

#ifndef MODEL_H
#define MODEL_H

// Generates domain from input file
class Domain
{

  public:
	int	N, Nion, Npsi, r_tSteps, i_tSteps,  EMpoints;
	double zMin, zMax, r_tLength, i_tLength;
	double *zDomPsi;
	double *zDomH, *zDomE;
	int NHy, NEx;
	int groundIons;
	int gVid;
	int autoSurf;
	double surfPos, dz;
	double i_dt_min, i_dt_max, r_dt;
	double fieldCoef;
	Domain(); // Standard Constructor

};

// Generates standard set of variable for ions
class Newton
{

  public:
	double* Rion; 
	double* Vion; 
	double* Mion;
	double* Zion;
	double* Fpsi;
	double* Fion;

	double surfPos;
	int Nion;
	int evolve;
	double c1;
	double c2;
	double c3;

	int maxEvals;
	double init, variation, yTol; 

	Newton( Domain dom ); // Standard constructor
};

// Generates standard set of variable for electrons
class Schrodinger
{

  public:
	CX ** PSI;
	CX ** k1, **k2, **k3, **k4;
	fftw_complex **PSItemp;
	CX *kvec1, *kvec2;
	CX* pmlCoef;
	double* kvec;
	double* rho;
	double* Vxc;
	double* Exc;
	double* coul;
	double* pCoul;
	double* zDomPsi;
	double epsFactor;
	double coreEffect;
	double coreRadius;
	double* A;
	double* Ax;
	double* Vx;
	double* VH;
	double mu;
	double VxcConst;
	double imToli;
	double imTole;
	double* coulPer;
	double* rhoPer;
	double dz;
	int Npsi;
	int N;
	int No;
	int NoPer;
	int pExt;
	int perSize;
	double gamma;
	int evolve;
	double eeCoef;
	string VxcSetting;
	string boundaryCond;
	string pseudo;
	string propagator;

	int PSIpml;
	double PSIPMLsigMax, PSIPMLgamma, PSIsigM;
	CX  b1, b2, b3, b4, b5, b6, b7;

	Schrodinger( Domain dom); // Standard constructor
	void VxcCalc();
	void ExcCalc();
	void VHCalc();

};

// Generates standard set of variable for electrical field solved via FDTD
class Maxwell
{
  private:
	void source( int R );

  public:
	int evolve;
	int EMpoints, NEx, NHy;
	int PML;
	int sourcePos;
	int addlpoints;
	int psiStart;
	int interpN;

	string pulseSetting;
	double maxT, dt, width, pAmp, omega;

	// PML ARRAYS
	double* bx;
	double* by;
	double* gx;
	double* gy;
	double* kx;
	double* ky;
	double* QEx;
	double* QHy;

	double dz;
	double* fieldDom;

	double* Ex, * Hy, * A, *Aold;

	double a1, a2, a3, a4, a5;
	Maxwell(Domain dom);
	void FDTDstep( Domain dom, Schrodinger elc, int r);

};

// Energy variables and functions (func are in modelFunctions.cpp)
class Energy
{
  public:
	double Ue, Uion, Ubond, Uxc, Te, Tion, Etot;
	Energy();
	void matterEnergy ( Domain dom, Schrodinger elc, Newton ions);

};

class Manager
{

	public:
		string psiBin;
		string rhoBin;
		string ionBin;
		string ionVelBin;
		string fdtdBin;
		string ionForceBin;
		string psiForceBin;
		string aBin;
		string hamBin;
		string potentialBin;
		string ionPosBin;
		double startTime;
		
		string saveFileElc;
		string saveFileIon;
		string importFileIon;
		string importFileElc;
		string readFile;
	
		double** Ham;

	int  frames, outputParameters, videoValue;

	//[Output]
	Manager(Domain dom);
	void Hamout( Domain dom, Schrodinger elc, Newton ions, Maxwell field );
	void HamoutUnit( Domain dom, Schrodinger elc, Newton ions, Maxwell field );
	void clearBins();
	void writeVideo(int r, Domain dom, Schrodinger elc, Newton ions, Maxwell field);
	void finalData( Domain dom, Schrodinger elc, Newton ions, Maxwell field );
	void writeForces(Domain dom, Newton ions);
	void startClock();
	void endClock();

};

#endif

