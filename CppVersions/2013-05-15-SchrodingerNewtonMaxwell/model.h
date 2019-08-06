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
	double surfPos, dz;
	double i_dt, r_dt;
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

	double surfPos;
	int Nion;
	int evolve;
	double c1;
	double c2;
	double c3;

	string groundFile;
	string readFile;

	int maxEvals;
	double init, variation, yTol; 

	Newton( Domain dom ); // Standard constructor
	void groundSave(); // save ground state file
	void groundRead();
};

// Generates standard set of variable for electrons
class Schrodinger
{

  public:
	CX ** PSI;
	CX* pmlCoef;
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
	double VxcConst;
	double imTol;
	int Npsi;
	int N;
	int No;
	double gamma;
	int evolve;
	double eeCoef;
	string VxcSetting;
	string boundaryCond;
	string pseudo;

	double pseudoCharge;
	int PSIpml;
	double PSIPMLsigMax, PSIPMLgamma, PSIsigM;
	CX  b1, b2, b3, b4, b5, b6;

	Schrodinger( Domain dom); // Standard constructor
	void VxcCalc();
	void ExcCalc();
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
	double Ue, Uion, Ubond, Uxc, Te, Tion;
	Energy();
	double totalEnergy();
	void matterEnergy ( Domain dom, Schrodinger elc, Newton ions);

};

class Manager
{
	public:
		string psiBin;
		string rhoBin;
		string ionBin;
		string fdtdBin;
		double startTime;

	int  frames, outputParameters, videoValue;

	//[Output]
	Manager();
	void clearBins();
	void writeVideo(int r, Domain dom, Schrodinger elc, Newton ions, Maxwell field);
	void finalData( Domain dom, Schrodinger elc, Newton ions, Maxwell field );
	void startClock();
	void endClock();

};

#endif

