#ifndef MODEL_H
#define MODEL_H
#include "3D_FDTD.h"

// Generates domain
class Domain
{

  public:
	int	xGrid, yGrid, zGrid;

	double xLength, yLength, zLength;

	double dx, dy, dz;

	int tSteps;

	double dt, maxTime, courant;

	Domain(); // Standard Constructor

};

// Generates Efields
class Efield
{

  public:
	  int xGrid, yGrid, zGrid;
	  double dt;
	  double *** field;
	  double *** sigma;
	  double *** eps;
	  double *** c1;
	  double *** c2;

	  Efield(int xIn, int yIn, int zIn, double dtIn); // Standard Constructor
};

// Generates Hfields
class Hfield
{

   public:
	  int xGrid, yGrid, zGrid;
	  double dt;
	  double *** field;
	  double *** sigma;
	  double *** mu;
	  double *** d1;
	  double *** d2;

	  Hfield(int xIn, int yIn, int zIn, double dtIn); // Standard Constructor
	
};

// Generates source conditions
class Source
{

   public:
	  int xPos, yPos, zPos;
	  string sourceType;
	  double freq;
	  double amplitude;
	  double width;
	  double waist;
	  double fLength;
	  double zo;
	  double vortexN;

	  Source();
	  void gaussian( Domain dom, Efield Ex, Efield Ey, Efield Ez, int r );
	  void gaussBeam( Domain dom, Efield Ex, Efield Ey, Efield Ez, int r );
	  void vortex( Domain dom, Efield Ex, Efield Ey, Efield Ez, int r );
};
#endif

