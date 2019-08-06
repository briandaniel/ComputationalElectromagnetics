
#include "3D_FDTD.h"
// Initializes classes defined in model.h

Domain::Domain()
{

	Config settings("settings.txt");
	
	settings.get("Domain/xGrid", xGrid);
	settings.get("Domain/yGrid", yGrid);
	settings.get("Domain/zGrid", zGrid);

	settings.get("Domain/xLength", xLength);
	settings.get("Domain/yLength", yLength);
	settings.get("Domain/zLength", zLength);

	settings.get("Domain/maxTime", maxTime);
	settings.get("Domain/courant", courant);

	dx = xLength/xGrid;
	dy = yLength/yGrid;
	dz = zLength/zGrid;

	dt = ( courant ) / ( CNOT * sqrt( 1/pow(dx,2.0) + 1/pow(dy,2.0) + 1/pow(dz,2.0) ) );

	tSteps = (int) ( maxTime/dt );

}

Efield::Efield( int xIn, int yIn, int zIn, double dtIn )
{
	xGrid = xIn; 
	yGrid = yIn;
	zGrid = zIn;
	dt = dtIn;

	// Field is constructed with the standard order x,y,z;
	field = new double ** [xGrid];  

	for (int i = 0; i < xGrid; i++)
		field[i] = new double * [yGrid]; 
	
	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
			field[i][j] = new double [zGrid]; 
	}


	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			for (int k = 0; k < zGrid; k++)
			{
				field[i][j][k] = 0.0;
			}
		}
	}


	// sigma initialized to zero, eps to one;
	sigma	= new double ** [xGrid];  
	eps		= new double ** [xGrid];  

	for (int i = 0; i < xGrid; i++)
	{
		sigma[i]	= new double * [yGrid]; 
		eps[i]		= new double * [yGrid]; 
	}

	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			sigma[i][j] = new double [zGrid]; 
			eps[i][j]	= new double [zGrid]; 
		}
	}


	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			for (int k = 0; k < zGrid; k++)
			{
				sigma[i][j][k]	= 0.0;
				eps[i][j][k]	= 1.0;
			}
		}
	}

	// Define constants
	c1 = new double ** [xGrid];  
	c2 = new double ** [xGrid];  

	for (int i = 0; i < xGrid; i++)
	{
		c1[i] = new double * [yGrid]; 
		c2[i] = new double * [yGrid]; 
	}

	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			c1[i][j] = new double [zGrid]; 
			c2[i][j] = new double [zGrid]; 
		}
	}


	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			for (int k = 0; k < zGrid; k++)
			{
				c1[i][j][k]	= ( 1.0 - (sigma[i][j][k]*dt) / (2.0*EPSNOT*eps[i][j][k]) ) / ( 1.0 + (sigma[i][j][k]*dt) / (2.0*EPSNOT*eps[i][j][k]) );
				c2[i][j][k]	= ( dt / (EPSNOT*eps[i][j][k]) ) / ( 1.0 + (sigma[i][j][k]*dt) / (2.0*EPSNOT*eps[i][j][k]) ); 
			}
		}
	}

}

Hfield::Hfield( int xIn, int yIn, int zIn, double dtIn )
{

	xGrid = xIn; 
	yGrid = yIn;
	zGrid = zIn;
	dt = dtIn;

	// Field is constructed with the standard order x,y,z;
	field = new double ** [xGrid];  

	for (int i = 0; i < xGrid; i++)
		field[i] = new double * [yGrid]; 
	
	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
			field[i][j] = new double [zGrid]; 
	}


	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			for (int k = 0; k < zGrid; k++)
			{
				field[i][j][k] = 0.0;
			}
		}
	}


	// sigma initialized to zero, eps to one;
	sigma	= new double ** [xGrid];  
	mu		= new double ** [xGrid];  

	for (int i = 0; i < xGrid; i++)
	{
		sigma[i]	= new double * [yGrid]; 
		mu[i]		= new double * [yGrid]; 
	}

	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			sigma[i][j] = new double [zGrid]; 
			mu[i][j]	= new double [zGrid]; 
		}
	}


	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			for (int k = 0; k < zGrid; k++)
			{
				sigma[i][j][k]	= 0.0;
				mu[i][j][k]	= 1.0;
			}
		}
	}

	// Define constants
	d1 = new double ** [xGrid];  
	d2 = new double ** [xGrid];  

	for (int i = 0; i < xGrid; i++)
	{
		d1[i] = new double * [yGrid]; 
		d2[i] = new double * [yGrid]; 
	}

	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			d1[i][j] = new double [zGrid]; 
			d2[i][j] = new double [zGrid]; 
		}
	}


	for (int i = 0; i < xGrid; i++)
	{
		for (int j = 0; j < yGrid; j++)
		{
			for (int k = 0; k < zGrid; k++)
			{
				d1[i][j][k]	= ( 1.0 - (sigma[i][j][k]*dt) / (2.0*MUNOT*mu[i][j][k]) ) / ( 1.0 + (sigma[i][j][k]*dt) / (2.0*MUNOT*mu[i][j][k]) );
				d2[i][j][k]	= ( dt / (MUNOT*mu[i][j][k]) ) / ( 1.0 + (sigma[i][j][k]*dt) / (2.0*MUNOT*mu[i][j][k]) ); 
			}
		}
	}


}

Source::Source()
{
	
	Config settings("settings.txt");
	settings.get("SourceCondition/xPosition", xPos);
	settings.get("SourceCondition/yPosition", yPos);
	settings.get("SourceCondition/zPosition", zPos);

	settings.get("SourceCondition/pulse", sourceType);
	settings.get("SourceCondition/frequency", freq);
	settings.get("SourceCondition/amplitude", amplitude);
	settings.get("SourceCondition/width", width);
	settings.get("SourceCondition/fLength", fLength);
	settings.get("SourceCondition/waist", waist);
	settings.get("SourceCondition/zo", zo);
	settings.get("SourceCondition/vortexN", vortexN);

	
}