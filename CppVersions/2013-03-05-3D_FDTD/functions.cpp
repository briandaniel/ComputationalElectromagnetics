#include "3D_FDTD.h"


/* FUNCTION LIST 

1. Hx Update
2. Hy Update
3. Hz Update

4. Ex Update
5. Ey Update
6. Ez Update

7. source

*/

// Updates Hx one time step
void HxUpdate(Domain dom, Hfield Hx, Efield Ey, Efield Ez)
{
	for (int i = 0; i < Hx.xGrid; i++)
	{
		for (int j = 0; j < Hx.yGrid; j++)
		{
			for (int k = 0; k < Hx.zGrid; k++)
			{
				Hx.field[i][j][k] = Hx.d1[i][j][k] * Hx.field[i][j][k]
				     + Hx.d2[i][j][k] * ( ( Ey.field[i+1][j][k+1] - Ey.field[i+1][j][k] ) / dom.dz
					 -					  ( Ez.field[i+1][j+1][k] - Ez.field[i+1][j][k] ) / dom.dy ) ;
			}
		}
	}
}

// Updates Hy one time step
void HyUpdate(Domain dom, Hfield Hy, Efield Ez, Efield Ex)
{
	for (int i = 0; i < Hy.xGrid; i++)
	{
		for (int j = 0; j < Hy.yGrid; j++)
		{
			for (int k = 0; k < Hy.zGrid; k++)
			{
				Hy.field[i][j][k] = Hy.d1[i][j][k] * Hy.field[i][j][k]
				     + Hy.d2[i][j][k] * ( ( Ez.field[i+1][j+1][k] - Ez.field[i][j+1][k] ) / dom.dx
					 -					  ( Ex.field[i][j+1][k+1] - Ex.field[i][j+1][k] ) / dom.dz ) ;
			}
		}
	}
}

// Updates Hz one time step
void HzUpdate(Domain dom, Hfield Hz, Efield Ex, Efield Ey)
{
	for (int i = 0; i < Hz.xGrid; i++)
	{
		for (int j = 0; j < Hz.yGrid; j++)
		{
			for (int k = 0; k < Hz.zGrid; k++)
			{
				Hz.field[i][j][k] = Hz.d1[i][j][k] * Hz.field[i][j][k]
				     + Hz.d2[i][j][k] * ( ( Ex.field[i][j+1][k+1] - Ex.field[i][j][k+1] ) / dom.dy
					 -					  ( Ey.field[i+1][j][k+1] - Ey.field[i][j][k+1] ) / dom.dx ) ;
			}
		}
	}
}

// Updates Ex one time step
void ExUpdate(Domain dom, Efield Ex, Hfield Hz, Hfield Hy)
{
	for (int i = 0; i < Ex.xGrid; i++)
	{
		for (int j = 1; j < ( Ex.yGrid - 1 ) ; j++)
		{
			for (int k = 1; k < ( Ex.zGrid - 1 ); k++)
			{
				Ex.field[i][j][k] = Ex.c1[i][j][k] * Ex.field[i][j][k]
				     + Ex.c2[i][j][k] * ( ( Hz.field[i][j][k-1] - Hz.field[i][j-1][k-1] ) / dom.dy
					 -					  ( Hy.field[i][j-1][k] - Hy.field[i][j-1][k-1] ) / dom.dz ) ;
			}
		}
	}
}

// Updates Ey one time step
void EyUpdate(Domain dom, Efield Ey, Hfield Hx, Hfield Hz)
{
	for (int i = 1; i < ( Ey.xGrid - 1 ); i++)
	{
		for (int j = 0; j < Ey.yGrid ; j++)
		{
			for (int k = 1; k < ( Ey.zGrid - 1 ); k++)
			{
				Ey.field[i][j][k] = Ey.c1[i][j][k] * Ey.field[i][j][k]
				     + Ey.c2[i][j][k] * ( ( Hx.field[i-1][j][k] - Hx.field[i-1][j][k-1] ) / dom.dz
					 -					  ( Hz.field[i][j][k-1] - Hz.field[i-1][j][k-1] ) / dom.dx ) ;
			}
		}
	}
}

// Updates Ez one time step
void EzUpdate(Domain dom, Efield Ez, Hfield Hy, Hfield Hx)
{
	for (int i = 1; i < ( Ez.xGrid - 1 ); i++)
	{
		for (int j = 1; j < ( Ez.yGrid -1 ) ; j++)
		{
			for (int k = 0; k < Ez.zGrid ; k++)
			{
				Ez.field[i][j][k] = Ez.c1[i][j][k] * Ez.field[i][j][k]
				     + Ez.c2[i][j][k] * ( ( Hy.field[i][j-1][k] - Hy.field[i-1][j-1][k] ) / dom.dx
					 -					  ( Hx.field[i-1][j][k] - Hx.field[i-1][j-1][k] ) / dom.dy ) ;
			}
		}
	}
}

// Source definition
void Source::gaussian ( Domain dom, Efield Ez, Efield Ex, Efield Ey, int r )
{
	double T = r*dom.dt;
	double source;
	double maxT = .5*width;
	double omega = 2*PI*freq;

	source = exp( -1.0*pow(maxT-T,2.0) / pow(width/5.0, 2.0) ) * amplitude  * sin(omega*T);

	int x = xPos;
	int y = yPos;
	int z = zPos;

	Ex.field[x][y][z] = source;
	Ez.field[x][y][z] = source;
	Ey.field[x][y][z] = source;
}

void Source::gaussBeam ( Domain dom, Efield Ez, Efield Ex, Efield Ey, int r )
{
	double T, timeComp, maxT, omega, lambda, zr, kNum, Rcurve, zNot, waistz; 

	T = r*dom.dt;
	maxT = .5*width;
	omega = 2*PI*freq;
	lambda = CNOT/freq;
	zr = PI*pow(waist,2.0)/lambda;
	kNum = 2*PI/lambda;
	zNot = fLength;
	Rcurve = zNot* ( 1 + pow(zr/zNot, 2.0) );
	waistz = waist*sqrt(1 + pow( zNot/zr, 2.0 ) );

	//timeComp = exp( -1.0*pow(maxT-T,2.0) / pow(width/5.0, 2.0) ) * amplitude  * sin(omega*T);
	timeComp = amplitude;

	int x = xPos;
	int y = yPos;
	int z = zPos;
	
	for (int i = 0; i < dom.xGrid; i++)
	{
		for (int j = 0; j < dom.yGrid; j++)
		{
					x = i - xPos;
					y = j - yPos;
					double r2 = pow(x*dom.dx, 2.0) + pow(y*dom.dx, 2.0) ;

					CX temp =  exp( (CX) - r2/ pow(waistz,2.0) - CX( 0.0, omega*r*dom.dt ) - CX( 0.0, kNum*r2/(2*Rcurve) )  - CX( 0.0, kNum*zo) );

					Ex.field[i][j][z] = ( waist/waistz ) * imag(temp);
					Ey.field[i][j][z] = ( waist/waistz ) * imag(temp);
		}
	}

}



void Source::vortex ( Domain dom, Efield Ez, Efield Ex, Efield Ey, int r )
{
	double T, timeComp, maxT, omega, lambda, zr, kNum, Rcurve, zNot, waistz; 

	T = r*dom.dt;
	maxT = .5*width;
	omega = 2*PI*freq;
	lambda = CNOT/freq;
	zr = PI*pow(waist,2.0)/lambda;
	kNum = 2*PI/lambda;
	zNot = fLength;
	Rcurve = zNot* ( 1 + pow(zr/zNot, 2.0) );
	waistz = waist*sqrt(1 + pow( zNot/zr, 2.0 ) );

	//timeComp = exp( -1.0*pow(maxT-T,2.0) / pow(width/5.0, 2.0) ) * amplitude  * sin(omega*T);
	timeComp = amplitude;

	int x = xPos;
	int y = yPos;
	int z = zPos;
	
	for (int i = 0; i < dom.xGrid; i++)
	{
		for (int j = 0; j < dom.yGrid; j++)
		{
					x = i - xPos;
					y = j - yPos;
					double r2 = pow(x*dom.dx, 2.0) + pow(y*dom.dx, 2.0) ;


					double phi;
					phi = atan2 ( (y*dom.dx) , (x*dom.dx) );
					CX polarize = exp( CX(0.0, vortexN*phi) );
					CX temp = polarize*exp( (CX) - r2/ pow(waistz,2.0) - CX( 0.0, omega*r*dom.dt ) - CX( 0.0, kNum*r2/(2*Rcurve) )  - CX( 0.0, kNum*zo) );

					Ex.field[i][j][z] = ( waist/waistz ) * imag(temp);
					Ey.field[i][j][z] = ( waist/waistz ) * imag(temp);
		}
	}

}