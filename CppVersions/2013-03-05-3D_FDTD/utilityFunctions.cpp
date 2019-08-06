#include "3D_FDTD.h"

void writeData( Efield aEfield, string name )
{
	ofstream data(name, ios::out | ios::binary); // Open file to store data

	double x = aEfield.xGrid;
	double y = aEfield.yGrid;
	double z = aEfield.zGrid;

	data.write((char *) &x, sizeof x);
	data.write((char *) &y, sizeof y);
	data.write((char *) &z, sizeof z);

	for (int i = 0; i < aEfield.xGrid; i++)
	{
		for (int j = 0; j < aEfield.yGrid; j++)
		{
			for (int k = 0; k < aEfield.zGrid; k++)
			{
			data.write((char *) &aEfield.field[i][j][k], sizeof aEfield.field[i][j][k]);
			}
		}
	}

	data.close();
}