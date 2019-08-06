// utility.cpp
#include "SN_3D.h"

Manager::Manager(Domain dom)
{
Config settings("settings.txt");

	settings.get("Other/frames", frames);

	Nt_re = dom.Nt_re;
	Nt_im =  dom.Nt_im;
	N3 = dom.N3;

	ofstream data("vid.bin", ios::out | ios::binary); 
	data.close();

	ofstream data2("vidCX.bin", ios::out | ios::binary | ios::app); 
	data2.close();
}

void Manager::vidWrite( int r, double* rho, CX* U, CX dt ){

	int Nt;

	if ( real(dt) == 0)
		Nt = Nt_im;
	else
		Nt = Nt_re;

	if ( mod(r, Nt/frames) == 0 )
	{
		
		ofstream data("vid.bin", ios::out | ios::binary | ios::app); 
		ofstream data2("vidCX.bin", ios::out | ios::binary | ios::app); 

		for (int k = 0; k < N3; k++)
			data.write((char *) & rho[k], sizeof rho[k]);

		
		for (int k = 0; k < N3; k++)
			data2.write((char *) & U[k], sizeof U[k]);

		cout << r << endl;
		
		cout << U[0] << endl;

		data2.close();
		data.close();
	}


}

void Manager::initalData( Domain dom, Schrodinger elc, Derivative D2 )
{
		
	ofstream doubleData("doubleData.bin", ios::out | ios::binary); 
	ofstream intData("intData.bin", ios::out | ios::binary); 
	ofstream data("data.bin", ios::out | ios::binary); 
	for (int k = 0; k < dom.N3; k++)
		data.write((char *) & elc.rho[k], sizeof elc.rho[k]);

	for (int k = 0; k < dom.N3; k++)
		data.write((char *) & elc.Vext[k], sizeof elc.Vext[k]);
	
	intData.write((char *) & dom.N3, sizeof dom.N3);
	intData.write((char *) & dom.Nx, sizeof dom.Nx);
	intData.write((char *) & dom.Ny, sizeof dom.Ny);
	intData.write((char *) & dom.Nz, sizeof dom.Nz);
	intData.write((char *) & frames, sizeof frames);

	for (int k = 0; k < D2.M; k++)
	{
		intData.write((char *) & D2.Di[0][k], sizeof D2.Di[0][k]);
		intData.write((char *) & D2.Di[1][k], sizeof D2.Di[1][k]);
		doubleData.write((char *) & D2.D[k], sizeof D2.D[k]);
	}

	doubleData.close();
	intData.close();
	data.close();
}

void Manager::startTime()
{
	beginT = clock();
	cout << "Starting simulation ... " << endl;
}

void Manager::midTime()
{
	
	double tempT = clock();
	elapTime = ((double)(tempT - beginT))/CLOCKS_PER_SEC;
	cout << "Current simulation time is " << elapTime << "[s]" << endl;
}


void Manager::endTime()
{
	endT = clock();

	elapTime = ((double)(endT - beginT))/CLOCKS_PER_SEC;

	cout << "Program finished." << endl;
	cout << "Total simulation time is " << elapTime << "[s]" << endl;
	cout << "Press enter to exit." << endl;
	cin.get();

}