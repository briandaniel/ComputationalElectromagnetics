#include "SNM.h"

void Manager::clearBins()
{
	ofstream data(psiBin, ios::out | ios::binary); 
	data.close();

	ofstream data1(rhoBin, ios::out | ios::binary); 
	data1.close();

	ofstream data2(ionBin, ios::out | ios::binary); 
	data2.close();

	ofstream data3(fdtdBin, ios::out | ios::binary); 
	data3.close();
}

void Manager::writeVideo(int r, Domain dom, Schrodinger elc, Newton ions, Maxwell field)
{

	ofstream data(psiBin, ios::out | ios::binary | ios::app); 
	ofstream data1(rhoBin, ios::out | ios::binary | ios::app); 
	ofstream data2(ionBin, ios::out | ios::binary | ios::app); 
	ofstream data3(fdtdBin, ios::out | ios::binary | ios::app); 

	if ( mod(r,round(dom.r_tSteps/frames)) == 0  && videoValue == 1)
		{
			for (int s=0; s< dom.Npsi; s++)
			{
				for (int k = 0; k < dom.N; k++)
				{
					double temp = real(elc.PSI[s][k]);
					data.write((char *) & temp, sizeof temp);
				}
			}
			
		for (int k = 0; k < dom.N; k++)
			{
				data1.write((char *) & elc.rho[k], sizeof elc.rho[k]);
			}

			for (int s=0; s< dom.Nion; s++)
			{
				data2.write((char *) & ions.Rion[s], sizeof ions.Rion[s]);
			}

			for (int s=0; s< field.NEx; s++)
			{
				data3.write((char *) & field.Ex[s], sizeof field.Ex[s]);
			}

		}

		data.close();
		data1.close();
		data2.close();
		data3.close();
}

void Manager::startClock()
{
	startTime = clock();
	srand( (unsigned) time(0)); 
}

void Manager::endClock()
{
	double endTime = clock();
	double programTime = (endTime - startTime)/CLOCKS_PER_SEC ;

	cout << "Simulation complete!" << endl;
	int hours = (int) ( programTime/3600);
	int minutes = (int) ( (programTime - hours*3600 )/60);
	programTime = ( programTime - 60* minutes - 3600*hours);
	cout << "Total Simulation Time = " << hours << " : "  << minutes <<  " : "  << (int) programTime << endl;
	cout << "" <<endl;

	cout << "Press enter to exit cmd window" << endl;
	
	std::cin.get();
}

void printData( int r, int steps,  double Etotal, double* Rion, int Nion, Energy eVals, double* rho, int N )
{
	
		if (mod(r,200) == 0)
		{
			std::cout << std::fixed;
			cout << endl << "Step: " << r << "/ " << steps 
				<< ",  Ion positions: " << endl;
			for (int k = 0; k < Nion; k++)
				cout << Rion[k]*1e9 << " [nm],   ";
			cout << endl << "Utot = " << Etotal << "  [eV]";
			cout << ",          Int(rho) = " << trapz(rho, N, 1) << endl;

			cout << 
				"Te = " << eVals.Te << "[eV],  " <<
				"Tion = " << eVals.Tion << "[eV],  " <<
				"Uee = " << eVals.Ue << "[eV],  " <<
				"Uxc = " << eVals.Uxc << "[eV],  " <<
				"Uii = " << eVals.Uion << "[eV],  " <<
				"Uie = " << eVals.Ubond << "[eV],  " <<
				"Utot = " << Etotal << "[eV],  " <<
				endl;
		}
}

void Manager::finalData( Domain dom, Schrodinger elc, Newton ions, Maxwell field )
{
	cout << endl << "Writing data..." << endl;
	
	ofstream data("PSI.bin", ios::out | ios::binary); // Open file to store data
	double doubleN = (double) dom.N;
	double doubleNpsi = (double) dom.Npsi;
	double doubleNEx = (double) field.NEx;
	double doubleEMpoints = (double) ( field.interpN );
	double psiStart = (double) field.psiStart;
	double doubleDT = (double) dom.r_tSteps/ (double) frames * dom.r_dt;
	double doubleFrames = (double)frames;

	data.write((char *) &doubleN, sizeof doubleN);
	data.write((char *) &doubleNpsi, sizeof doubleNpsi);
	data.write((char *) &doubleNEx, sizeof doubleNEx);
	data.write((char *) &doubleEMpoints, sizeof doubleEMpoints);
	data.write((char *) &psiStart, sizeof psiStart);
	data.write((char *) &doubleDT, sizeof doubleDT);
	data.write((char *) &doubleFrames, sizeof doubleFrames);
	data.write((char *) &field.pAmp, sizeof field.pAmp);
	data.write((char *) & dom.zMin, sizeof dom.zMin);
	data.write((char *) & dom.zMax, sizeof dom.zMax);
	
	for (int s = 0; s < elc.N; s++)
		data.write((char *) & elc.pCoul[s], sizeof elc.pCoul[s]);

	for (int s = 0; s < doubleEMpoints; s++)
		data.write((char *) & field.fieldDom[s], sizeof field.fieldDom[s]);


	for (int s=0; s< dom.Npsi; s++)
	{
		for (int k = 0; k < dom.N; k++)
		{
			data.write((char *) &elc.PSI[s][k], sizeof elc.PSI[s][k]);
		}
	}

	ofstream data2("RN.bin", ios::out | ios::binary); // Open file to store data
	for (int s=0; s < dom.Nion; s++)
	{
		data2.write((char *) &ions.Rion[s], sizeof ions.Rion[s]);
	}
	
	ofstream data3("domain.bin", ios::out | ios::binary); // Open file to store data
	for (int k=0; k < dom.N; k++)
	{
		data3.write((char *) &dom.zDomPsi[k], sizeof dom.zDomPsi[k]);
	}
	
	data.close();
	data2.close();
	data3.close();
}