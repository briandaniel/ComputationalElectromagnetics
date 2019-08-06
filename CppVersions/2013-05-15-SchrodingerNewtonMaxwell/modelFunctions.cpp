#include "SNM.h"

// Add up total Energy
double Energy::totalEnergy()
{
	double Etot = Te + Tion + Ubond + Ue + Uion + Uxc;
	return Etot;
}

// Calculates energy of the electron/ion models
void Energy::matterEnergy ( Domain dom, Schrodinger elc, Newton ions)
{

	CX* PSItemp = new CX [dom.N];
	double* integrand = new double [dom.N];
	double* Wr = new double [dom.N];
	double* rhoReg = new double [dom.N*2-1];

	//  [Te]  Calculate kinetic energy of electrons
	Te = 0;
	CX TeCX = (0.0, 0.0);
	for (int s = 0; s< dom.Npsi; s++)
	{
		CXderivative2(PSItemp, elc.PSI[s], dom.dz, dom.N);
		for (int k = 0; k<dom.N; k++)
		{
			PSItemp[k] = PSItemp[k]*conj(elc.PSI[s][k]);
		}

		double value = trapz(integrand, dom.N, dom.dz);
		CX integral = CXtrapz(PSItemp, dom.N, dom.dz);
		TeCX = TeCX - ( 2.0* pow(HBAR,2.0)/(2.0*ME) )* integral;
	}
	

	delete [] PSItemp; 

	Te = abs(TeCX);
	
	//  [Tion]  Calculate kinetic energy of ions
	Tion = 0.0;
	for (int s = 0; s < dom.Nion; s++ )
	{
		Tion = Tion + abs( .5*ions.Mion[s]* pow(ions.Vion[s],2.0) );
	}

	//  [Ue]  Calculate potential energy of e-e repulsion
 	fftConv (Wr, elc.rho, elc.coul, dom.N, elc.No);			// Calculate convolution
	density (elc.PSI, rhoReg, dom.Npsi, dom.N);				// Calculate density of electrons
	for (int k = 0; k < dom.N; k++)
	{
		rhoReg[k] = rhoReg[k]*Wr[k];
	}
	double coef = elc.eeCoef * pow(E,2.0) * elc.epsFactor;
	Ue = .5* coef *trapz( rhoReg, dom.N , dom.dz );

	//  [Uxc]  Caculate energy in exchange-correlation potential term
	elc.ExcCalc();
	density (elc.PSI, rhoReg, dom.Npsi, dom.N);		
	for (int k = 0; k < dom.N; k++)
	{
		rhoReg[k] = rhoReg[k]*elc.Exc[k];
	}
	coef = elc.VxcConst;
	Uxc = - coef *trapz( rhoReg, dom.N , dom.dz );

	/*
	if (elc.b6 == 0.0)
	{
		Ue = 0.0;
	}
	*/

	densitydz (elc.PSI, elc.rho, dom.Npsi, dom.N, dom.dz);	// Calculate density of electrons
	//  [Uie]  Calculate ion-e bonding potential energy
	Ubond = 0.0;
	VxCalc(dom, elc, ions);
	double modify;

	for (int k = 0; k<dom.N; k++)
	{
			integrand[k] = elc.rho[k] * elc.Vx[k];	
	}
	Ubond = pow(E,2.0)*elc.epsFactor*trapz(integrand,dom.N,1.0) ;

	//  [Uii]  Calculate ion-ion repulsion

	Uion = 0.0;
	for (int s = 0; s < dom.Nion; s++)
	{
		for (int ss = 0; ss < dom.Nion; ss++)
		{
			if (s != ss)
			{
				double val = ions.Zion[s]*ions.Zion[ss]*pow(E,2.0)*elc.epsFactor;
				double Enm = val/sqrt( pow(elc.gamma,2.0) + pow( (ions.Rion[ss] - ions.Rion[s]) , 2.0 ) );	
				Uion = Uion  + .5*Enm;
			}
		}
	}


	Te = Te*JEV;
	Ue = Ue*JEV;
	Uxc = Uxc*JEV;
	Tion= Tion*JEV;
	Uion = Uion*JEV;
	Ubond = Ubond*JEV;

	delete [] integrand;
	delete [] Wr;
	delete [] rhoReg;
}

// Saves ground positions of ions after runnin simplex search
void Newton::groundSave()
{
		ofstream data(groundFile, ios::out | ios::binary );

		data.write((char *) & Nion, sizeof Nion);

		for (int s=0; s< Nion - 1 ; s++)
		{
			double temp = Rion[s+1] - Rion[s];
			data.write((char *) & temp, sizeof temp);
		}
}

// Reads grounds positions when simplex search is not run
void Newton::groundRead()
{

		if (readFile == "file") {
			ifstream data(groundFile, ios::out | ios::binary );

			data.read((char *) & Nion, sizeof Nion);

			double temp;
			Rion[0] = surfPos;
			for (int s=1; s< Nion; s++)
			{
				data.read((char *) & temp, sizeof temp);
				Rion[s] = Rion[s-1] + temp;
			}
		}else if (readFile == "equalSpacing"){
			Rion[0] = surfPos;
			for (int s=1; s< Nion; s++)
			{
				Rion[s] = Rion[s-1] + init;
			}
		}
}