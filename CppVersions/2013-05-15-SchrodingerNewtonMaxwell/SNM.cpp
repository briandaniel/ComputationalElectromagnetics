#include "SNM.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void main () {

	cout << "Initializing domain..." << endl << endl;

	// Initialize class variables
	Domain dom;
	Newton ions( dom );
	Schrodinger elc ( dom );
	Maxwell field (dom);
	Energy eVals;
	Manager theBoss;
	theBoss.startClock();
	theBoss.finalData( dom, elc, ions, field);

	//Give initial sech shape
	sech(dom, elc);
	
	// Orthonormalize initial shape
	gramSchmidt(dom, elc);
	
	// Normalize initial shape
	sqNormalize(elc.PSI, dom.N, dom.Npsi, dom.dz);
	
//  | ------------- Calculate ground state -------------|
	double* initDis = new double [dom.N-1];
	for (int k = 0; k < dom.N; k++)
		initDis[k] = ions.init;

	for (int s = 0; s < ions.Nion; s++ )
		ions.Rion[s] = dom.surfPos + ions.init*s;

		
	cout << "Imaginary time step is: " << dom.i_dt*1e18 << " [as]" << endl;
	if (dom.groundIons == 1)
	{
		std::cout << "Calculating ground state of ions using simplex search... " << endl ;
		simplex( initDis, ions.variation, 0.0, ions.yTol, dom, elc, ions, eVals, ions.maxEvals );
		ground( dom, elc, ions, 1e-100 );	
		ions.groundSave();
	}
	else
	{
		ions.groundRead();

		std::cout << std::fixed;

		cout << "Ion positions (from ground file): " << endl;
		for (int k = 0; k < dom.Nion; k++)
			cout << ions.Rion[k]*1e9 << " [nm],     ";

		cout << endl << endl;
		std::cout << "Calculating ground state of electrons only... " << endl ;
		ground( dom, elc, ions, 1e-100 );	
		eVals.matterEnergy(dom, elc, ions);
	}
	
//  | ------------- Real Time Evolution -------------|
	theBoss.clearBins();
	cout << " "<< endl << "Calculating real time evolution ..." << endl;
	cout << " "<< endl << "Real time step is: " << dom.r_dt*1e18 << "[as]" << endl;

	for (int r = 0; r < dom.r_tSteps; r++)
	{

		// 1. Ion position update
		if ( ions.evolve == 1 )
			ionPosition( ions.Rion, ions.Vion, dom.r_dt, dom.Nion );
			
		// 2. PSI update
		if ( elc.evolve == 1 )
			psiStep( dom, elc, ions, field);
		
		// 3. Maxwell update
		if ( field.evolve == 1)
		{
			field.FDTDstep( dom, elc, r);
			interpolateA( dom, elc, field, r);
		}

		// 4. Ion velocity update
		if ( ions.evolve == 1 )
			ionVelocity( dom, elc, ions );

		if ( elc.evolve == 1 )
			eVals.matterEnergy( dom, elc, ions );
		double Etotal = eVals.totalEnergy();
		
		printData( r, dom.r_tSteps,  Etotal, ions.Rion, dom.Nion, eVals, elc.rho, dom.N );
		theBoss.writeVideo(r, dom, elc, ions, field);
		
		

	}

	theBoss.finalData( dom, elc, ions, field);
	theBoss.endClock();
	
}//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


