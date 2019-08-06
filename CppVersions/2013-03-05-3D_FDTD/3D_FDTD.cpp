#include "3D_FDTD.h"



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void main () {


	///////////////////////////////////

	double startTime = clock();
	
	// Initialize domain
	cout << "Initializing domain..." << endl;
	Domain dom;
	Source initial;

	// Initialize Efield variables
	Efield Ex(dom.xGrid     , dom.yGrid + 1 , dom.zGrid + 1 , dom.dt);
	Efield Ey(dom.xGrid + 1 , dom.yGrid     , dom.zGrid + 1 , dom.dt);
	Efield Ez(dom.xGrid + 1 , dom.yGrid + 1 , dom.zGrid     , dom.dt);

	// Initialize Hfield variables
	Hfield Hx(dom.xGrid - 1 , dom.yGrid     , dom.zGrid     , dom.dt);
	Hfield Hy(dom.xGrid     , dom.yGrid - 1 , dom.zGrid     , dom.dt);
	Hfield Hz(dom.xGrid     , dom.yGrid     , dom.zGrid - 1 , dom.dt);

	// Main solution loop
	cout << "Calculating main solution loop with " << dom.tSteps << " steps... " << endl;
	ofstream data("video.bin", ios::out | ios::binary);
	data.close();

	ofstream dataNew("video.bin", ios::out | ios::binary | ios::app);
	for (int r = 0; r < dom.tSteps; r++)
	{

		// Parallel processing commands
		#ifdef USE_OPENMP	
			omp_set_num_threads(OMP_NUM_THREADS);
			#pragma omp parallel
		#endif

		HxUpdate( dom, Hx, Ey, Ez );
		HyUpdate( dom, Hy, Ez, Ex );
		HzUpdate( dom, Hz, Ex, Ey );

		ExUpdate( dom, Ex, Hz, Hy );
		EyUpdate( dom, Ey, Hx, Hz );
		EzUpdate( dom, Ez, Hy, Hx );

		initial.vortex (dom, Ez, Ex, Ey, r );

	if ( mod(r,dom.tSteps/10) == 0 )
		cout << "Steps: " << r + 1 << " - " << r + dom.tSteps/10 << endl;

	}
	writeData( Ex, "Ex.bin" );
	writeData( Ey, "Ey.bin" );

	double endTime = clock();
	double programTime = (endTime - startTime)/CLOCKS_PER_SEC ;
	cout << endl << "Simulation complete!" << endl;
	cout << "Total Simulation Time = " <<  programTime  <<  " sec" <<endl << endl;
	cout << "Press enter to exit cmd window" << endl;
	std::cin.get();

}//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


