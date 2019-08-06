#include "SN_3D.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
int main () {


	CX dt;

	Domain dom ;
	Manager boss(dom);
	boss.startTime();
	Derivative D2( 2.0, 2.0, dom.boundary,  dom.Nx, dom.Ny, dom.Nz, dom.Lx, dom.Ly, dom.Lz);
	Schrodinger elc(dom);

	boss.initalData(dom, elc, D2);
	boss.midTime();


	// Imaginary time
	dt = CX( 0.0, - dom.dt_im);
		cout << dt << endl;
	RK4 ( elc, dom, D2, boss, dt);


	// Real time
	dt = CX(dom.dt_re, 0.0);
		cout << dt << endl;
	RK4 ( elc, dom, D2, boss, dt);
	
	boss.endTime();
	

	return(0);
	
}//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


