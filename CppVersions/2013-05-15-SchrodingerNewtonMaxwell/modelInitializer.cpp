#include "SNM.h"

Domain::Domain()
{

	Config settings("settings.txt");
	
	settings.get("Domain/zMin", zMin);
	settings.get("Domain/zMax", zMax);
	settings.get("Domain/Ne", Npsi);
	settings.get("Domain/Nion", Nion);
	settings.get("Domain/Ngrid", N);
	settings.get("Domain/Nem", EMpoints);
	settings.get("Domain/imaginary_tSteps", i_tSteps);
	settings.get("Domain/real_tSteps", r_tSteps);
	settings.get("Domain/real_tLength", r_tLength);
	settings.get("Domain/imaginary_tLength", i_tLength);
	settings.get("Domain/surfPos", surfPos);
	settings.get("Domain/fieldCoef", fieldCoef);
	settings.get("Ground/groundIons", groundIons);
	settings.get("Ground/groundVideo", gVid);
	double stepType;
	settings.get("Domain/stepType", stepType);
	
	dz = (zMax - zMin)/(N-1);

	if (stepType == 1){
		settings.get("Domain/imag_tStep", i_dt);
		settings.get("Domain/real_tStep", r_dt);
	} else if (stepType == 2){
		i_dt = i_tLength/i_tSteps;
		r_dt = r_tLength/r_tSteps;
	}
	

	zDomPsi = new double [N];
	for (int k = 0; k < N; k++)
	{
		zDomPsi[k] = dz*k;
	}

	NHy = N + EMpoints*2;
	NEx = N + EMpoints*2 + 1;

	zDomH = new double [NHy];
	zDomE = new double [NEx];

	for (int k = 0; k < NHy; k++)
		zDomH[k] = dz*(k - EMpoints);

	for (int k = 0; k < NEx; k++)
		zDomE[k] = dz*(k - EMpoints - .5);

}

Newton::Newton( Domain dom)
{

	Config settings("settings.txt");

	double Zn;
	double mass;
	settings.get("Newton/Zn", Zn);
	settings.get("Newton/Mion", mass);
	settings.get("Newton/evolve", evolve);
	settings.get("Ground/groundFile", groundFile);
	settings.get("Ground/initialDistance", init);
	settings.get("Ground/variation", variation);
	settings.get("Ground/yTol", yTol);
	settings.get("Ground/maxEvals", maxEvals);
	settings.get("Ground/readFile", readFile);

	if (evolve == 1)
		cout << "Newton time evolution on" << endl ;
	else
		cout << "Newton time evolution off" << endl ;

	Nion = dom.Nion;
	Vion = new double [Nion];
	Rion = new double [Nion];
	Zion = new double [Nion];
	Mion = new double [Nion];

	for (int s = 0; s<Nion; s++)
	{
		Vion[s] = 0;
		Rion[s] = 0;
		Zion[s] = Zn;
		Mion[s] = MP*Zn*mass;
	}

	// Real constants for Newton's equation	
	double epsFactor = 1/(4*PI*EPSNOT);
	c1 = epsFactor*pow(E,2.0);
	c2 = epsFactor*pow(E,2.0);
	c3 = epsFactor*E;

	surfPos = dom.surfPos;
}

Schrodinger::Schrodinger(Domain dom)
{
	Config settings("settings.txt");

	settings.get("PML/PSIpml", PSIpml);
	settings.get("PML/PSIsigM", PSIsigM);
	settings.get("PML/PSIPMLsigMax", PSIPMLsigMax);
	settings.get("PML/PSIPMLgamma", PSIPMLgamma);
	settings.get("Parameters/gamma", gamma);
	settings.get("Schrodinger/eeCoef", eeCoef);
	settings.get("Schrodinger/evolve", evolve);
	settings.get("Schrodinger/Vxc", VxcSetting);
	settings.get("Schrodinger/VxcConstant", VxcConst);
	settings.get("Schrodinger/boundary", boundaryCond);
	settings.get("Ground/imaginaryTolerance", imTol);
	settings.get("Schrodinger/pseudo", pseudo);
	settings.get("Schrodinger/coreEffect", coreEffect);
	settings.get("Schrodinger/coreRadius", coreRadius);
	settings.get("Schrodinger/charge", pseudoCharge);

	if (evolve == 1)
		cout << "Schrodinger real time evolution on" << endl << endl;
	else
		cout << "Schrodinger real time evolution off" << endl << endl;

	Npsi = dom.Npsi;
	N = dom.N;

	No = (int) floor(dom.N/2);

	PSI = new CX *[Npsi];
	for (int k = 0; k < Npsi; k++)
		PSI[k] = new CX [N];

	for (int s = 0; s < Npsi; s++)
	{
		for (int k = 0; k < N; k++)
			PSI[s][k] = 0.0;
	}

	
	A = new double [N];
	Ax = new double[N];
	Vx = new double[N];
	Vxc = new double[N];
	Exc = new double[N];
	rho = new double [2*N-1];
	zDomPsi = new double[N];
	coul = new double [2*N-1];
	pCoul = new double [2*N-1];
	pmlCoef = new CX [N];

	for (int k = 0; k<N; k++)
	{
		zDomPsi[k] = dom.zDomPsi[k];
		A[k]  = 0.0;
		Vxc[k] = 0.0;
		Exc[k] = 0.0;
		Ax[k] = 0.0;
		Vx[k] = 0.0;
		pmlCoef[k] = 0.0;
	}

	if (boundaryCond == "pml")
	{
		for (int k = 0; k <=  PSIpml; k++)
		{
			double pmlVal = PSIPMLsigMax * pow( (double) k/ ( (double) (PSIpml-1) ), PSIsigM);
			CX pmlGamma = CX( 0.0, pmlVal);
			pmlCoef[PSIpml - k] = 1.0/( 1.0 - pmlGamma ) - 1.0;
			pmlCoef[N - PSIpml - 1 + k] = 1.0/( 1.0 - pmlGamma ) - 1.0;
		}
	}

	double mid = ( zDomPsi[N-1] - zDomPsi[0] )/2 + zDomPsi[0];
	double modify;
	for (int k = 0; k<N; k++)
	{
		//double sg =  exp(- pow( ( zDomPsi[k] - dom.dz*(No-1)),8.0) / pow( (dom.zMax/4),8.0) );
		if (pseudo == "exp")
		{
			double zLoc = zDomPsi[k] - dom.dz*(No-1);
			modify = 1 - exp( -pow (zLoc,2.0) /  pow (coreRadius,2.0) )*coreEffect;
		}else{
		modify = 1;
		}
		coul[k] = 1/sqrt( pow(gamma,2.0) + pow(zDomPsi[k] - dom.dz*(No-1), 2.0) );
		pCoul[k] = modify/sqrt( pow(gamma,2.0) + pow(zDomPsi[k] - dom.dz*(No-1), 2.0) );
	}

	for (int k = N; k< 2*N-1; k++)
	{
		coul[k] = 0.0;
		pCoul[k] = 0.0;
	}
	for (int k = 0; k< 2*N-1; k++)
		rho[k] = 0.0;

	// Constants for Schrodinger equation
	epsFactor = 1/(4*PI*EPSNOT);
	// CPLX(          a              +                  b*i                     );
	b1 = CX(         0.0             ,             HBAR/(2*ME)					);
	b2 = CX(       -E / ME           ,                 0.0						)*dom.fieldCoef;
	b3 = CX(     -E / ( 2.0*ME )     ,                 0.0						)*dom.fieldCoef;
	b4 = CX( pow(E,2.0) / (2*HBAR*ME),                 0.0						)*pow( dom.fieldCoef, 2.0);
	b5 = CX(         0.0             , epsFactor*( - pow(E,2.0) / HBAR )		);
	b6 = CX(         0.0             , epsFactor*( - pow(E,2.0) / HBAR ) *eeCoef);
	

}

Maxwell::Maxwell(Domain dom)
{
	
	double dzApprox = dom.r_dt*CNOT/.5;
	double pointsApprox = ( dom.zMax - dom.zMin )/dzApprox;

	EMpoints = 1000;
	for (int k = (int) pointsApprox; k < (int) pointsApprox*2; k++)
	{
		if (mod(dom.N, k) == 0)
			EMpoints = k;
	}

	if ( EMpoints  == 1000 )
	{
		for (int k = 0; k < 1e3; k++)
			cout <<"ERROR 666: Domain error, reselect N" << endl;
	}
	double pmlWidth, sigM, sigMax, kappaMax, kappaM, aMax, aM;
	dz = ( dom.zMax - dom.zMin )/ ( EMpoints - 1);

	Config settings("settings.txt");
	settings.get("Maxwell/evolve", evolve);
	settings.get("Maxwell/pmlPoints", PML);
	settings.get("Maxwell/sigMax", sigMax);
	settings.get("Maxwell/sigM", sigM);
	settings.get("Maxwell/kappaMax", kappaMax);
	settings.get("Maxwell/kappaM", kappaM);
	settings.get("Maxwell/aMax", aMax);
	settings.get("Maxwell/aM", aM);

	double freq, ratio;
	settings.get("Maxwell/pulse", pulseSetting);
	settings.get("Maxwell/frequency", freq);
	settings.get("Maxwell/amplitude", pAmp);
	settings.get("Maxwell/width", width);
	settings.get("Maxwell/delayRatio", ratio);
	settings.get("Maxwell/addlpoints", addlpoints);
	addlpoints = addlpoints*2;
	omega = 2*PI*freq;
	dt = dom.r_dt;
	maxT = ratio*width;

	// Real constants for Maxwell's equations
	a1 = 1/MUNOT;
	a2 = 1/EPSNOT;
	a3 = (HBAR*E) / (ME*EPSNOT*2)*dom.fieldCoef;
	a4 = -(HBAR*pow(E,2.0))/(ME*EPSNOT);
	a5 = -E / EPSNOT;

	// Point arrangement:
	//  |---PML---|  buffer  |-----EMpts=----|  buffer  |---PML---|    
	//   x   o   x   o   x   o   x   o   x   o   x   o   x   o   x
	//   E   H   E   H   E   H   E   H   E   H   E   H   E   H   E
	
	sourcePos = PML+1;

	NEx = ( EMpoints - 1 ) + addlpoints + 2*PML + 2;
	NHy = EMpoints + addlpoints + 2*PML;
	
	interpN = 2*PML + ( EMpoints - 1 );
	psiStart = addlpoints/2 + 1;

	fieldDom = new double [interpN];
	fieldDom[0] =  dz/2.0 - dz* (double) PML;

	for (int k = 1; k < interpN; k++)
		fieldDom[k] = fieldDom[k-1] + dz;

	Ex   = new double [NEx];
	A    = new double [NEx];
	Aold = new double [NEx];
	Hy   = new double [NHy];
	QEx  = new double [NHy];
	QHy  = new double [NEx];
	
	
	bx = new double [NEx];
	gx = new double [NEx];
	kx = new double [NEx];

	by = new double [NHy];
	gy = new double [NHy];
	ky = new double [NHy];
	
	for (int k = 0; k < NEx; k++)
	{
		Ex[k] = 0.0;
		A[k] = 0.0;
		Aold[k] = 0.0;
		QHy[k] = 0.0;

		bx[k] = 0.0;
		gx[k] = 0.0;
		kx[k] = 0.0;
	}

	for (int k = 0; k < NHy; k++)
	{
		Hy[k] = 0.0;
		QEx[k] = 0.0;

		by[k] = 0.0;
		gy[k] = 0.0;
		ky[k] = 0.0;
	}	

	double xVec, sigVec, sigVecStag, kappaVec, aVec, aVecStag;		// PML calculation constants
	double kappaVecStag, tau, tauStag;							    // PML calculation constants

	// PML computations
	int m = PML - 1;
	pmlWidth = PML;

	double* tempbx = new double [PML];
	double* tempby = new double [PML];
	double* tempgx = new double [PML];
	double* tempgy = new double [PML];
	double* tempkx = new double [PML];
	double* tempky = new double [PML];
	for (int k = 0; k < PML; k++)
	{
		xVec = m + .5;
		m = m-1;

		sigVec		= sigMax * pow(xVec + .5, sigM)/(pow(pmlWidth, sigM));
		sigVecStag  = sigMax * pow(xVec, sigM)/(pow(pmlWidth, sigM));

		kappaVec	= 1 + (kappaMax-1)* pow(xVec + .5, kappaM)/(pow(pmlWidth, kappaM));
		kappaVecStag= 1 + (kappaMax-1)* pow(xVec , kappaM)/(pow(pmlWidth, kappaM));

		aVec		= aMax * pow(xVec + .5, aM)/(pow(pmlWidth, aM));
		aVecStag	= aMax * pow(xVec, aM)/(pow(pmlWidth, aM));

		tau			= (kappaVec*EPSNOT)/(kappaVec*aVec + sigVec);
		tauStag		= (kappaVecStag*EPSNOT)/(kappaVecStag*aVecStag + sigVecStag);
		 
		tempbx[k] = exp(-dom.r_dt/tau);
		tempby[k] = exp(-dom.r_dt/tauStag);

		tempgx[k] = (sigVec/(kappaVec*(kappaVec*aVec+sigVec)))*(1-tempbx[k]);
		tempgy[k] = (sigVecStag/(kappaVecStag*(kappaVecStag*aVecStag+sigVecStag)))*(1-tempby[k]);

		tempkx[k] = kappaVec;
		tempky[k] = kappaVecStag;
	}

	for (int k = 0; k < PML; k++)
	{
		bx[k] = tempbx[k];
		by[k] = tempby[k];
		gx[k] = tempgx[k];
		gy[k] = tempgy[k];
		kx[k] = tempkx[k];
		ky[k] = tempky[k];
	}

	int s = PML-1;
	for (int k = NHy-PML; k < NHy; k++)
	{
		by[k] = tempby[s];
		gy[k] = tempgy[s];
		ky[k] = tempky[s];
		s--;
	}

	s = PML-1;
	for (int k = NEx-PML; k < NEx; k++)
	{
		bx[k] = tempbx[s];
		gx[k] = tempgx[s];
		kx[k] = tempkx[s];
		s--;
	}

	cout << endl;
}

Energy::Energy()
{
	Ue		= 0.0;
	Uion	= 0.0;
	Ubond	= 0.0;
	Te		= 0.0;
	Tion	= 0.0;
}

Manager::Manager(){
	
	Config settings("settings.txt");
	settings.get("Output/frames", frames);
	settings.get("Output/outputParameters", outputParameters);
	settings.get("Output/video", videoValue);

	psiBin = "psiVid.bin";
	rhoBin = "rhoVid.bin";
	ionBin = "ionVid.bin";
	fdtdBin = "fdtd.bin";

}