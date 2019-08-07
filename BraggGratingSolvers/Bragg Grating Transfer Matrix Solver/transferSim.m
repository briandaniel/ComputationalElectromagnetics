close all

% Grating parameters
lambdaMax = 500.1e-9;
nH = 2.10;
nL = 1.59;
lH = lambdaMax/(4*nH);
lL = lambdaMax/(4*nL);
grateNum = 1;                   % Number of grating sets
grateWidth = lH*grateNum + lL*grateNum;

% Lossy
sigH = 5000;
sigL = 0;

lambdaWindow = [ 300, 800]*1e-9;
lambdaSteps = 1000;


% Grating in active region
numWell = 14;                   % Number of quantum wells
nQ = 3.7;                       % Refractive index of quantum wells
nCH = 3.40;                     % Refractive index of high spacers
nCL = 3.10;                     % Refractive index of low spacers
lWell = 8e-9;                   % Widths of quantum wells in meters (must be fixed)
lCladH = 69.8e-9;               % Widths of surrounding cladding (High refractive index)
lCladL = 76.6e-9;               % Widths of surrounding cladding (Low refractive index)
sigQ = 0 ;                      % Sigmas approximating gain of quantum well
sigmaQ = repmat([sigQ, 0, 0 ], 1, numWell);   
LActive = repmat([lWell, lCladH, lCladL ], 1, numWell);
nActive = repmat([nQ, nCH, nCL],1,numWell);

% Grating in DBR
lbdMaxRef = 950e-9;             % Wavelength of maximum reflection
grateNum = 22;                  % Number of bilayers in DBR
nH = 3.40;                      % High refractive index
nL = 3.10;                      % Low refractive index
lH = lbdMaxRef/(4*nH);          % Quarter wavelength widths of high refractive index
lL = lbdMaxRef/(4*nL);          % Quarter wavelength widths of low refractive index
LDBR = [lH, repmat([lL,lH], 1, grateNum)]; 
nDBR = [nH, repmat([nL,nH], 1, grateNum)];

% Total index and grating vectors
nInit = 1;                      % Refractive index before structure
nEnd = .5;                      % Refractive index after structure
L = [LActive, LDBR];
n = [nInit, nActive, nDBR, nEnd];
sig = [0 sigmaQ zeros(1,length(LDBR)), 0];

[lambdas, R, r, T, t] = transferMatrix(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sig);


lambdaStep = lambdas(end)-lambdas(end-1);

phi = phase(r);
% 
% ppPHI = pchip(lambdas(1:2:end), phi(1:2:end));
% 
% PHI = ppval(lambdas, ppPHI);
% phi = PHI
% 
phiDot = gradient(phi,lambdaStep);


for k = 3:lambdaSteps-20;
    
   meanPhase = abs(mean(phiDot(1:k+10)));
   
   
   if abs(phiDot(k)) > 5*meanPhase
       phiDot(k) = (phiDot(k+2) + phiDot(k-2))/2;
   end
        
    
    
end

% lambdaStep = lambdas(end)-lambdas(end-1);
% 
% ppPHIdot = pchip(lambdas(1:2:end), phiDot(1:2:end));
% 
% PHIdot = ppval(lambdas, ppPHIdot);
% phiDot = PHIdot
% 
% 
% 
% tau = -lambdas.^2./(2.*pi.*cNot).*phiDot;
% tauPico = tau.*10e12;

figure('outerposition', [200. 200, 1000, 1000])
subplot(2,1,1)
hold on;
plot(lambdas*1e9, R, lambdas*1e9, T)
axis([lambdaWindow(1)*1e9 lambdaWindow(2)*1e9 0 1.1])
xlabel('{\lambda}nm')
ylabel('Spectra')
subplot(2,1,2)
title('group delay')
ylabel('ps')
xlabel('{\lambda}nm')
legend('Reflected Wave', 'Transmitted Wave')