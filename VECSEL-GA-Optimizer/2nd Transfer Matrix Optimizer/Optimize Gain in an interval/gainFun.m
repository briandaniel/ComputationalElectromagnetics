function [lambdas, R, T, tauPicoR, tauPicoT, phaseR, phaseT] = gainFun(lambdaMin, lambdaMax, lambdaSteps, n, L, sigma)

%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

% solution with gain
[lambdas, R, r, T, t] = transferMatrix(lambdaMin, lambdaMax, lambdaSteps, n, L, sigma);

%% Calculations of group delay
lambdaStep = lambdas(end)-lambdas(end-1);
phi = phase(r);
phaseR = phi;
phiDot = gradient(phi,lambdaStep);

for k = 3:lambdaSteps-20;
    
   meanPhase = abs(mean(phiDot(1:k+10)));
      
   if abs(phiDot(k)) > 5*meanPhase
       phiDot(k) = (phiDot(k+2) + phiDot(k-2))/2;
   end
          
    
end
lambdaStep = lambdas(end)-lambdas(end-1);
ppPHIdot = pchip(lambdas(1:2:end), phiDot(1:2:end));
PHIdot = ppval(lambdas, ppPHIdot);
phiDot = PHIdot;
tau = -lambdas.^2./(2.*pi.*cNot).*phiDot;
tauPicoR = tau.*10e12;

%% Calculations of group delay for solution with no gain
lambdaStep2 = lambdas(end)-lambdas(end-1);
phi2 = phase(t);
phaseT = phi2;

phiDot2 = gradient(phi2,lambdaStep2);
for k = 3:lambdaSteps-20;
    
   meanPhase = abs(mean(phiDot2(1:k+10)));
      
   if abs(phiDot2(k)) > 5*meanPhase
       phiDot2(k) = (phiDot2(k+2) + phiDot2(k-2))/2;
   end
          
    
end
lambdaStep2 = lambdas(end)-lambdas(end-1);
ppPHIdot2 = pchip(lambdas(1:2:end), phiDot2(1:2:end));
PHIdot2 = ppval(lambdas, ppPHIdot2);
phiDot2 = PHIdot2;
tau2 = -lambdas.^2./(2.*pi.*cNot).*phiDot2;
tauPicoT = tau2.*10e12;


end