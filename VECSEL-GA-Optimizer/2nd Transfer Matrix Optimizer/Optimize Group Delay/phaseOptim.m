function [gainSumInv] = phaseOptim(lambdaMin, lambdaMax, lambdaSteps, n, L, sigma)

%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

% solution with gain
[lambdas, R, r, T, t] = transferMatrix(lambdaMin, lambdaMax, lambdaSteps, n, L, sigma);

%% Calculations of group delay
phi = phase(r);
phaseR = phi;
lambdaStep = lambdas(end)-lambdas(end-1);
phiDot = gradient(phi,lambdaStep);
ppPHIdot = pchip(lambdas(1:2:end), phiDot(1:2:end));
PHIdot = ppval(lambdas, ppPHIdot);
phiDot = PHIdot;
tau = -lambdas.^2./(2.*pi.*cNot).*phiDot;
tauPicoR = tau.*10e12;
lambdaStep = (lambdas(end)-lambdas(end-1))*1e9;
gainSumInv = length(lambdas)/(abs(sum(tauPicoR)));


end