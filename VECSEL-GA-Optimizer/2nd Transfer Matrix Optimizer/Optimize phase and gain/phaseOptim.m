function [gainSumInv, phaseCurveNorm] = phaseOptim(lambdaMin, lambdaMax, lambdaSteps, n, L, sigma)

%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

% solution with gain
[lambdas, R, r, T, t] = transferMatrix(lambdaMin, lambdaMax, lambdaSteps, n, L, sigma);

%% Calculations of group delay
lambdaStep = (lambdas(end)-lambdas(end-1))*1e9;
phi = phase(r);
phaseR = phi;
phaseCurve = (phi(3:end) - 2.*phi(2:end-1) + phi(1:end-2))./(lambdaStep).^2;
phaseCurveNorm = norm(phaseCurve,2);

gainSumInv = length(lambdas)/sum(R*100);



end