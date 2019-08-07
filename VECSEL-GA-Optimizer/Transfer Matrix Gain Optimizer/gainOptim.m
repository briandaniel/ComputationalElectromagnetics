function [gainInv] = gainOptim(lambda, n, L, sigma)

%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

% solution with gain
[R, r, T, t] = transferMatrixOne(lambda, n, L, sigma);

%% Calculations of group delay
% lambdaStep = lambdas(end)-lambdas(end-1);
% phi = phase(r);
% % 
% % ppPHI = pchip(lambdas(1:2:end), phi(1:2:end));
% % 
% % PHI = ppval(lambdas, ppPHI);
% % phi = PHI
% % 
% phiDot = gradient(phi,lambdaStep);
% for k = 3:lambdaSteps-20;
%     
%    meanPhase = abs(mean(phiDot(1:k+10)));
%       
%    if abs(phiDot(k)) > 5*meanPhase
%        phiDot(k) = (phiDot(k+2) + phiDot(k-2))/2;
%    end
%           
%     
% end
% lambdaStep = lambdas(end)-lambdas(end-1);
% ppPHIdot = pchip(lambdas(1:2:end), phiDot(1:2:end));
% PHIdot = ppval(lambdas, ppPHIdot);
% phiDot = PHIdot
% tau = -lambdas.^2./(2.*pi.*cNot).*phiDot;
% tauPico = tau.*10e12;

%% Calculations of group delay for solution with no gain
% lambdaStep2 = lambdas(end)-lambdas(end-1);
% phi2 = phase(r2);
% % 
% % ppPHI = pchip(lambdas(1:2:end), phi(1:2:end));
% % 
% % PHI = ppval(lambdas, ppPHI);
% % phi = PHI
% % 
% phiDot2 = gradient(phi2,lambdaStep2);
% for k = 3:lambdaSteps-20;
%     
%    meanPhase = abs(mean(phiDot2(1:k+10)));
%       
%    if abs(phiDot2(k)) > 5*meanPhase
%        phiDot2(k) = (phiDot2(k+2) + phiDot2(k-2))/2;
%    end
%           
%     
% end
% lambdaStep2 = lambdas(end)-lambdas(end-1);
% ppPHIdot2 = pchip(lambdas(1:2:end), phiDot2(1:2:end));
% PHIdot2 = ppval(lambdas, ppPHIdot);
% phiDot2 = PHIdot2
% tau2 = -lambdas.^2./(2.*pi.*cNot).*phiDot2;
% tauPico2 = tau2.*10e12;

%%
gain = (T+R-1).*100;
% normalGain = norm(gain,inf);
gainInv = 1/gain;
% 
% L;
% lambdas(index);
% 1/gainInv;

end