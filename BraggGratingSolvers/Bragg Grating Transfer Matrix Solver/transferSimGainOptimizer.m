close all
%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

% Grating parameters
lambdaMax = 500.1e-9;
nH = 2.10;
nL = 1.59;
lH = lambdaMax/(4*nH);
lL = lambdaMax/(4*nL);
grateNum = 1;                   % Number of grating sets
grateWidth = lH*grateNum + lL*grateNum;

% Lossy
sigH = -5000;
sigL = -5000;
sigH2 = 0;
sigL2 = 0;

lambdaWindow = [ 300, 800]*1e-9;
lambdaSteps = 1000;

N = grateNum;
n = [1, nH, repmat([nL,nH], 1, N), 1];        % indices for the layers A|H(LH)N |G 
L = [lH, repmat([lL,lH], 1, N)];                % Lengths of layers
sigma = [0, sigH, repmat([sigL,sigH], 1, N), 0];    
sigma2 = [0, sigH2, repmat([sigL2,sigH2], 1, N), 0];    

% solution with gain
[lambdas, R, r, T, t] = transferMatrix(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sigma);
% solution with no gain
[lambdas2, R2, r2, T2, t2] = transferMatrix(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sigma2);


%% Calculations of group delay
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
lambdaStep = lambdas(end)-lambdas(end-1);
ppPHIdot = pchip(lambdas(1:2:end), phiDot(1:2:end));
PHIdot = ppval(lambdas, ppPHIdot);
phiDot = PHIdot
tau = -lambdas.^2./(2.*pi.*cNot).*phiDot;
tauPico = tau.*10e12;

%% Calculations of group delay for solution with no gain
lambdaStep2 = lambdas(end)-lambdas(end-1);
phi2 = phase(r2);
% 
% ppPHI = pchip(lambdas(1:2:end), phi(1:2:end));
% 
% PHI = ppval(lambdas, ppPHI);
% phi = PHI
% 
phiDot2 = gradient(phi2,lambdaStep2);
for k = 3:lambdaSteps-20;
    
   meanPhase = abs(mean(phiDot2(1:k+10)));
      
   if abs(phiDot2(k)) > 5*meanPhase
       phiDot2(k) = (phiDot2(k+2) + phiDot2(k-2))/2;
   end
          
    
end
lambdaStep2 = lambdas(end)-lambdas(end-1);
ppPHIdot2 = pchip(lambdas(1:2:end), phiDot2(1:2:end));
PHIdot2 = ppval(lambdas, ppPHIdot);
phiDot2 = PHIdot2
tau2 = -lambdas.^2./(2.*pi.*cNot).*phiDot2;
tauPico2 = tau2.*10e12;


%%
figure('outerposition', [200. 200, 1000, 1000])
subplot(2,1,1)
hold on;
plot(lambdas*1e9, R, lambdas*1e9, T, lambdas*1e9, R2, lambdas*1e9, T2)
axis([lambdaWindow(1)*1e9 lambdaWindow(2)*1e9 0 1.1])
xlabel('{\lambda}nm')
ylabel('Spectra')
axis([lambdas(1)*1e9, lambdas(end).*1e9, 0, max(max(T),max(R))]);
subplot(2,1,2)
plot(lambdas*1e9,tauPico,'g', lambdas*1e9,tauPico2, 'r')
title('group delay')
ylabel('ps')
xlabel('{\lambda}nm')
legend('Reflected Gain Wave', 'Transmitted Gain Wave', 'Reflected Wave', 'Transmitted Wave')

figure('outerposition', [1000. 500, 1000, 500])
% plot(lambdas*1e9, R+T, lambdas*1e9, R2+T2)
gain = (R+T - 1)*100;
plot(lambdas*1e9, gain)
ylabel('% gain')
xlabel('{\lambda}nm')
title('Total gain (Transmitted+reflected)')


