clear all
close all

% Grating parameters
lambdaMax = 500.1e-9;
nH = 2.10;
nL = 1.59;
lH = lambdaMax/(4*nH);
lL = lambdaMax/(4*nL);
grateNum = 5;                   % Number of grating sets
grateWidth = lH*grateNum + lL*grateNum;
lambdaWindow = [ 450, 600]*1e-9;
lambdaSteps = 500;

% Lossy
sigH = -100;
sigL = -100;
sigH2 = 0;
sigL2 = 0;

% Vectors
N = grateNum;
n = [1, nH, repmat([nL,nH], 1, N), 1];        % indices for the layers A|H(LH)N |G 
L = [lH, repmat([lL,lH], 1, N)];                % Lengths of layers
sigma = [0, sigH, repmat([sigL,sigH], 1, N), 0];    
sigma2 = [0, sigH2, repmat([sigL2,sigH2], 1, N), 0];    



% Optimization calling function 'gainOptim.m'
gainFunNew = @(x) gainOptim(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, x, sigma, sigma2);
% options = optimset('tolx', 1e-6);

variables = length(L);
LB = ones(size(L))*1e-8;
UB = ones(size(L))*50e-8;
optimalLength = ga(gainFunNew, variables, [], [], [], [], LB, UB);


[gain, lambdas, R, T, R2, T2, tauPico, tauPico2] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, optimalLength, sigma, sigma2);
[gainB, lambdasB, RB, TB, R2B, T2B, tauPicoB, tauPico2B] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sigma, sigma2);


percentGainDifference = 1/gainFunNew(optimalLength) - 1/gainFunNew(L);
display(strcat('Average percent gain difference = ',  num2str(round(percentGainDifference*1e2)*1e-2), '%'))

%%
figure('outerposition', [200. 200, 1000, 1000])
subplot(2,1,1)
hold on;
plot(lambdas*1e9, R2, lambdas*1e9, T2, lambdas*1e9, R2B, lambdas*1e9, T2B)
axis([lambdaWindow(1)*1e9 lambdaWindow(2)*1e9 0 1.1])
xlabel('{\lambda}nm')
ylabel('Spectra')
axis([lambdas(1)*1e9, lambdas(end).*1e9, 0, max(max(T),max(R))]);
legend('Reflected Gain Wave', 'Transmitted Gain Wave', 'Original Reflected Gain', 'Original Transmitted Gain')

subplot(2,1,2)
plot(lambdas*1e9,tauPico2,'g', lambdas*1e9,tauPico2B, 'r')
axis([lambdas(1)*1e9, lambdas(end).*1e9, min(min(tauPico2),min(tauPico2B)), max(max(tauPico2),max(tauPico2B))]);
title('group delay')
ylabel('ps')
xlabel('{\lambda}nm')
legend('Optimized time delay', 'Original tiem delay')

figure('outerposition', [1000. 500, 1000, 500])
% plot(lambdas*1e9, R+T, lambdas*1e9, R2+T2)
plot(lambdas*1e9, gain, lambdas*1e9, gainB)
ylabel('% gain')
xlabel('{\lambda}nm')
title('Total gain (Transmitted+reflected)')
legend('Optimized Gain', 'Original Gain')


