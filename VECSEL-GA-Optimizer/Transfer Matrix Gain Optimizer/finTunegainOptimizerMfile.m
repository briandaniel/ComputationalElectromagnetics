close all

k = 0;
while k < inf;
    close all
    
    figure
    plot(1,1)
   k = k+1

%% SETUP
tic
% Grating parameters
lambdaMax = 721.42e-9;
nH = 4.00;
nL = 3.10;
lH = lambdaMax/(4*nH);
lL = lambdaMax/(4*nL);
grateNum = 10;                   % Number of grating sets
grateWidth = lH*grateNum + lL*grateNum;
lambdaWindow = [ 645, 655]*1e-9;
lambdaSteps = 3;
lambdas = linspace(lambdaWindow(1), lambdaWindow(2), lambdaSteps);
index = round(lambdaSteps/2);

% Lossy
sigH = 0;
sigL = -600;
sigH2 = 0;
sigL2 = 0;

% Vectors
N = grateNum;
n = [1, nH, repmat([nL,nH], 1, N), 1];        % indices for the layers A|H(LH)N |G 
L = [lH, repmat([lL,lH], 1, N)];                % Lengths of layers
sigma = [0, sigH, repmat([sigL,sigH], 1, N), 0];    
sigma2 = [0, sigH2, repmat([sigL2,sigH2], 1, N), 0];    

grateSize = sum(L);
cellWidth = grateSize/length(L); 

lambdaOptim = 650e-9;

optimalLength2 = 1.0e-06 * [0.0414
    0.0563
    0.0427
    0.0415
    0.0608
    0.0414
    0.0397
    0.0563
    0.0414
    0.0941
    0.0397
    0.1283
    0.0608
    0.0608
    0.0414
    0.0686
    0.0414
    0.0433
    0.0563
    0.0414
    0.0414];

%% OPTIMIZATION
% Optimization calling function 'gainOptim.m'
gainFunNew = @(x) gainOptim(lambdaOptim, n, x, sigma);

variables = length(L);
LB = ones(size(L))*1e-8;
UB = ones(size(L))*cellWidth*2.5;

options = gaoptimset('tolfun', 1e-2, 'populationsize', 5, 'timelimit', inf, ...
                     'display', 'off', 'stalltimelimit', inf, 'useparallel', 'always', 'crossoverfraction', .9, 'initialPopulation', optimalLength2');
A = ones(size(L));
optimalLength = ga (gainFunNew, variables, A, grateSize, [], [], LB, UB, [], options );

% optimalLength = fminsearch(gainFunNew, optimalLength2);
newGain = 1/gainFunNew(optimalLength)
newOptimum = optimalLength;

if currentGain < newGain
    
    currentGain = newGain
    currentOptimum = newOptimum;
    storage = [storage ; currentOptimum]
end
newGain
currentGain
toc
% end

end
%%
close all
optimalLength = currentOptimum;
% Perturb results
% L2 = optimalLength;
% L2 = L2+mean(L2).*(rand(1,length(L2))-.5).*.02;
% optimalLength = L2;

lambdaWindow = [ 645, 655]*1e-9;
lambdaSteps = 10000;

[gain, lambdas, R, T, R2, T2, tauPico, tauPico2] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, optimalLength, sigma, sigma2);
[gainB, lambdasB, RB, TB, R2B, T2B, tauPicoB, tauPico2B] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sigma, sigma2);

newGain = 1/gainFunNew(optimalLength);

% RUN WITH ONLY SINGLE LAYER
% Single layer is the same size as total optimized structure
newGrateSize = sum(optimalLength);

nC = [1, nH, 1];
LC = [newGrateSize];
sigmaC = [0, sigH, 0];
sigma2C = [0, sigH, 0];
[gainC, lambdasC, RC, TC, R2C, T2C, tauPicoC, tauPico2C] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, nC, LC, sigmaC, sigma2C);
display(' ')
display(' ')
oldGain = max(TB)*100;


display(strcat('Optimized wavelength = ', num2str(round(lambdaOptim*1e11)*1e-2), 'nm'))
display(' ')
display(strcat('Optimized Gain= ',  num2str(round(newGain*1e2)*1e-2), '%,             |  ', ...
               'Lambda/4 Gain= ',  num2str(round(oldGain*1e2)*1e-2), '%'))
           
percentGainDifference = newGain - oldGain;
relativeGainIncrease = (percentGainDifference)/(oldGain-100)*100;

display(strcat('Percent gain difference = ',  num2str(round(percentGainDifference*1e2)*1e-2), '%,     |  ', ...
               'Relative gain increase = ',  num2str(round(relativeGainIncrease*1e1)*1e-1), '%'))
% 
display(' ')
% PLOTS

figure('outerposition', [200. 200, 1000, 1000])
subplot(2,1,1)
hold on;
plot(lambdas*1e9, R, 'g-')
% plot(lambdas*1e9, T, 'b-')
plot(lambdas*1e9, RB, 'g-.')
% plot(lambdas*1e9, TB, 'b-.')
axis([lambdaWindow(1)*1e9 lambdaWindow(2)*1e9 0 2.0])
xlabel('{\lambda}nm')
ylabel('Spectra')
% axis([lambdas(1)*1e9, lambdas(end).*1e9, 0, max([max(T),max(R),max(R2),max(T2)])]);
legend('Optimized Reflected Wave', 'Optimized Transmitted Wave', '{\lambda}/4 Reflected Wave', '{\lambda}/4 Transmitted Wave')

subplot(2,1,2)
plot(lambdas*1e9,tauPico2,'g', lambdas*1e9,tauPico2B, 'r', lambdas(index)*1e9, tauPico2(index), 'ro')
axis([lambdas(1)*1e9, lambdas(end).*1e9, min(min(tauPico2),min(tauPico2B)), max(max(tauPico2),max(tauPico2B))]);
title('group delay')
ylabel('ps')
xlabel('{\lambda}nm')
legend('Optimized time delay', '{\lambda}/4 band edge time delay')

figure('outerposition', [1000. 500, 1000, 500])
% plot(lambdas*1e9, R+T, lambdas*1e9, R2+T2)
plot(lambdas*1e9, gain, lambdas*1e9, gainB)
ylabel('% gain')
xlabel('{\lambda}nm')
title('Total gain (Transmitted+reflected)')
legend('Optimized Gain', '{\lambda}/4 band edge Gain')
axis([lambdas(1)*1e9 lambdas(end)*1e9, 0, max(gain)*1.2])


figure('outerposition', [1500. 10, 1000, 500])
leg = 0;
stairs = zeros(1,2*length(optimalLength));
base = stairs;
base(1) = 0;
for k = 1:length(optimalLength);
    leg = optimalLength(k)+leg;
    base(2*k:2*k+1) = k;
    
    if k < length(optimalLength);
        stairs(2*k+1:2*k+2) = leg;    
    end
end

base(end) = length(optimalLength);

stairs(length(optimalLength)*2+1) = stairs(end) + optimalLength(end);

plot(stairs*1e6, base, 'k-')
ylabel('Layer Number')
xlabel('Widths ({\mu}m)')
axis([min(stairs*1e6), max(stairs*1e6),0, length(optimalLength)*1.2])
hold on


rectangle('position', [0, 0, optimalLength(1)*1e6, 1], 'facecolor', [0 .5 .7])

for k = 1:length(optimalLength)-1
    if mod(k,2) == 0;
    rectangle('position', [stairs(2*k+1)*1e6, 0, optimalLength(k+1)*1e6, k+1], 'facecolor', [0 .5 .7])
    else
    rectangle('position', [stairs(2*k+1)*1e6, 0, optimalLength(k+1)*1e6, k+1], 'facecolor', [0 0 .9])
    end
end

title('Profile of layering')

    
    
hold off
    
toc
