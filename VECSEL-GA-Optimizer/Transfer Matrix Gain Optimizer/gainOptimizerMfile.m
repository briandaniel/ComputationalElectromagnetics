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
lambdaMax = 721.22e-9;
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
sigL = -500;
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

%% OPTIMIZATION
% Optimization calling function 'gainOptim.m'
gainFunNew = @(x) gainOptim(lambdaOptim, n, x, sigma);

variables = length(L);
LB = ones(size(L))*1e-8;
UB = ones(size(L))*cellWidth*2.5;

options = gaoptimset('tolfun', 1e-2, 'populationsize', 10+size(storageBin), 'timelimit', inf, ...
                     'display', 'off', 'stalltimelimit', inf, 'useparallel', 'always', 'crossoverfraction', .9, ...
                     'initialPopulation', storageBin);
A = ones(size(L));
% optimalLength = ga (gainFunNew, variables, A, grateSize, [], [], LB, UB, [], options );

optimalLength = ga (gainFunNew, variables, [], [], [], [], LB, UB, @(x)consFun(x,grateSize), options );

%%
% optimalLength = 1.0e-06.*[0.0446
%     0.0435
%     0.0435
%     0.0444
%     0.0510
%     0.0444
%     0.0435
%     0.1066
%     0.0510
%     0.1281
%     0.0444
%     0.0435
%     0.0629
%     0.0444
%     0.0435
%     0.0435
%     0.0435
%     0.0435
%     0.0444
%     0.0629
%     0.0404]'

% optimalLength = 1.0e-07.*[    0.4679
%     0.4620
%     0.4620
%     0.4620
%     0.4620
%     0.4620
%     0.4711
%     0.462
%     0.462
%     0.4620
%     0.4620
%     0.4711
%     0.4620
%     0.4620
%     0.4620
%     0.4620
%     0.4620
%     0.4620
%     0.4620
%     0.4711
%     0.4620
%     0.4620
%     0.8239
%     0.4580
%     0.4620
%     0.4620
%     0.4620
%     0.4711
%     0.4620
%     0.4711
%     0.4620
%     0.4620
%     0.4711
%     0.4620
%     0.4620
%     0.4620
%     0.4620
%     0.4620
%     0.462
%     0.462
%     0.4620
%     0.462
%     0.462
%     0.8239
%     0.462
%     0.4620
%     0.8239
%     0.462
%     0.462
%     0.462
%     0.462
%     0.4620
%     0.4620
%     0.4620
%     0.4711
%     0.462
%     0.4711
%     0.4620
%     0.4711
%     0.8239
%     0.5035];

% 
% 
% optimalLength = 1.0e-07.*[0.0439
%     0.0439
%     0.0439
%     0.0439
%     0.0439
%     0.0439
%     0.0439
%     0.0439
%     0.0439
%     0.1282
%     0.0439
%     0.1282
%     0.0439
%     0.0439
%     0.0439
%     0.0456
%     0.0439
%     0.0439
%     0.0439
%     0.0439
%     0.0439]
%     
%%
newGain = 1/gainFunNew(optimalLength);
newOptimum = optimalLength;

if currentGain < newGain
    
    currentGain = newGain;
    currentOptimum = newOptimum;
    storage = [storage ; currentOptimum]
end

if length(storage(:,1)) > 2;
    storage = currentOptimum;
end


if k > 2000;
    storageBin = [storageBin; storage];
    gainStorage = [gainStorage; currentGain];
    
    storage = zeros(1,length(currentOptimum));
    currentOptimum = storage;
    currentGain = 0;
    k = 0;
    
end


newGain
currentGain
toc
end

%%
close all

% Perturb results
optimalLength = currentOptimum;
% L2 = optimalLength;
% L2 = L2+mean(L2).*(rand(1,length(L2))-.5).*.02;
% optimalLength = L2;

lambdaWindow = [ 500, 800]*1e-9;
lambdaSteps = 2000;

[gain, lambdas, R, T, R2, T2, tauPicoR, tauPicoT] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, optimalLength, sigma, sigma2);
[gainB, lambdasB, RB, TB, R2B, T2B, tauPicoRB, tauPicoTB] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sigma, sigma2);

newGain = 1/gainFunNew(optimalLength);

% RUN WITH ONLY SINGLE LAYER
% Single layer is the same size as total optimized structure
newGrateSize = sum(optimalLength);

nC = [1, nH, 1];
LC = [newGrateSize];
sigmaC = [0, sigH, 0];
sigma2C = [0, sigH, 0];
[gainC, lambdasC, RC, TC, R2C, T2C, tauPicoRC, tauPicoTC] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, nC, LC, sigmaC, sigma2C);
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
plot(lambdas*1e9, T, 'b-')
plot(lambdas*1e9, RB, 'g-.')
plot(lambdas*1e9, TB, 'b-.')
% axis([lambdaWindow(1)*1e9 lambdaWindow(2)*1e9 0 2.0])
xlabel('{\lambda}nm')
ylabel('Spectra')
% axis([lambdas(1)*1e9, lambdas(end).*1e9, 0, max([max(T),max(R),max(R2),max(T2)])]);
legend('Optimized Reflected Wave', 'Optimized Transmitted Wave', '{\lambda}/4 Reflected Wave', '{\lambda}/4 Transmitted Wave')

subplot(2,1,2)
hold on
plot(lambdas*1e9,tauPicoR,'g', lambdas*1e9,tauPicoRB, 'g--')
plot(lambdas*1e9,tauPicoT,'r', lambdas*1e9,tauPicoTB, 'r--')
hold off
% axis([lambdas(1)*1e9, lambdas(end).*1e9, min(min(tauPico2),min(tauPico2B)), max(max(tauPico2),max(tauPico2B))]);
title('group delay')
ylabel('ps')
xlabel('{\lambda}nm')
legend('Optimized time delay (reflection)', '{\lambda}/4 band edge time delay (reflection)', 'Optimized time delay (transmission)', '{\lambda}/4 band edge time delay (transmission)')

figure('outerposition', [1000. 500, 1000, 500])
% plot(lambdas*1e9, R+T, lambdas*1e9, R2+T2)
plot(lambdas*1e9, gain, lambdas*1e9, gainB)
ylabel('% gain')
xlabel('{\lambda}nm')
title('Total gain (Transmitted+reflected)')
legend('Optimized Gain', '{\lambda}/4 band edge Gain')
axis([lambdas(1)*1e9 lambdas(end)*1e9, 0, max(gain)*1.2])








% figure('outerposition', [1500, 0, 1000, 700])
% subplot(2,1,2)
% leg = 0;
% stairs = zeros(1,2*length(optimalLength));
% base = stairs;
% base(1) = 0;
% for k = 1:length(optimalLength);
%     leg = optimalLength(k)+leg;
%     base(2*k:2*k+1) = k;
%     
%     if k < length(optimalLength);
%         stairs(2*k+1:2*k+2) = leg;    
%     end
% end
% 
% base(end) = length(optimalLength);
% stairs(length(optimalLength)*2+1) = stairs(end) + optimalLength(end);
% plot(stairs*1e6, base, 'k-')
% ylabel('Layer Number')
% xlabel('Widths ({\mu}m)')
% 
% a = min(stairs*1e6);
% b = max(stairs*1e6);
% c = length(optimalLength)*1.1;
% 
% axis([a, b, 0, c])
% hold on
% 
% 
% 
% rectangle('position', [0, 0, optimalLength(1)*1e6, 1], 'facecolor', [0 .5 .7])
% for k = 1:length(optimalLength)-1
%     if mod(k,2) == 0;
%     rectangle('position', [stairs(2*k+1)*1e6, 0, optimalLength(k+1)*1e6, k+1], 'facecolor', [0 .5 .7])
%     else
%     rectangle('position', [stairs(2*k+1)*1e6, 0, optimalLength(k+1)*1e6, k+1], 'facecolor', [0 0 .9])
%     end
% end
% title('Profile of layering for {\lambda}/L')   
% hold off
% 
% 
% 
% subplot(2,1,1)
% leg = 0;
% stairs = zeros(1,2*length(optimalLength));
% base = stairs;
% base(1) = 0;
% for k = 1:length(optimalLength);
%     leg = L(k)+leg;
%     base(2*k:2*k+1) = k;
%     
%     if k < length(optimalLength);
%         stairs(2*k+1:2*k+2) = leg;    
%     end
% end
% 
% base(end) = length(optimalLength);
% stairs(length(optimalLength)*2+1) = stairs(end) + L(end);
% plot(stairs*1e6, base, 'k-')
% ylabel('Layer Number')
% xlabel('Widths ({\mu}m)')
% axis([a, b, 0, c])
% hold on
% 
% rectangle('position', [0, 0, L(1)*1e6, 1], 'facecolor', [0 .5 .7])
% for k = 1:length(optimalLength)-1
%     if mod(k,2) == 0;
%     rectangle('position', [stairs(2*k+1)*1e6, 0, L(k+1)*1e6, k+1], 'facecolor', [0 .5 .7])
%     else
%     rectangle('position', [stairs(2*k+1)*1e6, 0, L(k+1)*1e6, k+1], 'facecolor', [0 0 .9])
%     end
% end
% title('Profile of layering for optimized structure')   
% hold off
% 
% 
% 
% toc


% Alternative widths plot
figure('outerposition', [1500, 0, 1000, 500])
leg = 0;
stairs = zeros(1,2*length(optimalLength));
base = stairs;
base(1) = 0;
for k = 1:length(optimalLength);

    base(2*k:2*k+1) = k;
    
    if k < length(optimalLength);
        stairs(2*k-1:2*k+1) = optimalLength(k);    
    end
end

newBase = base(1:end-1);
base = newBase;
stairs(length(optimalLength)*2-1:length(optimalLength)*2) = optimalLength(end);
plot(base, stairs*1e6, 'k-')
xlabel('Layer Number')
ylabel('Widths ({\mu}m)')

axis([0, max(base), 0, 1.2*max(stairs*1e6)])
hold on

title('Layer Profile')   


% subplot(1,2,2)
% Alternative widths plot
leg = 0;
stairs = zeros(1,2*length(optimalLength));
base = stairs;
base(1) = 0;
for k = 1:length(optimalLength);

    base(2*k:2*k+1) = k;
    
    if k < length(optimalLength);
        stairs(2*k-1:2*k+1) = L(k);    
    end
end

newBase = base(1:end-1);
base = newBase;
stairs(length(optimalLength)*2-1:length(optimalLength)*2) = L(end);
plot(base, stairs*1e6, 'r--')
xlabel('Layer Number')
ylabel('Widths ({\mu}m)')

text(1, 1.1*max(optimalLength)*1e6, strcat('Total width of optimized layers: { }', num2str(1e-3*round(sum(optimalLength)*1e9)), '{\mu}m')) 
text(length(optimalLength)/2.5, 1.1*max(optimalLength)*1e6, strcat('Total width of {\lambda}/4 layers: { }', num2str(1e-3*round(sum(L)*1e9)), '{\mu}m')) 
% axis([0, max(base), 0, 1.1*max(stairs*1e6)])
legend('Optimized Profile', '{\lambda}/4 Profile')
