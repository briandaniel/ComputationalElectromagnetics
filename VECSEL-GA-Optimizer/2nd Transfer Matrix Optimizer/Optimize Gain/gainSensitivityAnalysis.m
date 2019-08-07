close all
clear all

%% SETUP

% Grating parameters
load('matlab.mat')

lambdaWindow = [ lbdMaxRef*1e9-5, lbdMaxRef*1e9+5] *1e-9;
lambdaSteps = 3;

[lambdas, R, T, tauPico, tauPico2, phaseR, phaseT] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, currentOptimum, sig);
[lambdasB, RB, TB, tauPicoB, tauPico2B, phaseRB, phaseTB] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sig);

gain = R;
gainB = RB;
%% PERTURB RESULTS
close all
optimalLength = currentOptimum;

maxIndex = 2;
maxIndexOpt = 2;

obs = 300;
scale = 1:50;
percentAxis = (scale-1).*.03/max(scale-1);

avgLoss = zeros(size(scale));
avgLossOpt = zeros(size(scale));


actualGain = zeros(size(scale));
actualGainOpt = zeros(size(scale));

for m = scale;
    

    avgGain = 0;
    avgGainOpt = 0;

    averageLoss = 0;
    averageLossOpt = 0;
    
    for k = 1:obs;

        % maxIndex = find(gainB >= max(gainB)*.99,100);
        oldGain = gainB(maxIndex);
        lamVar = lambdas(maxIndex);

        % maxIndexOpt = find(gain >= max(gain)*.99, 100);
        newGain = gain(maxIndexOpt);
        lamVarOpt = lambdas(maxIndexOpt);

        % Perturbation Value (same for both but scaled appropriately)
        perc = percentAxis(m);                         % Perturbation percentage
        del = rand(1,length(L))-.5;        
        delL = perc.*mean(L).*del;
        delOpt = perc.*mean(optimalLength).*del;

        perturbOpt = optimalLength + delOpt;
        perturbL = L + delL;

        [lambdas, R, T, tauPico, tauPico2, phaseR, phaseT] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, perturbOpt, sig);
        [lambdasB, RB, TB, tauPicoB, tauPico2B, phaseRB, phaseTB] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, perturbL, sig);

        gainPerturb = R;
        gainBPerturb = RB;
        
        oldGainPerturb = gainBPerturb(maxIndex);
        newGainPerturb = gainPerturb(maxIndexOpt);
        
        %Updates
        avgGain = avgGain + oldGainPerturb;
        avgGainOpt = avgGainOpt + newGainPerturb;
        
        lossL = oldGainPerturb - oldGain;
        lossOpt = newGainPerturb - newGain;

        percLoss = lossL./oldGain;
        percLossOpt = lossOpt./newGain;

        averageLoss = averageLoss + percLoss;
        averageLossOpt = averageLossOpt + percLossOpt;
        
        
    end

    avgLoss(m) = abs(averageLoss/obs);
    avgLossOpt(m) = abs(averageLossOpt/obs);
    actualGain(m) = (avgGain/obs);
    actualGainOpt(m) = (avgGainOpt/obs); 
    m
end

%%

    figure('outerposition', [500. 200, 600, 900])
    subplot(2,1,1)
    plot(percentAxis*100, avgLoss, 'r-', percentAxis*100, avgLossOpt, 'b-')
    legend('{\lambda}/4', 'Optimized')
    ylabel('% Gain loss at 650nm')
    xlabel('% Structure Perturbation')
    
    title('Sensitivity of optimized total gain (reflected+transmitted) active DBR')
    text(percentAxis(round(scale(end/10)))*100, avgLoss(end), strcat('Number of random observations:{ }', num2str(obs)))
    
    subplot(2,1,2)
    plot(percentAxis*100, actualGain, 'r-', percentAxis*100, actualGainOpt, 'b-', ...
         percentAxis*100, repmat(gainB(maxIndex),1, length(percentAxis)), 'r:', ...
         percentAxis*100, repmat(gain(maxIndex),1, length(percentAxis)), 'b:')
    legend('{\lambda}/4', 'Optimized',  'Max gain {\lambda}/4', 'Max gain optimized')
    ylabel('% Average Gain')
    xlabel('% Structure Perturbation')
   
    return
%%

% RUN WITH ONLY SINGLE LAYER
% Single layer is the same size as total optimized structure
close all

% Perturb results
% L2 = optimalLength;
% L2 = L2+mean(L2).*(rand(1,length(L2))-.5).*.02;
% optimalLength = L2;
tic
lambdaWindow = [ 500, 800]*1e-9;
lambdaSteps = 1000;

[gain, lambdas, R, T, R2, T2, tauPicoR, tauPicoT] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, optimalLength, sigma, sigma2);
toc
[gainB, lambdasB, RB, TB, R2B, T2B, tauPicoRB, tauPicoTB] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, perturbL, sigma, sigma2);


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
