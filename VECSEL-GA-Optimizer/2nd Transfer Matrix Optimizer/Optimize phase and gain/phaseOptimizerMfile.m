close all

%% "ESSENTIAL" PARAMETERS
lambdaWindowOpt = [ 972, 976]*1e-9;
lambdaStepsOpt = 16;
optNum = inf;                   % Number of loops that the optimizer conducts
optReset = 500;                   % Numer of loops before optimizer resets
tolerance = 1;

lambdaWindow = [ 900, 1020]*1e-9;
lambdaSteps = 1000;
%% STANDARD GRATING
% Grating before active region
lbdMaxRef = 974e-9;             % Wavelength of maximum reflection
grateNumAnt = 6;                % Number of bilayers in DBR
nH = 3.40;                      % High refractive index
nL = 3.10;                      % Low refractive index
lH = lbdMaxRef/(4*nH);          % Quarter wavelength widths of high refractive index
lL = lbdMaxRef/(4*nL);          % Quarter wavelength widths of low refractive index
LAnt = [repmat([lL,lH], 1, grateNumAnt), lL]; 
nAnt = [repmat([nL,nH], 1, grateNumAnt), nL];

% Grating after active region
grateNumPost = 6;                % Number of bilayers in DBR
nH = 3.40;                      % High refractive index
nL = 3.10;                      % Low refractive index
lH = lbdMaxRef/(4*nH);          % Quarter wavelength widths of high refractive index
lL = lbdMaxRef/(4*nL);          % Quarter wavelength widths of low refractive index
LPost = [repmat([lH, lL], 1, grateNumPost)]; 
nPost = [repmat([nH, nL], 1, grateNumPost)];

% Grating in active region
numWell = 14;                   % Number of quantum wells
nQ = 3.7;                       % Refractive index of quantum wells
nCH = 3.40;                     % Refractive index of high spacers
nCL = 3.10;                     % Refractive index of low spacers
lWell = 8e-9;                   % Widths of quantum wells in meters (must be fixed)
lCladH = lH;                    % Widths of surrounding cladding (High refractive index)
lCladL = lL;                    % Widths of surrounding cladding (Low refractive index)
sigQ = -300 ;                  % Sigmas approximating gain of quantum well
sigmaQ = repmat([sigQ, 0, 0 ], 1, numWell);   
LActive = repmat([lWell, lCladH, lCladL ], 1, numWell);
nActive = repmat([nQ, nCH, nCL],1,numWell);

% Grating in DBR
lbdMaxRef = 974e-9;             % Wavelength of maximum reflection
grateNum = 22;                  % Number of bilayers in DBR
nH = 3.40;                      % High refractive index
nL = 3.10;                      % Low refractive index
lH = lbdMaxRef/(4*nH);          % Quarter wavelength widths of high refractive index
lL = lbdMaxRef/(4*nL);          % Quarter wavelength widths of low refractive index
LDBR = [lH, repmat([lL,lH], 1, grateNum)]; 
nDBR = [nH, repmat([nL,nH], 1, grateNum)];

% Total index and grating vectors
nInit = 1;                      % Refractive index before structure
nEnd = 1;                       % Refractive index after structure
L = [LAnt, LActive, LPost, LDBR];
n = [nInit, nAnt, nActive, nPost, nDBR, nEnd];
sig = [0 zeros(1,length(LAnt)) sigmaQ zeros(1,length(LPost)) zeros(1,length(LDBR)), 0];
% Arbitrary widths for posterior and anterior materials
LInit = 100e-9;
LEnd = 100e-9;

%{
% Plot grating profile
figure('outerposition', [10, 250, 1000, 900])
subplot(2,1,2)

LGraph = [LInit, L, LEnd];
base = zeros(1,length(LGraph));
leg = 0;
for k = 1:length(LGraph);
    base(k) = leg;
    leg = LGraph(k)+leg;
end
hold on
base2 = [base base(end)+LEnd];
stairs(base2*1e9,[n n(end)], 'k-')
axis([0, max(base2)*1e9, 1, 1.1*max(n)])
text( 1e9*base(numWell), 2.02, 'Active Region', 'fontsize', 18, 'fontname', 'cordia new')
text( 1e9*base(3*numWell+grateNum), 2.02, 'DBR', 'fontsize', 18, 'fontname', 'cordia new')
annotation('doublearrow', [.15 .435], [1.09 , 1.09]/(1.1*max(n)), 'linestyle', '--', 'headstyle', 'cback1')
annotation('doublearrow', [.44 .89], [1.09 , 1.09]/(1.1*max(n)), 'linestyle', '--', 'headstyle', 'cback1')
hold off
title('Grating profile', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('nm', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Refractive Index', 'fontsize', 18, 'fontname', 'cordia new')

% Generate solution based on standard active grating
[lambdas, R, T, tauPicoR, tauPicoT, phaseR, phaseT] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sig);
subplot(2,1,1)
hold on;
plot(lambdas*1e9, R, 'r-')
plot(lambdas*1e9, T, 'b-')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('E_z', 'fontsize', 18, 'fontname', 'cordia new')
axis([lambdas(1)*1e9, lambdas(end).*1e9, 0, 1.2*max([max(T),max(R)])]);
ledger = ['{\lambda}/4 Reflected Wave  '; '{\lambda}/4 Transmitted Wave'];
hleg = legend(ledger, 'location', 'SouthWest');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')
title('Transfer Matrix Solution Spectrum', 'fontsize', 18, 'fontname', 'cordia new')

% Show time delay and phase in another graph
figure('outerposition', [1000, 50, 800, 800])
subplot(2,1,1)
plot(lambdas*1e9, phaseR, 'r-', lambdas*1e9, phaseT, 'b-')
title('Phase of Spectra', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Phase (rads)', 'fontsize', 18, 'fontname', 'cordia new')
hleg = legend(ledger, 'location', 'NorthEast');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')

subplot(2,1,2)
plot(lambdas*1e9, tauPicoR, 'r-', lambdas*1e9, tauPicoT, 'b-')
title('Group Time Delay', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Time delay (ps)', 'fontsize', 18, 'fontname', 'cordia new')
legend('{\lambda}/4 Reflected Wave', '{\lambda}/4 Transmitted Wave') 
hleg = legend(ledger, 'location', 'SouthWest');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')
%}

% Load previous data
fileName = input('Enter filename.mat for data storage: ', 's')
save(fileName)

%% OPTIMIZATION
% Assign values if the optimizing functions do not exist
resetVar = input('Reset variables? Type "0" if no, and "1" if yes.');

if exist('currentGain', 'var') == 0 || resetVar == 1
    currentGain = 1;
end

if exist('currentOptimum', 'var') == 0 || resetVar == 1
    currentOptimum = L;
end

if exist('storage', 'var') == 0 || resetVar == 1
    storage = L;
end

if exist('storageBin', 'var') == 0 || resetVar == 1
    storageBin = L;
end

if exist('storageGain', 'var') == 0 || resetVar == 1
    storageGain = 0;
end

% Set Optimization function 'gainOptim.m'
gainFunNew = @(x) phaseOptim(lambdaWindowOpt(1), lambdaWindowOpt(2), lambdaStepsOpt, n, x, sig);
gainFunNewOne = @(x) phaseOptimOne(lambdaWindowOpt(1), lambdaWindowOpt(2), lambdaStepsOpt, n, x, sig, tolerance);

variables = length(L);
options = gaoptimset('tolfun', 1e-4, 'populationsize', 20, 'timelimit', inf, ...
                     'display', 'off', 'stalltimelimit', inf, 'useparallel', 'always', 'crossoverfraction', .9, ...
                     'initialPopulation', storageBin);
% Set bounds

LB = ones(size(L))*.5*mean(L);
UB = ones(size(L))*2*mean(L);
% Sizes for quantum wells must remain static
for m = 1:length(L)
    if (sig(m+1) ~= 0)
        LB(m) = L(m);
        UB(m) = L(m);
    end
end

% Optimizing loop
k = 0;
grateSize = sum(L);
timeElapsed = 0;
tic

while k <= optNum
    tic
        k = k+1;
        optimalLength = gamultiobj (gainFunNew, variables, [], [], [], [], LB, UB, options );
%           optimalLength = ga(gainFunNewOne, variables, [], [], [], [], LB, UB, [], options );
    % Store previous solutions
    optimalLength = optimalLength(end,:);
    
    [newGainInv, newPhase] = gainFunNew(optimalLength)
    newGain = 1/newGainInv;
    newOptimum = optimalLength;

    
%     objective = gainFunNewOne(newOptimum);
    if newGainInv + newPhase < currentGain;
        currentGain = newGainInv + newPhase;
        currentGain2 = 1/newGainInv;
        currentPhase = newPhase;
        currentOptimum = newOptimum;
        storage = [storage ; currentOptimum];

    end

    if length(storage(:,1)) >= 10;
        storage = currentOptimum;
    end


    if k >= optReset;
        storageBin = [storageBin; storage];
        storageGain = [storageGain; currentGain];
        storage = zeros(1,length(currentOptimum));
        currentOptimum = storage;
        currentGain = 0;
        k = 0;
        save(fileName) 
    end
    
    runTime = toc;
    timeElapsed = timeElapsed + runTime;
    display(' ')
    display(strcat('Iterations: ', num2str(k)))
    display(strcat('Current Optimum: ', num2str(currentGain2)))
    display(strcat('Current Phase Curvature Norm: ', num2str(currentPhase)))
    display(strcat('Time elapsed: ', num2str(timeElapsed), '      Run: ', num2str(runTime)))
end

optimizeTime = toc;
display(strcat('Optimizing time: ', num2str(optimizeTime/3600), 'hr'))

%%
% Plot grating profile
figure('outerposition', [10, 250, 1000, 900])
% Arbitrary widths for posterior and anterior materials
subplot(2,1,2)
LGraph = [LInit, currentOptimum, LEnd];
base = zeros(1,length(LGraph));
leg = 0;
for k = 1:length(LGraph);
    base(k) = leg;
    leg = LGraph(k)+leg;
end
hold on
base2 = [base base(end)+LEnd];
stairs(base2*1e9,[n n(end)], 'k-')
axis([0, max(base2)*1e9, 0, 1.1*max(n)])
text( 1e9*base(numWell), 2.02, 'Active Region', 'fontsize', 18, 'fontname', 'cordia new')
text( 1e9*base(3*numWell+grateNum), 2.02, 'DBR', 'fontsize', 18, 'fontname', 'cordia new')
annotation('doublearrow', [.15 .435], [1.09 , 1.09]/(1.1*max(n)), 'linestyle', '--', 'headstyle', 'cback1')
annotation('doublearrow', [.44 .89], [1.09 , 1.09]/(1.1*max(n)), 'linestyle', '--', 'headstyle', 'cback1')
hold off
title('Grating profile', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('nm', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Refractive Index', 'fontsize', 18, 'fontname', 'cordia new')

% Generate solution based on standard active grating
[lambdas, R, T, tauPicoR, tauPicoT, phaseR, phaseT] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, currentOptimum, sig);
subplot(2,1,1)
hold on;
plot(lambdas*1e9, R, 'r-')
plot(lambdas*1e9, T, 'b-')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('E_z', 'fontsize', 18, 'fontname', 'cordia new')
axis([lambdas(1)*1e9, lambdas(end).*1e9, 0, 1.2*max([max(T),max(R)])]);
ledger = ['{\lambda}/4 Reflected Wave  '; '{\lambda}/4 Transmitted Wave'];
hleg = legend(ledger, 'location', 'SouthWest');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')
title('Transfer Matrix Solution Spectrum', 'fontsize', 18, 'fontname', 'cordia new')

% Show time delay and phase in another graph
figure('outerposition', [1000, 50, 800, 800])
subplot(2,1,1)
plot(lambdas*1e9, phaseR, 'r-', lambdas*1e9, phaseT, 'b-')
title('Phase of Spectra', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Phase (rads)', 'fontsize', 18, 'fontname', 'cordia new')
hleg = legend(ledger, 'location', 'NorthEast');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')

subplot(2,1,2)
plot(lambdas*1e9, tauPicoR, 'r-', lambdas*1e9, tauPicoT, 'b-')
title('Group Time Delay', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Time delay (ps)', 'fontsize', 18, 'fontname', 'cordia new')
legend('{\lambda}/4 Reflected Wave', '{\lambda}/4 Transmitted Wave') 
hleg = legend(ledger, 'location', 'SouthWest');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')

