close all
clear all
%% "ESSENTIAL" PARAMETERS
lambdaWindow = [ 800, 1050]*1e-9;
lambdaSteps = 1000;
optNum = inf;                   % Number of loops that the optimizer conducts
optReset = 2000;                   % Numer of loops before optimizer resets

%% STANDARD GRATING

% Grating in active region
lbdMaxRef = 974e-9;             % Wavelength of maximum reflection
numWell = 10;                    % Number of quantum wells
nQ = 3.7;                       % Refractive index of quantum wells
nCH = 3.40;                     % Refractive index of high spacers
nCL = 3.10;                     % Refractive index of low spacers
lWell = 8e-9;                   % Widths of quantum wells in meters (must be fixed)
lH = lbdMaxRef/(4*nCH);          % Quarter wavelength widths of high refractive index
lL = lbdMaxRef/(4*nCL);          % Quarter wavelength widths of low refractive index
lCladH = lH;                    % Widths of surrounding cladding (High refractive index)
lCladL = lL;                    % Widths of surrounding cladding (Low refractive index)
sigQ = -500 ;                   % Sigmas approximating gain of quantum well
sigmaQ = repmat([0, sigQ, 0 ], 1, numWell);   
LActive = repmat([ lCladH , lWell, lCladL], 1, numWell);
nActive = repmat([ nCH, nQ, nCL],1,numWell);


% Total index and grating vectors
nInit = 1;                      % Refractive index before structure
nEnd = 1;                       % Refractive index after structure
L = [ LActive ];
n = [nInit, nActive, nEnd];
sig = [0 sigmaQ 0];
LInit = 100e-9;
LEnd = 100e-9;

%{
% Plot grating profile
figure('outerposition', [10, 250, 1000, 900])
% Arbitrary widths for posterior and anterior materials
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

%% OPTIMIZATION
scaleFactor = 1e8;

% Set Optimization function 'gainOptim.m'
gainFunNew = @(x) gainOptimL(lbdMaxRef, n, x, lWell, numWell, sig, scaleFactor);
variables = 2;

% Set bounds
LB = ones(size(L))*1e-9;
UB = ones(size(L))*200e-9;

values = linspace(min(LB)*scaleFactor, max(UB)*scaleFactor, 200);
[X Y] = meshgrid(values);

space = size(X);
Z = zeros(space);

for k = 1:space(1);
    for m = 1:space(2);
        Z(k,m) = -gainFunNew([X(k,m), Y(k,m)]);
    end
end



% Sizes for quantum wells must remain static
for m = 1:length(L)
    if (sig(m+1) ~= 0)
        LB(m) = L(m);
        UB(m) = L(m);
    end
end

%%
% Optimizing 


global pop;
global storage
global storageInit
storage = [];
storageInit = [];
pop = 300;

k = 0;
grateSize = sum(L);
timeElapsed = 0;
optionsGA = gaoptimset('TolFun', 1e-16, ...
                   'populationsize', pop, ...
                   'crossoverfraction', .5, ...
                   'initialpenalty', 10, ...
                   'display', 'iter', ...
                   'stalltimelimit', inf, ...
                   'useparallel', 'always', ...
                   'timelimit', inf, ...
                   'generation', 500, ...
                   'stallgenlimit', 50, ...
                   'elitecount', 0);

optionsGA = gaoptimset('TolFun', 1e-4, ...
                   'tolcon', 1e-4, ...
                   'populationsize', pop, ...
                   'crossoverfraction', .3, ...
                   'initialpenalty', 10, ...
                   'display', 'iter', ...
                   'stalltimelimit', inf, ...
                   'useparallel', 'always', ...
                   'timelimit', inf, ...
                   'generation', 500, ...
                   'stallgenlimit', 50, ...
                   'PopInitRange', [min(LB); max(UB)]*scaleFactor, ...
                   'elitecount', 2);         
tic

optimalLength = ga (gainFunNew, variables, [], [], [], [], [min(LB), min(LB)]*scaleFactor, [max(UB), max(UB)]*scaleFactor, [], optionsGA );

%%
display(' ')
currentGain = -gainFunNew( optimalLength );
display(strcat('Current Optimum: ', num2str(currentGain)))

optionsSIMPLEX = optimset('tolX', 1e-16, ...
                          'tolfun', 1e-16, ...
                          'tolcon', 1e-16, ...
                          'maxfunevals', 1e6, ...
                          'maxIter', 1e6);
% optimalLength = fmincon(gainFunNew,optimalLength, [],[],[],[],LB, UB,[], optionsSIMPLEX);
optimalLength = fminsearch(gainFunNew,optimalLength,optionsSIMPLEX);

currentGain = -gainFunNew( optimalLength );
runTime = toc;
display(' ')
display(strcat('Current Optimum: ', num2str(currentGain)))

optimizeTime = toc;
display(strcat('Optimizing time: ', num2str(optimizeTime/3600), 'hr'))

%%
figure('outerposition', [ 100, 100, 1700, 1000])
hold on

mesh(X/scaleFactor, Y/scaleFactor, Z);
plot3(optimalLength(1)/scaleFactor, optimalLength(2)/scaleFactor, -gainFunNew(optimalLength), 'g.', 'markersize', 40);

Zmax = max(max(Z))*ones(length(storage),1);
Zmax2 = max(max(Z))*ones(length(storageInit),1);
plot3(storage(:,1), storage(:,2), Zmax, 'k.')
plot3(storageInit(:,1), storageInit(:,2), Zmax2, 'b.')
return


%%
close all

lambdaWindow = [ 800, 1050]*1e-9;
lambdaSteps = 10000;

% Plot grating profile
figure('outerposition', [10, 250, 1000, 900])
% Arbitrary widths for posterior and anterior materials
subplot(2,1,2)
LGraph = [LInit, optimalLength, LEnd];
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
% text( 1e9*base(3*numWell+grateNum), 2.02, 'DBR', 'fontsize', 18, 'fontname', 'cordia new')
annotation('doublearrow', [.15 .435], [1.09 , 1.09]/(1.1*max(n)), 'linestyle', '--', 'headstyle', 'cback1')
annotation('doublearrow', [.44 .89], [1.09 , 1.09]/(1.1*max(n)), 'linestyle', '--', 'headstyle', 'cback1')
hold off
title('Grating profile', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('nm', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Refractive Index', 'fontsize', 18, 'fontname', 'cordia new')

% Generate solution based on standard active grating
[lambdas, R, T, tauPicoR, tauPicoT, phaseR, phaseT] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, optimalLength, sig);
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

