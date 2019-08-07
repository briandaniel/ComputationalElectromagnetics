close all

%% "ESSENTIAL" PARAMETERS
lambdaWindow = [ 800, 1050]*1e-9;
lambdaSteps = 1000;
tLim1 = 10;
tLim2 = 2;

%% STANDARD GRATING
% Grating before active region
lbdMaxRef = 974e-9;             % Wavelength of maximum reflection
grateNumAnt = 10;                % Number of bilayers in DBR
nH = 3.40;                      % High refractive index
nL = 2.90;                      % Low refractive index
lH = lbdMaxRef/(4*nH);          % Quarter wavelength widths of high refractive index
lL = lbdMaxRef/(4*nL);          % Quarter wavelength widths of low refractive index
LAnt = [repmat([lL,lH], 1, grateNumAnt), lL]; 
nAnt = [repmat([nL,nH], 1, grateNumAnt), nL];

% Grating after active region
grateNumPost = 0;                % Number of bilayers in DBR
nH = 3.40;                      % High refractive index
nL = 2.90;                      % Low refractive index
lH = lbdMaxRef/(4*nH);          % Quarter wavelength widths of high refractive index
lL = lbdMaxRef/(4*nL);          % Quarter wavelength widths of low refractive index
LPost = [repmat([lH, lL], 1, grateNumPost)]; 
nPost = [repmat([nH, nL], 1, grateNumPost)];

% Grating in active region
lbdMaxTrans = 823.6e-9;
numWell = 24;                   % Number of quantum wells
nQ = 3.7;                       % Refractive index of quantum wells
nCH = 3.40;                     % Refractive index of high spacers
nCL = 2.90;                     % Refractive index of low spacers
lWell = 8e-9;                   % Widths of quantum wells in meters (must be fixed)
lH = lbdMaxTrans/(4*nH);          % Quarter wavelength widths of high refractive index
lL = lbdMaxTrans/(4*nL);          % Quarter wavelength widths of low refractive index
sigQ = -50 ;                  % Sigmas approximating gain of quantum well
sigmaQ = repmat([sigQ, 0, 0 ], 1, numWell);   
LActive = repmat([lWell, lH, lL ], 1, numWell);
nActive = repmat([nQ, nCH, nCL],1,numWell);

% Grating in DBR
lbdMaxRef = 974e-9;             % Wavelength of maximum reflection
grateNum = 22;                  % Number of bilayers in DBR
nH = 3.40;                      % High refractive index
nL = 2.90;                      % Low refractive index
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
LInit = 100e-9;
LEnd = 100e-9;


%% OPTIMIZATION
scale = 1e9;

% Set Optimization function 'gainOptim.m'
gainFunNew = @(x) gainOptim(lbdMaxRef, n, x, sig, scale);
variables = length(L);

% Set bounds
LB = ones(size(L))*1e-9;
UB = ones(size(L))*200e-9;
% Sizes for quantum wells must remain static
for m = 1:length(L)
    if (sig(m+1) ~= 0)
        LB(m) = L(m);
        UB(m) = L(m);
    end
end

LB = LB*scale;
UB = UB*scale;

% Optimizing 
k = 0;
grateSize = sum(L);
timeElapsed = 0;
pop = 3*variables;
optionsGA = gaoptimset('TolFun', 1e-8, ...
                   'populationsize', pop, ...
                   'crossoverfraction', .5, ...
                   'initialpenalty', 10, ...
                   'display', 'iter', ...
                   'stalltimelimit', inf, ...
                   'useparallel', 'always', ...
                   'timelimit', 60*tLim1, ...
                   'generation', 2000, ...
                   'stallgenlimit', 50, ...
                   'PopInitRange', [LB; UB], ...
                   ...'hybridFcn', @patternsearch,...
                   'elitecount', round(pop/variables), ...
                   'initialpopulation', [L*scale; optimalLength]);
        
tic

optimalLength = ga (gainFunNew, variables, [], [], [], [], LB, UB, [], optionsGA );

%% Refine search
display(' ')
currentGain = -gainFunNew( optimalLength );
display(strcat('Current Optimum: ', num2str(currentGain)))

optionsPattern = psoptimset('display', 'iter', ...
                            'maxiter', 1e5 ,...
                            'timelimit', tLim2*60);
optimalLength = patternsearch(gainFunNew, optimalLength, [], [], [], [], LB, UB, [], optionsPattern );

currentGain = -gainFunNew( optimalLength );
runTime = toc;
display(' ')
display(strcat('Current Optimum: ', num2str(currentGain)))

optimizeTime = toc;
display(strcat('Optimizing time: ', num2str(optimizeTime/3600), 'hr'))

%%
close all
lambdaWindow = [ 900, 1050]*1e-9;
lambdaSteps = 10000;

currentOptimum = optimalLength./scale;

% Total index and grating vectors
nInit = 1;                      % Refractive index before structure
nEnd = 1;                       % Refractive index after structure
L = [LAnt, LActive, LPost, LDBR];
n = [nInit, nAnt, nActive, nPost, nDBR, nEnd];
sig = [0 zeros(1,length(LAnt)) sigmaQ zeros(1,length(LPost)) zeros(1,length(LDBR)), 0];

% Plot grating profile
figure('outerposition', [10, 250, 1000, 900])
% Arbitrary widths for posterior and anterior materials
subplot(2,1,2)
LInit = 100e-9;
LEnd = 100e-9;
LGraph = [LInit, L, LEnd];
currentOptimumGraph = [LInit, currentOptimum, LEnd];
base = zeros(1,length(LGraph));
baseOpt = zeros(1,length(currentOptimumGraph));

leg = 0;
legOpt = 0;
for k = 1:length(LGraph);
    base(k) = leg;
    baseOpt(k) = legOpt;
    legOpt = currentOptimumGraph(k)+legOpt;
    leg = LGraph(k)+leg;
end
hold on
base2 = [base base(end)+LEnd];
base2Opt = [baseOpt baseOpt(end)+LEnd];

stairs(base2*1e9,[n n(end)], 'k:')
stairs(base2Opt*1e9, [n n(end)]-.5, 'r-')
axis([0, max(base2)*1e9, 1, 1.1*max(n)])
text( 1e9*base(numWell), 2.02, 'Active Region', 'fontsize', 18, 'fontname', 'cordia new')
text( 1e9*base(3*numWell+grateNum), 2.02, 'DBR', 'fontsize', 18, 'fontname', 'cordia new')
annotation('doublearrow', [.15 .435], [1.09 , 1.09]/(1.1*max(n)), 'linestyle', '--', 'headstyle', 'cback1')
annotation('doublearrow', [.44 .89], [1.09 , 1.09]/(1.1*max(n)), 'linestyle', '--', 'headstyle', 'cback1')
hold off
title('Grating profile', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('nm', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Refractive Index', 'fontsize', 18, 'fontname', 'cordia new')
ledger = ['{\lambda}/4 Grating Profile'; 'Optimized grating profile  '];
hleg = legend(ledger, 'location', 'SouthEast');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')

% Generate solution based on standard active grating
[lambdas, R, T, tauPicoR, tauPicoT, phaseR, phaseT] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L, sig);
[lambdas, Ropt, Topt, tauPicoRopt, tauPicoTopt, phaseRopt, phaseTopt] = gainFun(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, currentOptimum, sig);

subplot(2,1,1)
hold on;
plot(lambdas*1e9, R, 'k--')
plot(lambdas*1e9, Ropt, 'r-')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('E_z', 'fontsize', 18, 'fontname', 'cordia new')
axis([lambdas(1)*1e9, lambdas(end).*1e9, 0, 1.2*max([max(T),max(R)])]);
ledger = ['{\lambda}/4 Reflected Spectrum'; 'Optimized Reflected Spectrum  '];
hleg = legend(ledger, 'location', 'SouthWest');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')
title('Transfer Matrix Solution Spectrum', 'fontsize', 18, 'fontname', 'cordia new')

% Show time delay and phase in another graph
figure('outerposition', [1000, 50, 800, 800])
subplot(2,1,1)
plot(lambdas*1e9, phaseR, 'k--', lambdas*1e9, phaseRopt, 'r-')
title('Phase of Spectra', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Phase (rads)', 'fontsize', 18, 'fontname', 'cordia new')
hleg = legend(ledger, 'location', 'NorthEast');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')

subplot(2,1,2)
plot(lambdas*1e9, tauPicoR, 'k--', lambdas*1e9, tauPicoRopt, 'r-')
title('Group Time Delay', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('{\lambda} [nm]', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('Time delay (ps)', 'fontsize', 18, 'fontname', 'cordia new')
legend('{\lambda}/4 Reflected Wave', '{\lambda}/4 Transmitted Wave') 
hleg = legend(ledger, 'location', 'SouthWest');
set(hleg, 'fontsize', 18, 'fontname', 'cordia new')

figure('outerposition', [950, 750, 900, 400])
stairs((currentOptimum-L).*1e9);
title('Structure adjustments made during optimization', 'fontsize', 18, 'fontname', 'cordia new')
ylabel('{\Delta}width [nm]', 'fontsize', 18, 'fontname', 'cordia new')
xlabel('Layer number', 'fontsize', 18, 'fontname', 'cordia new')

