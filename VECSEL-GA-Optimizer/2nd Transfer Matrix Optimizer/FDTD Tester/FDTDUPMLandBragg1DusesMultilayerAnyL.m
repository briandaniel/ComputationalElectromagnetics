clear all;
close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Info                                                            %
% Simulates FDTD in one dimension                                         %
%                                                                         %
% Note: Electric field, E, undergoes a change of variables for            %
% simplification: E_prgrm = sqrt(eps_not/mu_not)E_exact                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physical Constants
epsNot = 8.854187817620e-12;    % Permittivity of Free Space
muNot = 1.25663706e-6;          % Magnetic Constant
cNot = 2.99792458e8;            % Speed of light
etaNot = sqrt(muNot/epsNot);    % Characteristic impedence of free space

%% Parameters

% Define FDTD physical space
spaceSize = 11.2e-6;            % Space size (1 = 1meter)                  [USER DEFINED]
gridSize = 8000;               % Total number of cells                    [USER DEFINED]
dx = spaceSize/gridSize;       % Cell Size (1 = 1meter)                   
stopTime = 3200;
% Define FDTD temporal space
maxIter = gridSize*250;          % Maximum number of iterations + 1         [USER DEFINED]
dt = .5*dx/(cNot);              % Time step (1 = 1sec)                     
Sc = cNot*dt/dx;                % Courant Number (only works in Vacuum)
stoppingNorm = .0002;
disConv = dx*10^6;
skip = 20;
start = maxIter/2;

% sVid = videowriter('newVid.mp4')
% sVid.FrameRate = 15;
% sVid.Quality = 99;
% open(sVid)
% Print current settings
display(strcat('Spacial domain =  ',  num2str(round(spaceSize*1e7)*1e-1), ' microns'))
display(strcat('dx =  ',  num2str(round(dx*1e10)*1e-1), ' nm'))
display(strcat('Time domain =  ',  num2str(round(maxIter*dt*1e16)*1e-1), ' fs'))
display(strcat('Time step =  ',  num2str(round(dt*1e19)*1e-4), ' fs'))

% Define important positions
widthPML = 50;                 % Width of PML in grid points              [USER DEFINED]
sourcePos = 1e-6 ;              % Location of source in meters             [USER DEFINED]
sourcePos = round(sourcePos/dx);
posMed = 1.2e-6;                  % Position of medium relative to source    [USER DEFINED]
                                % in meters
posMed = round(posMed/dx)+sourcePos;

% Grating parameters
load newData3.mat
grateSize1 = round(lH/dx);
grateSize2 = round(lL/dx); 
refractiveIndices = n;
L = currentOptimum;
positions = round(L./dx);
grateEnd = posMed + sum(positions);



%% Initialize solution vectors
xDom = 1:gridSize;
Ex = zeros(1,gridSize);
Dx = Ex;
DxStore = Dx;
Hy = zeros(1,gridSize-1);
By = Hy;
ByStore = By;

ExStoreEnd = Ex(end-1);
ExStoreStart = Ex(2);

%% Hard source parameters

% Gaussian Pulse parameters
width = 320;                    % Width of Gaussian                        [USER DEFINED]
maxT = width*5;

% Sine Wave Parameters
f = 308e12;                     % Frequency                                 [USER DEFINED]
amp = 5;                        % Amplitude of Wave                        [USER DEFINED]
lambda =cNot/f;                 % Wavelength
omega = 2*pi*f;                 % Angular Frequency

% Ricker Wavelet Parameters
fp = 5e8;                       % Peak Frequency                           [USER DEFINED]
md = 1;                         % Temporal Delay Multiple                  [USER DEFINED]
dr = md/fp;                     % Temporal Delay

%% Parameters associated with the dielectric 
%UPML Parameters
xVec = widthPML:-1:0;
sigMax = .12*3.5/(etaNot*dx);
kappaMax = 1;
sigmaPML = ((xVec+.5)./widthPML).^3*sigMax;
sigmaPML2 = (xVec./widthPML).^3*sigMax;
kappaPML = 1+(kappaMax-1).*(xVec+.5).^3;
kappaPML2 = 1+(kappaMax-1).*(xVec).^3;


sigmaPML = [ sigmaPML, zeros(1,gridSize-2*widthPML-2), sigmaPML(length(sigmaPML):-1:1)];
kappaPML = [ kappaPML, ones(1,gridSize-2*widthPML-2), kappaPML(length(kappaPML):-1:1)];
sigmaPML2 = [ sigmaPML2, zeros(1,gridSize-2*widthPML-3), sigmaPML2(length(sigmaPML2):-1:1)];
kappaPML2 = [ kappaPML2, ones(1,gridSize-2*widthPML-3), kappaPML2(length(kappaPML2):-1:1)];

c1 = (2.*epsNot.*kappaPML-sigmaPML.*dt)./(2.*epsNot.*kappaPML+sigmaPML.*dt);
c2 = (2.*epsNot)./(2.*epsNot.*kappaPML+sigmaPML.*dt);
c3 = (2.*epsNot.*kappaPML-sigmaPML.*dt)./(2.*epsNot.*kappaPML+sigmaPML.*dt);
c4 = 1./((2.*epsNot.*kappaPML+sigmaPML.*dt));
c5 = (2.*epsNot.*kappaPML+sigmaPML.*dt);
c6 = (2.*epsNot.*kappaPML+sigmaPML.*dt);

d1 = (2.*epsNot.*kappaPML2-sigmaPML2.*dt)./(2.*epsNot.*kappaPML2+sigmaPML2.*dt);
d2 = (2.*epsNot)./(2.*epsNot.*kappaPML2+sigmaPML2.*dt);
d3 = (2.*epsNot.*kappaPML2-sigmaPML2.*dt)./(2.*epsNot.*kappaPML2+sigmaPML2.*dt);
d4 = 1./((2.*epsNot.*kappaPML2+sigmaPML2.*dt));
d5 = (2.*epsNot.*kappaPML2+sigmaPML2.*dt);
d6 = (2.*epsNot.*kappaPML2+sigmaPML2.*dt);

c1a = c1;
c2a = c2;
medVecAlpha = ones(1,gridSize);
medVecBeta = ones(1,gridSize);

cPos = posMed;
for k = 1:length(positions);
    sigma = sig(k+1);                  
    dielC = refractiveIndices(k+1)^2;     
    factor = dt*sigma/(2*epsNot*dielC);
    alpha = (1-factor)/(1+factor);
    beta = 1/(dielC*(1+factor));

    medVecAlpha(cPos:cPos + positions(k)-1) = alpha*ones(1,positions(k));
    medVecBeta(cPos:cPos + positions(k)-1) = beta*ones(1,positions(k));
    cPos = cPos + positions(k);
end

c1a(widthPML+1:gridSize -widthPML-1) = medVecAlpha(widthPML+1:gridSize -widthPML-1);
c2a(widthPML+1:gridSize -widthPML-1) = medVecBeta(widthPML+1:gridSize -widthPML-1);
%% Parameters associated with Debye Material
% dEps = 3;                       % Change in rel. Permittivity due to pole  [USER DEFINED]
% epsInf = 1;                     % Rel. Permittivity at infinite freq.      [USER DEFINED]  
% tau  = .001e-6  ;               % Pole relaxation time                     [USER DEFINED]
% chiNot = dEps*(1-exp(-dt/tau));                    % Initial X, susceptibility
% xiNot = dEps*tau/dt*(1-(dt/tau+1)*exp(-dt/tau));  % Inital Zeta, Convolution term
% dXnot = chiNot*(1-exp(-dt/tau));
% dZnot = xiNot*(1-exp(-dt/tau));
% psi = Ex(posMed:end-1);
% ExStore = Ex;
% alpha = (epsInf-xiNot)/(epsInf-xiNot+chiNot);
% beta = 1/(epsInf-xiNot+chiNot);
% kappa = 1/(epsInf-xiNot+chiNot);


%% Parameters associated with Lorentz Material

% dEps = 3;                       % Change in rel. Permittivity due to pole  [USER DEFINED]
% epsInf = 1.5;                   % Rel. Permittivity at infinite freq.      [USER DEFINED]
% freqL = .5e9;                   % Undamped Resonance Frequency of med      [USER DEFINED]         
% omegaMed = 2*pi*freqL;          % Angular freq equivalent of freqL
% sigmaMed = .1*omegaMed  ;       % Conductivity of Lorentz Medium           [USER DEFINED]
% 
% % Parameters associated with lorentz medium
% alphaL = sigmaMed;
% betaL = sqrt(omegaMed^2-sigmaMed^2);
% gamma = dEps*omegaMed^2/betaL
% 
% chiNot = -j*gamma/(alphaL-j*betaL)*(1-exp((-alphaL+j*betaL)*dt));     % Initial X, susceptibility
% xiNot = (-j*gamma/dt)/(alphaL-j*betaL)^2*...
%        (((alphaL-j*betaL)*dt+1)*exp((-alphaL+j*betaL)*dt)-1)  ;     % Inital Zeta, Convolution term
%    
% dXnot = chiNot*(1-exp((-alphaL+j*betaL)*dt))
% dZnot = xiNot*(1-exp((-alphaL+j*betaL)*dt));
% 
% chiNot = real(chiNot)
% xiNot = real(xiNot)
% 
% psi = Ex(posMed:end-1);
% ExStore = Ex;
% alpha = (epsInf-xiNot)/(epsInf-xiNot+chiNot);
% beta = 1/(epsInf-xiNot+chiNot);
% kappa = 1/(epsInf-xiNot+chiNot);


%% Initialization of Detector for frequency domain

% Specify detector locations
detInitial = sourcePos+1;             % Position of initial wave detector  [USER DEFINED]
% detPost = posMed+100;                % Position of latter wave detector  [USER DEFINED]
detPost = grateEnd+20;                % Position of latter wave detector  [USER DEFINED]
detRef = detInitial+10;                 % Position of reflection detector

detectorInitial = zeros(1,maxIter-1);
detectorFinal = zeros(1,maxIter-1);
detectorRef = detectorFinal;

%%

%%%%%% MAIN SOLUTION LOOP %%%%%%
figure('outerposition', [210   200   1300   700])
set(1,'color',[.9 .9 .9])

for n = 1:maxIter

% Generate Gaussian Pulse
% if n < maxT*2;
% source = exp(-.5*((maxT*dx-n*dt*cNot/Sc)/(width*dx))^2);
% end

% Generate wave packet
source = exp(-.5*((maxT*dx-n*dt*cNot/Sc)/(width*dx))^2) * amp*sin(omega*n*dt);
    
% Generate sine hard source
% if n < maxT*2;
% source = amp*sin(omega*n*dt);
% end
% Generate Ricker Wavelet
% source = (1-2*(pi*fp*(n*dt-dr))^2)*exp(-(pi*fp*(n*dt-dr))^2);


%% Without dielectric

% 
% % E-field Update
% Ex(2:end-1) = Ex(2:end-1) - Sc*( Hy(2:end) - Hy(1:end-1) );
% Ex(sourcePos)  = source;
% 
% % MUR boundary
% Ex(end) = ExStoreEnd + ((cNot*dt-dx)/(cNot*dt+dx))*(Ex(end-1)-Ex(end));
% Ex(1) = ExStoreStart + ((cNot*dt-dx)/(cNot*dt+dx))*(Ex(2)-Ex(1));
% 
% % Store Ex
% ExStoreEnd = Ex(end-1);
% ExStoreStart = Ex(2);
% 
% % H-field Update
% Hy          =      Hy     - Sc*( Ex(2:end) - Ex(1:end-1) ); 

%% With Dielectric and UPML


% E-field Update
DxStore = Dx;
Dx(2:end-1) = c1a(2:end-1).*Dx(2:end-1) + Sc.*c2a(2:end-1).*(Hy(2:end) - Hy(1:end-1));
Ex = c3.*Ex+c4.*(c5.*Dx - c6.*DxStore);

% E-field Boundary Update
Ex(1)   = 0 ;
Ex(end) = 0 ;

% Hard source update
if n < stopTime;
Ex(sourcePos)    = source;
end

% H-field Update
% Hy = Hy - Sc.*( Ex(2:end)  - Ex (1:end-1) + psiH(2:end)); 
ByStore = By;

By = d1.*By + Sc.*d2.*(Ex(2:end) - Ex(1:end-1));


Hy = d3.*Hy+d4.*(d5.*By - d6.*ByStore);


%% With Debye Media
% % E-field Boundary Update
% Ex(1)   = 0 ;
% Ex(end) = 0 ;
% 
% % E-field Update
% % Before media is reached
% Ex(2:posMed-1)   = Ex(2:posMed-1) - Sc*( Hy(2:posMed-1) - Hy(1:posMed-2) );
% % Update after media is reached
% Ex(posMed:end-1) = alpha*Ex(posMed:end-1) - ...
%                    beta*Sc*( Hy(posMed:end) - Hy(posMed-1:end-1) ) + ...
%                    kappa*psi;
% Ex(sourcePos)    = source;
% 
% % Update Psi
% psi = (dXnot-dZnot)*Ex(posMed:end-1) + dZnot*ExStore(posMed:end-1) + psi*exp(-dt/tau);
% ExStore = Ex;
% 
% 
% % H-field Update
% Hy  =  Hy - Sc*( Ex(2:end) -  Ex(1:end-1) ); 


%% With Lorentz Media
% E-field Boundary Update
% Ex(1)   = 0 ;
% Ex(end) = 0 ;
% 
% % E-field Update
% % Before media is reached
% Ex(2:posMed-1)   = Ex(2:posMed-1) - Sc*( Hy(2:posMed-1) - Hy(1:posMed-2) );
% % Update after media is reached
% Ex(posMed:end-1) = alpha*Ex(posMed:end-1) - ...
%                    beta*Sc*( Hy(posMed:end) - Hy(posMed-1:end-1) ) + ...
%                    kappa*real(psi);
% Ex(sourcePos)    = source;
% 
% % Update Psi
% psi = (dXnot-dZnot)*Ex(posMed:end-1) + dZnot*ExStore(posMed:end-1) + psi*exp((-alphaL+j*betaL)*dt);
% ExStore = Ex;
% 
% % H-field Update
% Hy  =  Hy - Sc*( Ex(2:end) -  Ex(1:end-1) ); 



%% Update detector

if n < stopTime
    detectorInitial(n) = Ex(detInitial);
end

if n > 1
    detectorFinal(n) = Ex(detPost);
end

if n > stopTime
    detectorRef(n) = Ex(detRef);
end

%%
if mod(n,skip) == 0 && n >= 1 ;
if n < maxIter-2  
% Generate plots/movie
hold on

% Plot grating rectangles
cPos = posMed;
maxN = max(refractiveIndices);

for k = 1:length(positions);
    if sig(k+1) == 0;
    rectangle('position',[cPos*disConv, -2, positions(k)*disConv, 10],...
    'facecolor', [refractiveIndices(k+1) refractiveIndices(k+1) refractiveIndices(k+1)]/maxN, ...
    'edgecolor', 'none')
    else
    rectangle('position',[cPos*disConv, -2, positions(k)*disConv, 10],...
    'facecolor', [1 0 0], 'edgecolor', [1 0 0], 'linewidth', 2)
    end
    cPos = cPos + positions(k);
end

plot(xDom*disConv, Ex, 'b-', 'linewidth', 2)
axis([0 max(xDom*disConv) -.2 .2])
title('Wave packet entering perfect dielectric')
xlabel('{\mu}m')
ylabel('E_z')

plot(detInitial*disConv, 0, 'ro', detPost*disConv, 0, 'go' , detRef*disConv, 0, 'b+')

end

timeStep = num2str(n);
time = strcat('Time Step: { }', timeStep);
text(gridSize/1.5*disConv, .2, time);

courant = num2str(Sc);
infoC = strcat('Courant number {\nu}= ', courant);
text(gridSize/1.5*disConv, .3, infoC); 

timeElaps = num2str(round(n*dt*10^16)*10^-1);
timeElapsed = strcat('Time elapsed:{ } ', timeElaps, 'fs');
text(gridSize/1.5*disConv, .15, timeElapsed); 
s = getframe;

if n > start
writeVideo(sVid, s) 
display('Writing Video')
end

end
if n < maxIter-1
    hold off
    clf
end

if n == maxIter-2 | (n > 500 && abs(norm(Ex(widthPML:end-widthPML), inf)) < stoppingNorm)    
close(1)
    %%
    maxIter = length(detectorInitial);
    
    figure('outerposition', [1000   800   1000   800])
    
    % Plot detected waves
    subplot(2,1,2)
    observations = length(detectorInitial);
    timeDom = (0:observations-1).*dt;
    plot(timeDom, detectorInitial, timeDom, detectorFinal, timeDom, detectorRef)
    legend('Initial wavelet', 'Transmitted wavelet', 'Reflected wavelet')
    title('Detected waves in time domain')
    xlabel('t, time')
    ylabel('Amplitude')
    
    
    % Fourier analysis with interpolation of recorded time domain waves
    scaledDom = timeDom*10^15;
    ppInc = spline(scaledDom, detectorInitial);
    ppTrans = spline(scaledDom, detectorFinal);
    ppRef = spline(scaledDom, detectorRef);
    
    maxT = max(scaledDom);
    samples = maxIter*1;
    samples = round(maxIter)
    samples = 2^nextpow2(samples*10);
    samples2 = maxIter;
    Fs = 1/dt;

    
    freqDom = (0:1:round(samples/2)-1)*Fs/samples/(maxIter)*samples2;
    fourierInit = fft(ppval(ppInc, linspace(0,maxT,samples2)), samples);
    fourierFinal = fft(ppval(ppTrans, linspace(0,maxT,samples2)), samples);
    fourierRef = fft(ppval(ppRef, linspace(0,maxT,samples2)), samples);
    magInit = (abs(fourierInit(1:round(samples/2)))).^2;
    magFinal = (abs(fourierFinal(1:round(samples/2)))).^2;
    magRef = (abs(fourierRef(1:round(samples/2)))).^2;
    
    scaledFreqDom = freqDom.*10^-12;
    [A indexMax] = max(magInit);
    primaryFreq = freqDom(indexMax);
    wavelengthPrimary = cNot/primaryFreq;
    pointsPerWave = wavelengthPrimary/dx;
    
    % Plot freq. spectrum
    subplot(2,1,1)
    hold on
    plot(freqDom./10^12, magInit, 'y')
    plot(freqDom./10^12, magFinal, 'g')
    plot(freqDom./10^12, magRef, 'r')
    plot(freqDom./10^12, magRef+magFinal, 'k:')
    hold off
    ylabel('Power') 
    xlabel('Frequency (THz)')
    legend('Inital Spectrum', 'Transmitted Spectrum', 'Reflected spectrum', 'Sum of transmiited and reflected spectra')
    title('Frequency Spectrum')
    axis([200 1400 0 1.2*max(magInit)])
    
    % Plot transmission graph
    figure('outerposition', [1400   250   1000   500])
    trans = zeros(length(freqDom),1);
    refl = trans;
    
    for k = 1:round(samples/2);
        if magInit(k) > 1e-12
             trans(k) = magFinal(k)/magInit(k);
             refl(k) = magRef(k)/magInit(k);
        end
    end
    
    trans = sqrt(trans);
    refl = sqrt(refl);
    hold on
    plot(freqDom./10^12, refl, freqDom./10^12, trans)
    ylabel('Transmitted/Reflected') 
    xlabel('Frequency (THz)')
    legend('Reflectance(E^r/E^i)', 'Transmittance(E^t/E^i)')
    title('Transmission/Reflection')
    
peakFreq = num2str(round(primaryFreq/10^11)*10^-1);
ppLambda = num2str(round(pointsPerWave));
waveLeg = num2str(round(wavelengthPrimary*10^10)*10^-1);
info1 = strcat('Points per wavelength = { } ', ppLambda, '{ } at peak frequency:{ }', peakFreq, 'THz');
info2 = strcat('Primary wavelength = { }', waveLeg, '{ }nm');

    axis([1 2400 0 1.2])
    text(150, .8, info1)
    text(150, .6, info2)
    
% Wavelength domain    
    figure('outerposition', [1600   50   1000   500])
    lambdaDom = cNot./freqDom;
    plot(lambdaDom*10^9, trans.^2, lambdaDom*10^9, refl.^2, lambdaDom*10^9, refl.^2+trans.^2)
    xlabel('{\lambda} (nm)')
    ylabel('Transmitted/Reflected') 
    legend('Transmittance (E^t/E^i)', 'Reflectance(E^r/E^i)', 'Sum')
    title('Transmission/Reflection')    
    axis([300 1800 0 1.2])

%% Gain Plot    
    figure('outerposition', [50   50   600   600])
    finalSpec = magFinal+magRef;
    gain = (finalSpec./magInit);
   
    plot(lambdaDom*10^9, gain.^1-1)
    xlabel('{\lambda} (nm)')
    ylabel('Gain') 
    legend('E^f/E^i - 1')
    title('Optimized for gain at 650nm')    
    axis([500 800 0 2000])
    
    return
end


end
    
%% Call brag function

lambdaWindow = [ 920, 1050]*1e-9;
lambdaSteps = 10000;
correction = 0;
ppR = pchip(lambdaDom+ correction, refl.^2);
compareR = ppval(ppR, linspace(lambdaWindow(1), lambdaWindow(2), lambdaSteps));
ppT = pchip(lambdaDom+ correction, trans.^2);
compareT = ppval(ppT, linspace(lambdaWindow(1), lambdaWindow(2), lambdaSteps));
newL = positions*dx;
[lambdas, R, r, T, t] = transferMatrix(lambdaWindow(1), lambdaWindow(2), lambdaSteps, refractiveIndices, newL, sig);
[lambdas, R2, r2, T2, t2] = transferMatrix(lambdaWindow(1), lambdaWindow(2), lambdaSteps, refractiveIndices, currentOptimum, sig);

%%
% error = compare'-R;
figure('outerposition', [100   500   800   500])
hold on
plot(lambdas*1e9, R, 'r-', lambdas*1e9, T, 'b-')
plot(lambdas*1e9, R2, 'r--', lambdas*1e9, T2, 'b--')
plot(lambdas*1e9, compareR, 'g--')
plot(lambdas*1e9, compareT, 'k')
axis([lambdaWindow(1)*1e9 lambdaWindow(2)*1e9 0 2.1])

hold on
% plot(lambdas*1e9, error+.5, 'g-.', lambdas*1e9, .5*ones(length(lambdas),1), 'k-.');
%  plot(lambdas.*10^9, TP)
% legend('Reflection calculated', 'Transmittance calculated', 'Reflection from simulation', 'Transmission from simulation', 'Error')
ylabel('Power')
xlabel('Wavelength(nm)')



%% 
close (sVid);
