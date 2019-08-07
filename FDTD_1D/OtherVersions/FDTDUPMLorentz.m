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
spaceSize = 4.5e-6;            % Space size (1 = 1meter)                  [USER DEFINED]
gridSize = 1000;               % Total number of cells                    [USER DEFINED]
dx = spaceSize/gridSize;       % Cell Size (1 = 1meter)                   

% Define FDTD temporal space
maxIter = 2000;          % Maximum number of iterations + 1         [USER DEFINED]
dt = .5*dx/(cNot);              % Time step (1 = 1sec)                     
Sc = cNot*dt/dx;                % Courant Number (only works in Vacuum)
stoppingNorm = .005;

% Print current settings
display(strcat('Spacial domain =  ',  num2str(round(spaceSize*1e7)*1e-1), ' microns'))
display(strcat('dx =  ',  num2str(round(dx*1e10)*1e-1), ' nm'))
display(strcat('Time domain =  ',  num2str(round(maxIter*dt*1e16)*1e-1), ' fs'))
display(strcat('Time step =  ',  num2str(round(dt*1e19)*1e-4), ' fs'))

% Define important positions
widthPML = 50;                 % Width of PML in grid points              [USER DEFINED]
sourcePos = .5e-6 ;              % Location of source in meters             [USER DEFINED]
sourcePos = round(sourcePos/dx);
posMed = 1.5e-6;                  % Position of medium relative to source    [USER DEFINED]
                                % in meters
posMed = round(posMed/dx)+sourcePos;
sizeMed = round(1.5e-6/dx);
% Grating parameters
lambdaMax = 500e-9;
nH = 1.82;
nL = 1.38;
lH = lambdaMax/(4*nH);
lL = lambdaMax/(4*nL);
grateNum = 4;                   % Number of grating sets
grateWidth = lH*grateNum + lL*grateNum;

% grateWidth = 2e-6;             % Width of entire grating
% relativeSize = lH/lL;          % Relative widths of the alternating grates

grateEnd = posMed+sizeMed;
grateSize1 = round(lH/dx);
grateSize2 = round(lL/dx); 




%% Initialize solution vectors
xDom = 1:gridSize;
Ex = zeros(1,gridSize);
Dx = Ex;
DxStore = Dx;
Hy = zeros(1,gridSize-1);
By = Hy;
ByStore = By;
dField = Dx;
ExStoreEnd = Ex(end-1);
ExStoreStart = Ex(2);

%% Hard source parameters

% Gaussian Pulse parameters
width = 50;                    % Width of Gaussian                        [USER DEFINED]
maxT = width*5;

% Sine Wave Parameters
f = 650e12;                     % Frequency                                 [USER DEFINED]
amp = 1;                        % Amplitude of Wave                        [USER DEFINED]
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
kappaPML2 = 1+(kappaMax-1).*(xVec).^3

% plot((xVec+.5)./widthPML,sigmaPMLx, 'ro', (xVec./widthPML), sigmaPMLx2 ,'b*')

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

dEps = 3;                       % Change in rel. Permittivity due to pole  [USER DEFINED]
epsInf = 2;                    % Rel. Permittivity at infinite freq.      [USER DEFINED]
freqL = 500e2;                   % Undamped Resonance Frequency of med      [USER DEFINED]         
omegaMed = 2*pi*freqL;          % Angular freq equivalent of freqL
deltaMed = 0*omegaMed  ;       % Conductivity of Lorentz Medium           [USER DEFINED]

% Parameters associated with lorentz medium
alphaL = deltaMed;
betaL = sqrt(omegaMed^2-deltaMed^2);
gamma = dEps*omegaMed^2/betaL

chiNot = -j*gamma/(alphaL-j*betaL)*(1-exp((-alphaL+j*betaL)*dt));     % Initial X, susceptibility
xiNot = (-j*gamma/dt)/(alphaL-j*betaL)^2*...
       (((alphaL-j*betaL)*dt+1)*exp((-alphaL+j*betaL)*dt)-1)  ;     % Inital Zeta, Convolution term
   
dXnot = chiNot*(1-exp((-alphaL+j*betaL)*dt));
dZnot = xiNot*(1-exp((-alphaL+j*betaL)*dt));

chiNot = real(chiNot);
xiNot = real(xiNot);

psi = Ex;
ExStore = Ex;
alpha = (epsInf-xiNot)/(epsInf-xiNot+chiNot);
beta = 1/(epsInf-xiNot+chiNot);
kappa = 1/(epsInf-xiNot+chiNot);

c1a = c1;
c2a = c2;

medVecAlpha = ones(1,gridSize);
medVecBeta = ones(1,gridSize);

% simple dielectric equations
medVecAlpha(posMed:grateEnd) = alpha*ones(1,grateEnd-posMed+1);
medVecBeta(posMed:grateEnd) = beta*ones(1,grateEnd-posMed+1);

c1a(widthPML+1:gridSize -widthPML-1) = medVecAlpha(widthPML+1:gridSize -widthPML-1);
c2a(widthPML+1:gridSize -widthPML-1) = medVecBeta(widthPML+1:gridSize -widthPML-1);


% Generate plots of constants of lorentz medium
figure('outerposition',[1 200 500 500])
fstep = 10^3;
fdom = linspace(omegaMed/10,omegaMed*10,fstep);
chi = (dEps*omegaMed^2)./(omegaMed^2 + 2*j*deltaMed.*fdom - fdom.^2);
epsOmega = epsInf*ones(1,fstep) + chi;
semilogx(fdom/2/pi, epsOmega, omegaMed/2/pi, 2, 'ro')
axis([omegaMed/10/2/pi, omegaMed*10/2/pi, min(epsOmega), max(epsOmega)])
title('Dielectric Constant')
xlabel('Frequency')

%% Initialization of Detector for frequency domain

% Specify detector locations
detInitial = sourcePos+1;             % Position of initial wave detector  [USER DEFINED]
detPost = posMed+1;                % Position of latter wave detector  [USER DEFINED]
detRef = posMed-150;                 % Position of reflection detector

detectorInitial = zeros(1,maxIter-1);
detectorFinal = zeros(1,maxIter-1);
detectorRef = detectorFinal;
detectordField = detectorFinal;
%%

%%%%%% MAIN SOLUTION LOOP %%%%%%
figure('outerposition', [710   500   1000   700])
set(1,'color',[.9 .9 .9])

for n = 1:maxIter

% Generate Gaussian Pulse
% if n < maxT*2;
% source = exp(-.5*((maxT*dx-n*dt*cNot/Sc)/(width*dx))^2);
% end

% Generate wave packet
source = exp(-.5*((maxT*dx-n*dt*cNot/Sc)/(width*dx))^2) * amp*sin(omega*n*dt);
    
% Generate sine hard source
% source = amp*sin(omega*n*dt);

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

% E-field Boundary Update
Ex(1)   = 0 ;
Ex(end) = 0 ;

% E-field Update
DxStore = Dx;
Dx(2:end-1) = c1a(2:end-1).*Dx(2:end-1) + Sc.*c2a(2:end-1).*(Hy(2:end) - Hy(1:end-1));
Ex = c3.*Ex + c4.*(c5.*Dx - c6.*DxStore) + kappa*real(psi);
dField(2:end-1) = dField(2:end-1) + Sc.*(Hy(2:end) - Hy(1:end-1));

% Hard source update
if n < maxT*2;
Ex(sourcePos)    = source;
end

% Recursive accumlator for Lorentz medium
psi(posMed:grateEnd) = (dXnot-dZnot)*Ex(posMed:grateEnd) + dZnot*ExStore(posMed:grateEnd) + psi(posMed:grateEnd)*exp((-alphaL+j*betaL)*dt);
ExStore = Ex;

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
if n <= maxIter-2;
    if n < 850
        detectorInitial(n) = Ex(detInitial);
    end

    if n < 2000
        detectorFinal(n) = Ex(detPost);
        detectordField(n) = dField(detPost);
    end

    if n > 850
        detectorRef(n) = Ex(detRef);
    end
end
%%
if mod(n,16) == 0;
if n < maxIter-2    
% Generate plots/movie
hold on
disConv = dx*10^6;

% Plot grating rectangles
rectangle('position',[(posMed)*disConv, -2, (grateEnd-posMed)*disConv, 4], 'facecolor', [0 .5 .7])

plot(xDom*disConv, Ex, 'r-')
axis([0 max(xDom*disConv) -1.2 1.2])
title('Wave packet entering perfect dielectric')
xlabel('{\mu}m')
ylabel('E_z')

plot(detInitial*disConv, 0, 'ro', detPost*disConv, 0, 'go' , detRef*disConv, 0, 'b+')
end

timeStep = num2str(n);
time = strcat('Time Step: { }', timeStep);
text(gridSize/1.5*disConv, .5, time);

courant = num2str(Sc);
infoC = strcat('Courant number {\nu}= ', courant);
text(gridSize/1.5*disConv, .3, infoC); 

timeElaps = num2str(round(n*dt*10^16)*10^-1);
timeElapsed = strcat('Time elapsed:{ } ', timeElaps, 'fs');
text(gridSize/1.5*disConv, .15, timeElapsed); 


M(n/16) = getframe(gcf);
end
if n < maxIter-1
    hold off
    clf
end

if n == maxIter-2 | (n > 500 && norm(Ex) < stoppingNorm)
    
close(2)
    %%
    maxIter = length(detectorInitial);
    
    figure('outerposition', [1000   500   1000   500])

    % Plot detected waves
    subplot(2,1,1)
    observations = length(detectorInitial);
    timeDom = (0:observations-1).*dt;
    plot(timeDom, detectorInitial, timeDom, detectorFinal, timeDom, detectorRef)
    legend('Initial wavelet', 'Transmitted wavelet', 'Reflected wavelet')
    title('Detected waves in time domain')
    xlabel('t, time')
    ylabel('Amplitude')
        
    dfield = detectordField;
    
    % Fourier analysis with interpolation of recorded time domain waves
    scaledDom = timeDom*10^15;
    
    ppInc = spline(scaledDom, detectorInitial);
    ppTrans = spline(scaledDom, detectorFinal);
    ppRef = spline(scaledDom, detectorRef);
    ppdField = spline(scaledDom, dfield);
    
    maxT = max(scaledDom);
    samples = maxIter*10;
    samples = 2^nextpow2(samples);
    samples2 = maxIter;
    Fs = 1/dt;
    
    freqDom = (0:1:round(samples/2)-1)*Fs/samples/(maxIter)*samples2;
    fourierInit = fft(ppval(ppInc, linspace(0,maxT,samples2)), samples);
    fourierFinal = fft(ppval(ppTrans, linspace(0,maxT,samples2)), samples);
    fourierRef = fft(ppval(ppRef, linspace(0,maxT,samples2)), samples);
    fourierdField = fft(ppval(ppdField, linspace(0,maxT,samples2)), samples);
    
    magInit = (abs(fourierInit(1:round(samples/2)))).^2;
    magFinal = (abs(fourierFinal(1:round(samples/2)))).^2;
    magRef = (abs(fourierRef(1:round(samples/2)))).^2;
    magdField = (abs(fourierdField(1:round(samples/2)))).^2;
    
    scaledFreqDom = freqDom.*10^-12;
    [A indexMax] = max(magInit);
    primaryFreq = freqDom(indexMax);
    wavelengthPrimary = cNot/primaryFreq;
    pointsPerWave = wavelengthPrimary/dx;
    
    % Plot freq. spectrum
    subplot(2,1,2)
    hold on
    plot(freqDom./10^12, magInit, 'k')
    plot(freqDom./10^12, magFinal, 'b')
    plot(freqDom./10^12, magRef, 'g')
    hold off
    ylabel('Power') 
    xlabel('Frequency (THz)')
    legend('Inital Spectrum', 'Transmitted Spectrum', 'Reflected spectrum')
    title('Frequency Spectrum')
    axis([300 1200 0 1.2*max(magInit)])
    
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
    
    epsApprox = magdField./(magFinal);
    
    
    trans = sqrt(trans);
    refl = sqrt(refl);
    hold on
    plot(freqDom./10^12, trans, freqDom./10^12, refl)
    ylabel('Transmitted/Reflected') 
    xlabel('Frequency (THz)')
    legend('Transmitted(E^t/E^i)', 'Reflectance(E^r/E^i)')
    title('Transmission/Reflection')
    
peakFreq = num2str(round(primaryFreq/10^11)*10^-1);
ppLambda = num2str(round(pointsPerWave));
waveLeg = num2str(round(wavelengthPrimary*10^10)*10^-1);
info1 = strcat('Points per wavelength = { } ', ppLambda, '{ } at peak frequency:{ }', peakFreq, 'THz');
info2 = strcat('Primary wavelength = { }', waveLeg, '{ }nm');

    axis([300 1200 0 1.2])
    text(150, .8, info1)
    text(150, .6, info2)
    
% Wavelength domain    
    figure('outerposition', [1600   50   1000   500])
    lambdaDom = cNot./freqDom;
    plot(lambdaDom*10^9, trans, lambdaDom*10^9, refl)
    xlabel('{\lambda} (nm)')
    ylabel('Transmitted/Reflected') 
    legend('Transmittance (E^t/E^i)', 'Reflectance(E^r/E^i)', 'Sum')
    title('Transmission/Reflection')    
    axis([300 800 0 1.2])
%%
 


% % Generate plots of constants of lorentz medium
figure('outerposition',[1 200 500 500])
fstep = 10^4;
fdom = linspace(100e12,1200e12,fstep);
chi = (dEps*omegaMed^2)./(omegaMed^2 + 2*j*deltaMed.*fdom - fdom.^2);
epsOmega = epsInf*ones(1,fstep) + chi;
plot(fdom/2/pi, epsOmega, omegaMed/2/pi, 2, 'ro', freqDom, epsApprox)
axis([200e12/(2*pi), 1500e12/(2*pi), 0, 4])
title('Dielectric Constant')
xlabel('Frequency')


    return
end


end
    
%% Call brag function
close all

lambdaWindow = [ 300, 800];
lambdaSteps = 500;

ppR = pchip(lambdaDom*1e9, refl.^2);
compare = ppval(ppR, linspace(lambdaWindow(1), lambdaWindow(2), lambdaSteps));

[lambdas, R] = multiLayer(506.5, lambdaWindow(1), lambdaWindow(2), lambdaSteps, nH, nL, grateNum);
figure('outerposition', [100   500   800   500])

R = R./100;
error = compare-R;

plot(lambdas, R, lambdas, compare, 'r-')
hold on;
plot(lambdas, error+.8, 'g-.');

%  plot(lambdas.*10^9, TP)
legend('Reflection calculated', 'Reflection from simulation', 'Error')
ylabel('Power')
xlabel('Wavelength(nm)')
axis([lambdaWindow(1) lambdaWindow(2) 0 1.1])

%% 
%movie2avi(M, 'Gauss Pulse 1D Debye Media', 'fps', 40);

