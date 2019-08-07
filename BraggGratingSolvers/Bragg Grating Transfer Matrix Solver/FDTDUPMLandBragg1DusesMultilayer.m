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
spaceSize = 5.5e-6;            % Space size (1 = 1meter)                  [USER DEFINED]
gridSize = 500;               % Total number of cells                    [USER DEFINED]
dx = spaceSize/gridSize;       % Cell Size (1 = 1meter)                   

% Define FDTD temporal space
maxIter = gridSize*100;          % Maximum number of iterations + 1         [USER DEFINED]
dt = .5*dx/(cNot);              % Time step (1 = 1sec)                     
Sc = cNot*dt/dx;                % Courant Number (only works in Vacuum)
stoppingNorm = .005;

% Print current settings
display(strcat('Spacial domain =  ',  num2str(round(spaceSize*1e7)*1e-1), ' microns'))
display(strcat('dx =  ',  num2str(round(dx*1e10)*1e-1), ' nm'))
display(strcat('Time domain =  ',  num2str(round(maxIter*dt*1e16)*1e-1), ' fs'))
display(strcat('Time step =  ',  num2str(round(dt*1e19)*1e-4), ' fs'))

% Define important positions
widthPML = 20;                 % Width of PML in grid points              [USER DEFINED]
sourcePos = .5e-6 ;              % Location of source in meters             [USER DEFINED]
sourcePos = round(sourcePos/dx);
posMed = 1e-6;                  % Position of medium relative to source    [USER DEFINED]
                                % in meters
posMed = round(posMed/dx)+sourcePos;

% Grating parameters
lambdaMax = 500e-9;
nH = 2.10;
nL = 1.59;
lH = lambdaMax/(4*nH);
lL = lambdaMax/(4*nL);
grateNum = 8;                   % Number of grating sets
grateWidth = lH*grateNum + lL*grateNum;

% grateWidth = 2e-6;             % Width of entire grating
% relativeSize = lH/lL;          % Relative widths of the alternating grates

grateEnd = posMed+round(grateWidth/dx);
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

ExStoreEnd = Ex(end-1);
ExStoreStart = Ex(2);

%% Hard source parameters

% Gaussian Pulse parameters
width = 50;                    % Width of Gaussian                        [USER DEFINED]
maxT = width*5;

% Sine Wave Parameters
f = 590e12;                     % Frequency                                 [USER DEFINED]
amp = 1;                        % Amplitude of Wave                        [USER DEFINED]
lambda =cNot/f;                 % Wavelength
omega = 2*pi*f;                 % Angular Frequency

% Ricker Wavelet Parameters
fp = 5e8;                       % Peak Frequency                           [USER DEFINED]
md = 1;                         % Temporal Delay Multiple                  [USER DEFINED]
dr = md/fp;                     % Temporal Delay

%% Parameters associated with the dielectric 
% First dielectric
sigma = 0;                      % Conductivity of dielectric               [USER DEFINED]
dielC = nH^2;                   % Dielectric Constant (rel. Permittivity)  [USER DEFINED]    
epsDiel = dielC*epsNot;         % Permittivity of dielectric
factor = dt*sigma/(2*epsNot*dielC);
alpha = (1-factor)/(1+factor);
beta = 1/(dielC*(1+factor));

% Second dielectric
sigma2 = 0;                    % Conductivity of dielectric               [USER DEFINED]
dielC2 = nL^2;                 % Dielectric Constant (rel. Permittivity)  [USER DEFINED]    
epsDiel2 = dielC2*epsNot;      % Permittivity of dielectric

factor2 = dt*sigma2/(2*epsNot*dielC2);
alpha2 = (1-factor2)/(1+factor2);
beta2 = 1/(dielC2*(1+factor2));

dn = sqrt(dielC2) - sqrt(dielC)
bandwidth = .99*dn/pi*620e-9*1e9

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

c1a = c1;
c2a = c2;
medVecAlpha = ones(1,gridSize);
medVecBeta = ones(1,gridSize);

cPos = posMed;
for k = 1:grateNum;
    
    medVecAlpha(cPos:cPos + grateSize1-1) = alpha*ones(1,grateSize1);
    medVecBeta(cPos:cPos + grateSize1-1) = beta*ones(1,grateSize1);
    cPos = cPos + grateSize1;
    
    medVecAlpha(cPos:cPos + grateSize2-1) = alpha2*ones(1,grateSize2);
    medVecBeta(cPos:cPos + grateSize2-1) = beta2*ones(1,grateSize2);
        
    cPos = cPos + grateSize2;
end

medVecAlpha(cPos:cPos + grateSize1-1) = alpha*ones(1,grateSize1);
medVecBeta(cPos:cPos + grateSize1-1) = beta*ones(1,grateSize1);
cPos = cPos + grateSize1;

% simple dielectric equations
% medVecAlpha(posMed:grateEnd) = alpha2*ones(1,grateEnd-posMed+1);
% medVecBeta(posMed:grateEnd) = beta2*ones(1,grateEnd-posMed+1);

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
detPost = grateEnd +15;                % Position of latter wave detector  [USER DEFINED]
detRef = detInitial+10;                 % Position of reflection detector

detectorInitial = zeros(1,maxIter-1);
detectorFinal = zeros(1,maxIter-1);
detectorRef = detectorFinal;

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


% E-field Update
DxStore = Dx;
Dx(2:end-1) = c1a(2:end-1).*Dx(2:end-1) + Sc.*c2a(2:end-1).*(Hy(2:end) - Hy(1:end-1));
Ex = c3.*Ex+c4.*(c5.*Dx - c6.*DxStore);

% E-field Boundary Update
Ex(1)   = 0 ;
Ex(end) = 0 ;

% Hard source update
if n < maxT*2;
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
if n <= maxIter-2;
    if n < 500
        detectorInitial(n) = Ex(detInitial);
    end

    if n > 1
        detectorFinal(n) = Ex(detPost);
    end

    if n > 500
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
cPos = posMed;
for k = 1:grateNum
    rectangle('position',[cPos*disConv, -2, grateSize1*disConv, 10], 'facecolor', [0 .5 .7])
    cPos = cPos + grateSize1;
    rectangle('position',[cPos*disConv, -2, grateSize2*disConv, 10], 'facecolor', [0 0  .9])
    cPos = cPos + grateSize2;
end

rectangle('position',[cPos*disConv, -2, grateSize1*disConv, 10], 'facecolor', [0 .5 .7])
cPos = cPos + grateSize1;
    
plot(xDom*disConv, Ex, 'r-')
axis([0 max(xDom*disConv) -1.2 1.2])
title('Wave packet entering perfect dielectric')
xlabel('{\mu}m')
ylabel('E_z')

plot(detInitial*disConv, 0, 'ro', detPost*disConv, 0, 'go' , detRef*disConv, 0, 'b+')
end

% Display pertinent parameters
% waveFreq = num2str(f/10^6);
% wave = strcat('Peak Frequency: ', waveFreq, 'MHz');
% text(gridSize/10, 1.05, wave)
% 
timeStep = num2str(n);
time = strcat('Time Step: { }', timeStep);
text(gridSize/1.5*disConv, .5, time);

courant = num2str(Sc);
infoC = strcat('Courant number {\nu}= ', courant);
text(gridSize/1.5*disConv, .3, infoC); 

timeElaps = num2str(round(n*dt*10^16)*10^-1);
timeElapsed = strcat('Time elapsed:{ } ', timeElaps, 'fs');
text(gridSize/1.5*disConv, .15, timeElapsed); 

% dielEps = num2str(dEps);
% dielectricCons = strcat('{\Delta} {\epsilon}_r = { } ', dielEps);
% text(posMed+5, .05, dielectricCons)
% 
% dielEps = num2str(epsInf);
% dielectricCons = strcat('{\Delta} {\epsilon}_{\infty} = { } ', dielEps);
% text(posMed+5, .1, dielectricCons)

% dielCons = num2str(dielC);
% dielectricCons = strcat('Dielectric Constant {\epsilon}_r = ', dielCons);
% text(posMed+1, 1.0, dielectricCons)
% 
% conduct = num2str(sigma);
% conductivity = strcat('Conductivity {\sigma} = ', conduct);
% text(posMed+1, .8, conductivity)
% 
% subplot(2,1,2)
% hold on
% rectangle('position',[posMed, -1.2, gridSize, 4], 'facecolor', [.8 .8 .8])
% plot(xDom(1:end-1),Hy)
% axis([0 gridSize -1.2 1.2])  
% xlabel('x/4 (mm)')
% ylabel('H_y')
%  
M(n/16) = getframe(gcf);
end
if n < maxIter-1
    hold off
    clf
end

if n == maxIter-2 | (n > 500 && norm(Ex) < stoppingNorm)
    
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
    samples = maxIter*100;
    samples = 2^nextpow2(samples);
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
    plot(freqDom./10^12, refl.^2)
    ylabel('Transmitted/Reflected') 
    xlabel('Frequency (THz)')
    legend('Reflectance(E^r/E^i)')
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
    plot(lambdaDom*10^9, refl.^2)
    xlabel('{\lambda} (nm)')
    ylabel('Transmitted/Reflected') 
    legend('Transmittance (E^t/E^i)', 'Reflectance(E^r/E^i)', 'Sum')
    title('Transmission/Reflection')    
    axis([300 800 0 1.2])

    return
end


end
    
%% Call brag function

lambdaWindow = [ 300, 800]*1e-9;
lambdaSteps = 500;
correction = 0;

ppR = pchip(lambdaDom+ correction, refl.^2);
compare = ppval(ppR, linspace(lambdaWindow(1), lambdaWindow(2), lambdaSteps));

lL = grateSize2*dx;
lH = grateSize1*dx;
N = grateNum;
n = [1, nH, repmat([nL,nH], 1, N), 1];        % indices for the layers A|H(LH)N |G 
L = [lH, repmat([lL,lH], 1, N)];                % Lengths of layers

[lambdas, R, r, T, t] = transferMatrix(lambdaWindow(1), lambdaWindow(2), lambdaSteps, n, L);

fq = cNot./lambdas;
freqDomNew = linspace(100e12,2000e12, 1e5)
ppr = pchip(fq, r);
rOmega = ppval(ppr,freqDomNew);

figure
plot(freqDomNew./1e12, abs(rOmega).^2)

thetaR = phase(r);
thetaT = phase(t);
dthetaR = gradient(theta);
dthetaT = gradient(theta);


% 
% delay3 = -((lambdas.^2./(2*pi*cNot))).*dtheta1;
% delay = delay3*1e12;
% plot(delay)

error = compare'-R;
figure('outerposition', [100   500   800   500])

plot(lambdas*1e9, R, lambdas*1e9, T, lambdas*1e9, compare, lambdaDom*1e9+correction*1e9, trans.^2)
axis([lambdaWindow(1)*1e9 lambdaWindow(2)*1e9 0 1.1])

hold on;
plot(lambdas*1e9, error+.5, 'g-.', lambdas*1e9, .5*ones(length(lambdas),1), 'k-.');
return
%  plot(lambdas.*10^9, TP)
legend('Reflection calculated', 'Reflection from simulation', 'Error')
ylabel('Power')
xlabel('Wavelength(nm)')



%% 
%movie2avi(M, 'Gauss Pulse 1D Debye Media', 'fps', 40);

