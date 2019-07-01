clear all;
close all;

% Physical Constants
epsNot = 8.854187817620e-12;        % Permittivity of Free Space
muNot  = 1.25663706e-6;             % Magnetic Constant
cNot   = 2.99792458e8;              % Speed of light in vacuum
etaNot = sqrt(muNot/epsNot);        % Characteristic impedence of free space
planck = 6.62606957e-34;            % Plancks constant
dirac  = planck/(2*pi);             % Dirac constant (Reduced Plancks constant)
e = 1.60217646e-19;                 % Charge of an electron (Coulombs)
me = 9.10938188e-31;                % Electron Mass (kg)

%% Program Info: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Simulates wave propagation in 1D using the FDTD method                  %
% 1. Utilizes Ex and Hy fields propagating in z direction                 %
% 2. CPML utilized to truncate spatial domain                             %
% 3. Source condition can be varied in the settings                       %
% 4. Coefficient vectors allow simple implementation of refractive        %
%    indices and/or simple lossy media                                    %
% 5. Model interacts with an inhomogeneous medium constructed using       %
%    a set of 2-level atoms models represented by Maxwell-Bloch equations %
%                                                                         %
%                                                                         %
%  This Source Code Form is subject to the terms of the Mozilla Public    %
%  License, v. 2.0. If a copy of the MPL was not distributed with this    %
%  file, You can obtain one at http://mozilla.org/MPL/2.0/.               %
%                                                                         %
%  Copyright (C) 2012-2013 Brian D. Hong                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



%% Parameters and Settings
%%%%%%%%%%%%%%%%%%%%%%% Define FDTD solution space %%%%%%%%%%%%%%%%%%%%%%%% 
spaceSize    = 23.0e-6;         % Spactial Domain  (1 = 1 meter)           [USER DEFINED]
timeDom      = 200e-15;         % Time domain      (1 = 1 sec)             [USER DEFINED]
gridSize     = 600;            % Total number of spatial cells            [USER DEFINED]
dx           = spaceSize/gridSize;                 
dt           = .7*dx/(cNot);  
maxIter      = round(timeDom/dt);    
skipFrames   = 300;

%%%%%%%%%%%%%%%%%%%%%%%% Material Model Parameters %%%%%%%%%%%%%%%%%%%%%%%%
locM  = [2 22]*1e-6;            % Starting and ending location of medium
T1    = 100e-12;                % Excited state lifetime
T2    = 3e-12;                  % Dephasing time
TNatom = 10e24;                 % Number of atoms per cubic meter
feLimits = [100e12 300e12];     % Atomic transition resonance frequency limits
numFreq = 301;                  % Number of transition frequencies
numFreq2 = 60;
feLimits2 = [100e12 250e12];
fe = linspace(feLimits(1), feLimits(2), numFreq)';
u = linspace(-numFreq/2,numFreq/2,numFreq); 
v = 1.3*(exp(-(u).^2/50^2))'; 
Natom = v.*TNatom/numFreq;
gammaN = 1e-29;
rho30 = -1;                     % Initial excitation state in medium
normalizer = 1;
normalizer3 = 1;
amp = 4.2186e9;                 % Amplitude of Wave                        [USER DEFINED]

%%%%%%%%%%%%%%%%%%%%%%%% Hard source parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian = 1,        Modulated Gaussian = 2,  Sine Wave = 3, 
% Ricker Wavelet = 4;  SECH pulse = 5;          pump = 6;

pulseSetting = 5;               % Choose which source condition to use     [USER DEFINED]
sourcePos    = 1e-6 ;           % Location of source in meters             [USER DEFINED]
f            = 200e12;          % Frequency                                [USER DEFINED]
sourcePos    = round(sourcePos/dx);

% Gaussian Pulse parameters
width        = 10e-15;         % Approximate width of pulse in time (s)   [USER DEFINED]
maxT         = width;
% Sine Wave Parameters
lambda       = cNot/f;          % Wavelength
omega        = 2*pi*f;          % Angular Frequency
% Ricker Wavelet Parameters
fp           = 5e8;             % Peak Frequency                           [USER DEFINED]
md           = 1;               % Temporal Delay Multiple                  [USER DEFINED]
dr           = md/fp;           % Temporal Delay
% SECH pulse Parameters
tpS           = 20/f;            % Approximate width of pulse in time (s)   [USER DEFINED]
% Pump pulse parameters
tp           = 5e-15;
td           = 2*tp;
betaPump     = 1e-4;
startUp      = 5;               % Cycles before startup

%%%%%%%%%%%%%%%%%%%%%%% Frequency domain parameters %%%%%%%%%%%%%%%%%%%%%%%
% lambdaWindow = [ 100, 300]*1e-9;
freqMin = 100e12;
freqMax = 300e12;
stepsF       = 1000;
detInitial   = sourcePos+5;     % Position of initial detector             [USER DEFINED]
detPost      = round(locM(2)/dx)+2;    % Position of latter detector              [USER DEFINED]
detRef       = detInitial+10;   % Position of reflection detector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PML Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmlWidth     = 20;              % Width of PML in grid points              [USER DEFINED]
sigM         = 3;               % Power of sigma in PML                    [USER DEFINED]
sigMax       = 5*(sigM+1)*...   % Max value of sigma in PML
              .8/(etaNot*dx);
kappaMax     = 1;               % Max value of kappa in PML                [USER DEFINED]
kappaM       = 3;               % Power of kappa in PML                    [USER DEFINED]
aMax         = 1;               % Max value of kappa in PML                [USER DEFINED]
aM           = 3;               % Power of a in PML                        [USER DEFINED]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Other Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stopTime = round(tpS/dt);
posMed = 3e-6;                  % Position of medium relative to source    [USER DEFINED]                              
posMed = round(posMed/dx)+sourcePos;
visible = [0 14.0 -1 1];
visible2 = [0 14 -.2 .2];
stoppingNorm = .0002;

% Print current settings
display(strcat('Spacial domain =  ',  num2str(round(spaceSize*1e7)*1e-1), ' microns'))
display(strcat('dx =  ',  num2str(round(dx*1e10)*1e-1), ' nm'))
display(strcat('Time domain =  ',  num2str(round(maxIter*dt*1e16)*1e-1), ' fs'))
display(strcat('Time step =  ',  num2str(round(dt*1e19)*1e-4), ' fs'))


%% Frequency domain parameters
% freqMin = cNot/lambdaWindow(1);
% freqMax = cNot/lambdaWindow(2);
testFreq = linspace(freqMin, freqMax, stepsF);
testWL = cNot./testFreq;
EfInit  = zeros(size(testFreq));
EfFinal = zeros(size(testFreq));
EfRef   = zeros(size(testFreq));

%% Initialize solution vectors
xDom = 1:gridSize;

% Electric Field
Ex = zeros(1,gridSize);
Jx = zeros(1,gridSize);

% Magnetic Field
Hy = zeros(1,gridSize-1);
My = zeros(1,gridSize-1);

% Relative permittivity/permeability
sigx = zeros(1,gridSize);
epsx = ones(1,gridSize);
muy = ones(1,gridSize-1);

% PML Auxilliary variables
QEx = zeros(size(Hy));
QHy = zeros(size(Ex));

%% Compute coefficient vectors
% E-Field coeff.
cx1 = (1-(sigx.*dt)./(2*epsNot.*epsx))./(1+(sigx.*dt)./(2*epsNot.*epsx));
cx2 = (dt./(epsNot.*epsx))./(1+(sigx.*dt)./(2*epsNot.*epsx));

% D-Field coeff.
dy1 = dt./(muNot*muy);

%% Compute PML Coefficients
xVec = (1:pmlWidth)-.5;
xVecShift = xVec+.5;

% Basic Vectors
sigVec = (abs(xVec+.5).^sigM./pmlWidth.^sigM).*sigMax;
sigVecStag = (abs(xVec).^sigM./pmlWidth.^sigM).*sigMax;

kappaVec = 1+(abs(xVec+.5).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVecStag = 1+(abs(xVec).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);

aVec = (abs(xVec+.5).^aM./pmlWidth.^aM).*aMax;
aVecStag = (abs(xVec).^aM./pmlWidth.^aM).*aMax;

% Full vector set spanning grid space
sigVecD = [sigVec(end:-1:1) zeros(1,gridSize-2*length(sigVec)) sigVec];
sigVecStagD = [sigVecStag(end:-1:1) zeros(1,gridSize-2*length(sigVec)-1) sigVecStag];

kappaVecD = [kappaVec(end:-1:1) ones(1,gridSize-2*length(kappaVec)) kappaVec];
kappaVecStagD = [kappaVecStag(end:-1:1) ones(1,gridSize-2*length(kappaVec)-1) kappaVecStag];

aVecD = [aVec(end:-1:1) zeros(1,gridSize-2*length(aVec)) aVec];
aVecStagD = [aVecStag(end:-1:1) zeros(1,gridSize-2*length(aVec)-1) aVecStag];

tau = (kappaVecD.*epsNot)./(kappaVecD.*aVecD+sigVecD);
tauStag = (kappaVecStagD.*epsNot)./(kappaVecStagD.*aVecStagD+sigVecStagD);

% Compute coefficients for use in update
bx = exp(-dt./tau);
by = exp(-dt./tauStag);

gx = (sigVecD./(kappaVecD.*(kappaVecD.*aVecD+sigVecD))).*(1-exp(-dt./tau));
gy = (sigVecStagD./(kappaVecStagD.*(kappaVecStagD.*aVecStagD+sigVecStagD))).*(1-exp(-dt./tauStag));

kx = kappaVecD;
ky = kappaVecStagD;

for k = 1:length(gx)
    if isnan(gx(k))
        gx(k) = 0;
    end
end

for k = 1:length(gy)
    if isnan(gy(k))
        gy(k) = 0;
    end
end

%% Initialization and settings for microscopic medium

%%%%%%%%%%%%%%%%%%%%%% Solution vector initilization %%%%%%%%%%%%%%%%%%%%%% 
% Auxiliary u vectors -- roughly equivalent to :
%                        rho_n = exp(-t/T_m)*u_n

% Each matrix, u, contains the same number of rows as there are response
% frequencies, and the same number of columns as the gridsize.
u1 = zeros(numFreq,gridSize);
u2 = zeros(numFreq,gridSize);
u3 = zeros(numFreq,gridSize);

% Coefficient vectors
A = zeros(numFreq,1);
B = zeros(numFreq,1);
Cp = 0;
Cm = 0;
D = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omge = 2*pi*fe;         
medPos = round(locM/dx);

time = [];
record = [];
intenseR = [];
intenseL = [];
recordRho = [];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
figure('outerposition', [70   70   1200   1000])
set(1,'color',[.9 .9 .9])

for n = 1:maxIter
    
    % Source Condition  
    if pulseSetting == 1            % Gaussian Source
        source = amp*exp(-(maxT-n*dt)^2/(width/5)^2);
    elseif pulseSetting == 2        % Modulated Gaussian Source
        source = amp*exp(-(maxT-n*dt)^2/(width/5)^2) * amp*sin(omega*n*dt);
    elseif pulseSetting == 3        % CW source
        source = amp*sin(omega*n*dt);
    elseif pulseSetting == 4        % Ricker Wavelet source
        source = amp*(1-2*(pi*fp*(n*dt-dr))^2)*exp(-(pi*fp*(n*dt-dr))^2);
    elseif pulseSetting == 5
        source = amp*sech((tpS*f/2)*(dt*n-tpS/2)/(tpS/10)) *sin(omega*n*dt);
        
    elseif pulseSetting == 6
        ft = 0;
        smoothFactor = (1-(dt*n/(tp*startUp)-1).^2).^4;
        if dt*n < tp*startUp
            gt = sin(omega*(n*dt-td))*smoothFactor;
        elseif dt*n >= tp*startUp
            gt = sin(omega*(n*dt-td));
        end
        
        source = amp*(ft + gt);

    end
    intenseR(n) = max(abs(Ex(medPos(2)+1:end)));
    intenseL(n) = max(abs(Ex(1:medPos(1))));
    time(n) = dt*n;
    
    if isnan(Ex(round(gridSize/2)))
        return
    end
    
    % Field Updates with microscopic medium
    tHalf = dt*(n+1/2);
    t     = dt*(n+1);
    display(n)
    % Compute new coefficient values for time dependent coefficients
    % These are at n+1/2
    A  = (  (Natom.*gammaN)/(epsNot*T2)  ) * exp( -(tHalf)/T2 );
    B  = ( (Natom.*gammaN.*omge)/(epsNot) ) * exp( -(tHalf)/T2 );
    Cm = ( (2*gammaN)/dirac ) * exp( -(tHalf)*(1/T2-1/T1) );
    
    Cp = ( (2*gammaN)/dirac ) * exp( -(t)*(1/T1-1/T2) );
    D  = ( (2*gammaN*rho30)/dirac )* exp( t/T2 );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% H-field Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    QEx = by.*QEx - gy.*(Ex(2:end) - Ex(1:end-1))/dx;
    Hy = Hy - dy1.*ky.*((Ex(2:end) - Ex(1:end-1))/dx + QEx);

    u1Store = u1;
    u2Store = u2;
    u3Store = u3;
    ExStore  = Ex;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% E-field Update %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Sum the bloch vectors to create the polarization interaction
    P1 = A'*u1;
    P2 = B'*u2;
    
    QHy(2:end-1) = bx(2:end-1).*QHy(2:end-1) - gx(2:end-1).*(Hy(2:end) - Hy(1:end-1))/dx;
    Ex(2:end-1)  = cx1(2:end-1).*Ex(2:end-1) - cx2(2:end-1).*kx(2:end-1).*( ((Hy(2:end) - Hy(1:end-1))/dx) + QHy(2:end-1) ) ...
                   + dt.* ( - P1(2:end-1) +  P2(2:end-1) );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% BLOCH Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make temporary Ex vectors to be used for multiplication with bloch
    % matrix.
    
    ExN = repmat(ExStore,numFreq,1);
    ExN1 = repmat(Ex,numFreq,1);
    
    u3 = u3 - ( dt*Cm/2 ).*( (ExN1 + ExN).*(u2) );
    
    % Implicit solution to u1 which allows u2 to be updated
    Q = ( Cp.*u3.*ExN1 + D.*ExN1 );
    
    % Create a matrix of omegas to use during computation
    omgeM = repmat(omge,1,gridSize);
    u1 = ( ( 1 - (omgeM.^2.*dt^2/4) ).*u1 + omgeM.*dt.*u2 + ( (omgeM.*dt^2)/2 ).*Q )./( 1 + (omgeM.^2*dt.^2)/4 );
    u2 = u2 - ( (omgeM.*dt)/2 ).*( u1 + u1Store ) + dt.*Q;

    % Remove values if not part of the medium
    u1(:,1:medPos(1)) = zeros(numFreq,medPos(1));
    u2(:,1:medPos(1)) = zeros(numFreq,medPos(1));
    u3(:,1:medPos(1)) = zeros(numFreq,medPos(1));
    u1(:,medPos(2):end) = zeros(numFreq,gridSize - medPos(2)+1);
    u2(:,medPos(2):end) = zeros(numFreq,gridSize - medPos(2)+1);
    u3(:,medPos(2):end) = zeros(numFreq,gridSize - medPos(2)+1);
    rho3 = (rho30 + exp(-dt*n/T1)*u3)/normalizer3;
    
    % Source Update
    Ex(sourcePos) = source;


    % E-field Boundary Update
    Ex(1)   = 0 ;
    Ex(end) = 0 ;
    
    % Frequency Domain via Discrete Fourier Transform
    if n < stopTime
        EfInit = EfInit + Ex(detInitial).*exp(-1i.*2.*pi.*dt*n.*testFreq);
    else
        EfRef = EfRef + Ex(detRef).*exp(-1i.*2.*pi.*dt*n.*testFreq);
    end
        EfFinal = EfFinal + Ex(detPost).*exp(-1i.*2.*pi.*dt*n.*testFreq);

    % Continuous plotting
    
    if mod(n,skipFrames) == 0
        %%%%%%%%%%%%%%%%%%%%%% Generate plots/movie %%%%%%%%%%%%%%%%%%%%%%
        clf
        disConv = dx*10^6;
        recColor = [.7 .7 .7];
        lineColor1 = [.8 0 0];
        
        subplot(2,2,1)
        hold on
        % E-field plot
        rectangle('position', [medPos(1)*disConv  -1.7 (medPos(2)*disConv-medPos(1)*disConv) 4], 'facecolor', recColor , 'edgecolor', 'none')
                    line([0 gridSize*disConv], [1 1] )
        plot(xDom*disConv, Ex/amp, 'color', lineColor1, 'linewidth', 1);
        title('1D FDTD Simulation', 'fontsize', 18, 'fontname', 'cordia new')
        xlabel('[{\mu}m]', 'fontsize', 18, 'fontname', 'cordia new')
        ylabel('E_x', 'fontsize', 18, 'fontname', 'cordia new')
        axis([ 0 gridSize*disConv -1.5 1.5])

        % Plot detectors
        plot(detInitial*disConv, 0, 'ro', detPost*disConv, 0, 'go' , detRef*disConv, 0, 'b+')
        plot(detInitial*disConv, 0, 'ro', detPost*disConv, 0, 'go' , detRef*disConv, 0, 'b+')
        % Display current time
        timeElaps = num2str(round(n*dt*10^16)*10^-1);
        timeElapsed = strcat('Time elapsed:{ } ', timeElaps, 'fs');
        text(gridSize/1.3*disConv, .8, timeElapsed); 
        gainVal = num2str(max(Ex)/amp-1);
        gainVal2 = strcat('Max(Ex)-1:{ } ', gainVal);
        text(gridSize/1.3*disConv, .5, gainVal2); 
        
        hold off

            subplot(2,2,3)
            
            hold on
            rectangle('position', [medPos(1)*disConv  -1.2 (medPos(2)*disConv-medPos(1)*disConv) 4], 'facecolor', recColor , 'edgecolor', 'none')

            rho3 =(rho30 + exp(-t/T1)*u3)/normalizer3;
            rho3(:,1:medPos(1)) = zeros(numFreq,medPos(1));
            rho3(:,medPos(2):end) = zeros(numFreq,gridSize - medPos(2) + 1);
            plot(xDom*disConv, rho3, 'color', lineColor1, 'linewidth', 1)
            axis([ 0 gridSize*disConv -1.1 1.1])
            title('{\rho}_3', 'fontsize', 18, 'fontname', 'cordia new')
            hold off

            subplot(2,2,2)
            axis( [ 0 1e12*time(n) .99 1.06])
            hold on
            rectangle('position', [medPos(1)*disConv  -1.2 (medPos(2)*disConv-medPos(1)*disConv) 4], 'facecolor', recColor , 'edgecolor', 'none')
            rho1 = exp(-t/T2)*u1/normalizer;
            rho2 = exp(-t/T2)*u2/normalizer;
            plot(xDom*disConv, rho1, 'color', lineColor1, 'linewidth', 1)
            plot(xDom*disConv, rho2, 'b-', 'linewidth', 1)
            title('{\rho}_1 / {\rho}_2', 'fontsize', 18, 'fontname', 'cordia new')
            axis([ 0 gridSize*disConv -1.1 1.1])
            hold off

            subplot(2,2,4) 
            axis( [ 0 1e12*time(n) 0 1])
            hold on
            rectangle('position', [medPos(1)*disConv  -1.2 (medPos(2)*disConv-medPos(1)*disConv) 4], 'facecolor', recColor , 'edgecolor', 'none')
            plot (xDom*disConv, (rho1.^2 +rho2.^2 + rho3.^2), 'g', 'linewidth', 1)
            title('{\rho}_1^2 + {\rho}_2^2 + {\rho}_3^2', 'fontsize', 18, 'fontname', 'cordia new')
            axis([ 0 gridSize*disConv -1.1 1.1])
            hold off

            drawnow
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Generation and plotting of frequency domain results 

% Amplitude in frequency domain
magI = abs(EfInit).^2;
magF = abs(EfFinal).^2;
magR = abs(EfRef).^2;

% Plot freq. spectrum
figure('outerposition', [50   50   700   400])
hold on
plot(testFreq*1e-12, magI, 'b')
plot(testFreq*1e-12, magF, 'g')
plot(testFreq*1e-12, magR, 'r')
plot(testFreq*1e-12, magR+magF, 'k:')
ylabel('Power') 
xlabel('{\lambda} [nm]')
legend('Inital Spectrum', 'Transmitted Spectrum', 'Reflected spectrum', 'Sum of transmiited and reflected spectra')
title('Spectrum')
axis([min(testFreq)*1e-12, max(testFreq)*1e-12,  0 1.4*max([max(magI), max(magR), max(magI)])]) 

% Plot Transmission
trans = magF./magI;
refl = magR./magI;
figure('outerposition', [50   450   700   400])
plot(testFreq*1e-12, refl, testFreq*1e-12, trans)
ylabel('Transmitted/Reflected') 
xlabel('{\lambda} [nm]')
legend('Reflectance (E^r/E^i)^2', 'Transmittance (E^t/E^i)^2')
title('Transmission/Reflection')
axis([min(testFreq)*1e-12, max(testFreq)*1e-12,  0 1.2*max([max(trans), max(refl)])]) 


figure('outerposition', [50   850   700   400])
plot(testFreq*1e-12, 1-trans)
ylabel('Abosrbance [normalized]') 
xlabel('Frequency [THz]')
title('Medium absorbance')
