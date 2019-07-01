clear all;
close all;

% Physical Constants
epsNot = 8.854187817620e-12;    % Permittivity of Free Space
muNot = 1.25663706e-6;          % Magnetic Constant
cNot = 2.99792458e8;            % Speed of light
etaNot = sqrt(muNot/epsNot);    % Characteristic impedence of free space

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Info:                                                           %
% Simulates wave propagation in 1D using the FDTD method                  %
% 1. Utilizes Ex and Hy fields propagating in z direction                 %
% 2. CPML utilized to truncate spacial domain                             %
% 3. Source condition can be varied in the settings                       %
% 4. Coefficient vectors allow simple implementation of refractive        %
%    indices and/or simple lossy media                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters and Settings
%%%%%%%%%%%%%%%%%%%%%%% Define FDTD solution space %%%%%%%%%%%%%%%%%%%%%%%% 
spaceSize    = 14.0e-6;         % Spactial Domain  (1 = 1 meter)           [USER DEFINED]
timeDom      = 70e-15;           % Time domain      (1 = 1 sec)             [USER DEFINED]
gridSize     = 1000;            % Total number of spatial cells            [USER DEFINED]
dx           = spaceSize/gridSize;                 
dt           = .5*dx/(cNot);  
maxIter      = round(timeDom/dt);    

%%%%%%%%%%%%%%%%%%%%%%%% Hard source parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian = 1,  Modulated Gaussian = 2,  Sine Wave = 3, 
% Ricker Wavelet = 4;
pulseSetting = 2;               % Choose which source condition to use     [USER DEFINED]
sourcePos    = 2e-6 ;           % Location of source in meters             [USER DEFINED]
sourcePos    = round(sourcePos/dx);
% Gaussian Pulse parameters
width        = 10e-15;          % Approximate width of pulse in time (sec) [USER DEFINED]
maxT         = width*2;
% Sine Wave Parameters
f            = 460e12;          % Frequency                                [USER DEFINED]
amp          = 1;               % Amplitude of Wave                        [USER DEFINED]
lambda       = cNot/f;          % Wavelength
omega        = 2*pi*f;          % Angular Frequency
% Ricker Wavelet Parameters
fp           = 5e8;             % Peak Frequency                           [USER DEFINED]
md           = 1;               % Temporal Delay Multiple                  [USER DEFINED]
dr           = md/fp;           % Temporal Delay

%%%%%%%%%%%%%%%%%%%%%%% Frequency domain parameters %%%%%%%%%%%%%%%%%%%%%%%
lambdaWindow = [ 500, 800]*1e-9;
stepsF       = 1000;
detInitial   = sourcePos+1;     % Position of initial detector             [USER DEFINED]
detPost      = gridSize-120;    % Position of latter detector              [USER DEFINED]
detRef       = detInitial+10;   % Position of reflection detector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PML Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmlWidth     = 10;              % Width of PML in grid points              [USER DEFINED]
sigM         = 3;               % Power of sigma in PML                    [USER DEFINED]
sigMax       = 5*(sigM+1)*...   % Max value of sigma in PML
              .8/(etaNot*dx);
kappaMax     = 1;               % Max value of kappa in PML                [USER DEFINED]
kappaM       = 3;               % Power of kappa in PML                    [USER DEFINED]
aMax         = 1;               % Max value of kappa in PML                [USER DEFINED]
aM           = 3;               % Power of a in PML                        [USER DEFINED]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Other Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stopTime = 2000;
skipFrames = 50;
posMed = 3e-6;                  % Position of medium relative to source    [USER DEFINED]                              
posMed = round(posMed/dx)+sourcePos;
visible = [0 14.0 -1 1];
visible2 = [0 14 -.2 .2];
stoppingNorm = .0002;
switchAxis = 125000;

% Print current settings
display(strcat('Spacial domain =  ',  num2str(round(spaceSize*1e7)*1e-1), ' microns'))
display(strcat('dx =  ',  num2str(round(dx*1e10)*1e-1), ' nm'))
display(strcat('Time domain =  ',  num2str(round(maxIter*dt*1e16)*1e-1), ' fs'))
display(strcat('Time step =  ',  num2str(round(dt*1e19)*1e-4), ' fs'))

%% Frequency domain parameters
freqMin = cNot/lambdaWindow(1);
freqMax = cNot/lambdaWindow(2);
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




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
figure('outerposition', [110   100   1100   600])
set(1,'color',[.9 .9 .9])

for n = 1:maxIter
    %% Source Condition  
    if pulseSetting == 1            % Gaussian Source
        source = exp(-(maxT-n*dt)^2/(width/5)^2);
    elseif pulseSetting == 2        % Modulated Gaussian Source
        source = exp(-(maxT-n*dt)^2/(width/5)^2) * amp*sin(omega*n*dt);
    elseif pulseSetting == 3        % CW source
        source = amp*sin(omega*n*dt);
    elseif pulseSetting == 4        % Ricker Wavelet source
        source = (1-2*(pi*fp*(n*dt-dr))^2)*exp(-(pi*fp*(n*dt-dr))^2);
    end

    %% Field Updates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% H-field Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    QEx = by.*QEx - gy.*(Ex(2:end) - Ex(1:end-1))/dx;
    Hy = Hy + dy1.*ky.*((Ex(2:end) - Ex(1:end-1))/dx + QEx);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% E-field Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    QHy(2:end-1) = bx(2:end-1).*QHy(2:end-1) - gx(2:end-1).*(Hy(2:end) - Hy(1:end-1))/dx;
    Ex(2:end-1)  = cx1(2:end-1).*Ex(2:end-1) + cx2(2:end-1).*kx(2:end-1).*( ((Hy(2:end) - Hy(1:end-1))/dx) + QHy(2:end-1) );

    % Source Update
    if n*dt < maxT*2
    Ex(sourcePos) = source;
    end

    % E-field Boundary Update
    Ex(1)   = 0 ;
    Ex(end) = 0 ;

    %% Frequency Domain via Discrete Fourier Transform

    if n < stopTime
        EfInit = EfInit + Ex(detInitial).*exp(-1i.*2.*pi.*dt*n.*testFreq);
    else
        EfRef = EfRef + Ex(detRef).*exp(-1i.*2.*pi.*dt*n.*testFreq);
    end
        EfFinal = EfFinal + Ex(detPost).*exp(-1i.*2.*pi.*dt*n.*testFreq);

    %% Continuous plotting
    
    if mod(n,skipFrames) == 0 
        %%%%%%%%%%%%%%%%%%%%%% Generate plots/movie %%%%%%%%%%%%%%%%%%%%%%
        clf
        hold on
        disConv = dx*10^6;
        
        % E-field plot
        plot(xDom*disConv, Ex, 'color', [ 0 0 .6], 'linewidth', 2) 
        title('1D FDTD Simulation', 'fontsize', 18, 'fontname', 'cordia new')
        xlabel('[{\mu}m]', 'fontsize', 18, 'fontname', 'cordia new')
        ylabel('E_x', 'fontsize', 18, 'fontname', 'cordia new')
        axis([ 0 gridSize*disConv -1 1])

        % Plot detectors
        plot(detInitial*disConv, 0, 'ro', detPost*disConv, 0, 'go' , detRef*disConv, 0, 'b+')
        plot(detInitial*disConv, 0, 'ro', detPost*disConv, 0, 'go' , detRef*disConv, 0, 'b+')
       

        % Display current time
        timeElaps = num2str(round(n*dt*10^16)*10^-1);
        timeElapsed = strcat('Time elapsed:{ } ', timeElaps, 'fs');
        text(gridSize/1.2*disConv, .105, timeElapsed); 

        % Makes run slower but shows graph each time
        M = getframe(gcf);
        % Movie recording goes here

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Generation and plotting of frequency domain results 

% Amplitude in frequency domain
magI = abs(EfInit).^2;
magF = abs(EfFinal).^2;
magR = abs(EfRef).^2;

figure('outerposition', [10   600   600   600])
% Plot freq. spectrum
hold on
plot(testWL*1e9, magI, 'y')
plot(testWL*1e9, magF, 'g')
plot(testWL*1e9, magR, 'r')
plot(testWL*1e9, magR+magF, 'k:')
hold off
ylabel('Power') 
xlabel('{\lambda} [nm]')
legend('Inital Spectrum', 'Transmitted Spectrum', 'Reflected spectrum', 'Sum of transmiited and reflected spectra')
title('Spectrum')
axis([1e9*min(testWL), 1e9*max(testWL),  0 1.4*max([max(magI), max(magR), max(magI)])]) 

% Plot Transmission
trans = magF./magI;
refl = magR./magI;
figure('outerposition', [700   600   600   600])
plot(testWL*1e9, refl, testWL*1e9, trans)
ylabel('Transmitted/Reflected') 
xlabel('{\lambda} [nm]')
legend('Reflectance (E^r/E^i)^2', 'Transmittance (E^t/E^i)^2')
title('Transmission/Reflection')
axis([1e9*min(testWL), 1e9*max(testWL),  0 1.2*max([max(trans), max(refl)])]) 
