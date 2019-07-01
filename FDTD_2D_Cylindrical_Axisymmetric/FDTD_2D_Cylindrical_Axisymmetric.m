clear all;
close all;

run colorScheme.m;

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
% Simulates wave propagation in Axisymmetric Cylindrical coordinates      %
% using the FDTD method in two dimensions (r,z)                           %
%                                                                         %
% 1. Utilizes all six field components (Er, Ez, Ephi and Hr, Hz, Hphi)    %
% 2. Symmetry is governed by azimuthal harmonic index "m"                 %
%    Commonly m = 0 or m = 1 (See Mittra articles)                        %
% 3. CPML utilized to truncate spacial domain                             %
% 4. Source conditions can be varied                                      %
% 5. Coefficient vectors allow simple implementation of refractive        %
%    indices and/or simple lossy media                                    %
% 6. PEC and other waveguides/optical devices have been implemented,      %
%    but are not included here to prevent cluttering the code.            %
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
rangeR      = 4e-6;             % r Space size (1 = 1meter)                [USER DEFINED]
rangeZ      = 10e-6;            % z Space size (1 = 1meter)                [USER DEFINED]
timeDom     = 70e-15;           % Time domain  (1 = 1 sec)                 [USER DEFINED]
m = 1;                          % Azimuthal Harmonic Index                 [USER DEFINED]                               
gridR       = 100;              % Number of cells along r                  [USER DEFINED]       
gridZ       = 300;              % Number of cells along z                  [USER DEFINED]

dr          = rangeR/gridR;       
dz          = rangeZ/gridZ;     
dt          = .3*min([dr, dz])/(cNot);   
maxIter     = round(timeDom/dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
skipframes  = 10;               % Number of iters to skip while graphing   [USER DEFINED]
startFrame  = 1;                % Number of iters to skip before graphing  [USER DEFINED]

% "surfaceView" value determines which fields to display
  surfaceView = 1;
% 0 = All Fields
% 1 = Ephidouble
% 2 = Erdouble
% 3 = Ezdouble
% 4 = Beam profile plot
% 5 = Total Efield (i.e. E_T = Ephi + Ez + Er)

%%%%%%%%%%%%%%%%%%%%%%%% Hard source parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian = 1,        Modulated Gaussian = 2,      Sine Wave = 3, 
% SECH pulse = 4;      Gaussian Beam = 5;

pulseSetting = 5;               % Choose which source condition to use     [USER DEFINED]
sourcePosR   = [1e-6 2e-6];     % r Location of source in meters           [USER DEFINED]
sourcePosZ   = [1e-6 1e-6];     % z Location of source in meters           [USER DEFINED]        
f            = 200e12;          % Frequency                                [USER DEFINED]
amp          = 1;               % Amplitude of Wave                        [USER DEFINED]        
sourceR      = round(sourcePosR/dr);
sourceZ      = round(sourcePosR/dz);

%%%%%%% Gaussian Pulse parameters %%%%%%%
width        = 50e-15;         % Approximate width of pulse in time (s)   [USER DEFINED]
maxT         = width/2; 
%%%%%%%%% Sine Wave Parameters %%%%%%%%%%
lambda       = cNot/f;          % Wavelength
omega        = 2*pi*f;          % Angular Frequency
%%%%%%%%% SECH pulse Parameters %%%%%%%%%
tp           = 20/f;            % Approximate width of pulse in time (s)   [USER DEFINED]
%%%%%%% Gaussian Beam parameters %%%%%%%%
widthBeamFocus2 = 1e-6;          % Width of Gaussian                       [USER DEFINED]
z = -2e-6;                       % Starting position along beam axis       [USER DEFINED]
pulsedBeam = 1;                  % 1 = pulsed beam, 0 = CW beam            [USER DEFINED]

%%%%%%%%%%%%%%%%%%%%%%% Frequency domain parameters %%%%%%%%%%%%%%%%%%%%%%%

%%% INCOMPLETE %%%
lambdaWindow = [ .1, 10]*1e-6;
freqMin = cNot/lambdaWindow(1);
freqMax = cNot/lambdaWindow(2);
steps = 100;
detZ = [sourceZ(1)+100, gridZ - 40];
detZ = [40, 120];
%%% INCOMPLETE %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PML Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmlWidth     = 20;              % Width of PML in grid points              [USER DEFINED]
sigM         = 3;               % Power of sigma in PML                    [USER DEFINED]
sigMax       = 5*(sigM+1)*...   % Max value of sigma in PML
              .8/(etaNot*dr);
kappaMax     = 1;               % Max value of kappa in PML                [USER DEFINED]
kappaM       = 3;               % Power of kappa in PML                    [USER DEFINED]
aMax         = 1;               % Max value of kappa in PML                [USER DEFINED]
aM           = 3;               % Power of a in PML                        [USER DEFINED]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Other Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
planePos = 8e-6;               % Position of cross sectional viewing plane
planePos = round(planePos/dz);
zMax1 = 1;
zMax2 = 1;
cMin = -0.5;
cMax = 0.5;
cVar = [cMin cMax];
axiconPos = 2e-6;
displayC = 0;    % Exits before FDTD loop and Displays PML coefficients, 
                 % 0 = no display, 1 = display
                

                 
%% Calculations for Gaussian Beam and Plots
widthBeamFocus = widthBeamFocus2./dr;
totalW = gridR-2;
zR = pi.*widthBeamFocus2.^2/lambda;
rVec = linspace(0,rangeR,gridR);
R = z.*(1+(zR./z).^2);
zVec = -sqrt(R.^2-(rVec).^2)-R;
rScaled = round(rVec./dr);
zScaled = round(zVec./dr);

positionsB = zeros(gridR,2);
for x = 1:gridR-1
       positionsB(x,1) = rScaled(x+1);
       positionsB(x,2) = zScaled(x+1)+sourceZ(1);
       if imag(positionsB(x,2)) ~= 0
           positionsB(x,2) = positionsB(x-1,2);
       end
end

% Make Phi space for plot
gridPhi = gridR-1;
phiEx = linspace(0, 2*pi, gridPhi);
phi = phiEx(1:end);
X = [];
Y = [];

for k = 1:gridPhi;
    
    [T S] = pol2cart(phi(k).*ones(1,gridR-1), linspace(0, rangeR, gridR-1));
    X = [X ; T];
    Y = [Y ; S];
    
end



%% Initialize solution vectors
% General Cylindrical coordinate form: f(r,phi,z)
% Reduced here to: f(r,z)
Hr   = zeros(gridR  , gridZ-1 );
Hz   = zeros(gridR  , gridZ   );
Hphi = zeros(gridR  , gridZ-1 );

Er   = zeros(gridR  , gridZ  );
Ez   = zeros(gridR  , gridZ-1);
Ephi = zeros(gridR  , gridZ  );

Jr   = zeros(gridR  , gridZ  );
Jz   = zeros(gridR  , gridZ-1);
Jphi = zeros(gridR  , gridZ  );

sigr   = zeros(gridR  , gridZ  );
sigz   = zeros(gridR  , gridZ-1);
sigphi = zeros(gridR  , gridZ  );

epsr   = ones(gridR  , gridZ  );
epsz   = ones(gridR  , gridZ-1);
epsphi = ones(gridR  , gridZ  );

mur    = ones(gridR  , gridZ-1);
muz    = ones(gridR  , gridZ  );
muphi  = ones(gridR  , gridZ-1);

QEzphi = zeros(size(Hr));
QEzr   = zeros(size(Hphi));
QErz   = zeros(size(Hphi));
QErphi = zeros(size(Hz));
QHzphi = zeros(size(Er));
QHrphi = zeros(size(Ez)); 
QHzr   = zeros(size(Ephi));
QHrz   = zeros(size(Ephi));

% Frequency Domain Initilization
testFreq = linspace(freqMin, freqMax, steps);
testWL = cNot./testFreq;
detR = 1:gridR-pmlWidth;
detL = length(detR);
ErInit  = zeros(detL, length(testFreq));
ErFinal = zeros(detL, length(testFreq));
EphiInit  = zeros(detL, length(testFreq));
EphiFinal = zeros(detL, length(testFreq));
EzInit  = zeros(detL, length(testFreq));
EzFinal = zeros(detL, length(testFreq));



%% Compute coefficient vectors
% E-Field coeff.
cr1 = (1-(sigr.*dt)./(2*epsNot.*epsr))./(1+(sigr.*dt)./(2*epsNot.*epsr));
cr2 = (dt./(epsNot.*epsr))./(1+(sigr.*dt)./(2*epsNot.*epsr));
cr3 = ((m*dt)./(epsNot.*epsr))./(1+(sigr.*dt)./(2*epsNot.*epsr));
cr4 = ((4*m*dt)./(epsNot.*epsr(1,:)))./(1+(sigr(1,:).*dt)./(2*epsNot.*epsr(1,:)));

cz1 = (1-(sigz.*dt)./(2*epsNot.*epsz))./(1+(sigz.*dt)./(2*epsNot.*epsz));
cz2 = ((m*dt)./(epsNot.*epsz))./(1+(sigz.*dt)./(2*epsNot.*epsz));
cz3 = ((dt)./(epsNot.*epsz))./(1+(sigz.*dt)./(2*epsNot.*epsz));

cphi1 = (1-(sigphi.*dt)./(2*epsNot.*epsphi))./(1+(sigphi.*dt)./(2*epsNot.*epsphi));
cphi2 = (dt./(epsNot.*epsphi))./(1+(sigphi.*dt)./(2*epsNot.*epsphi));
cphi3 = (dt./(epsNot.*epsphi))./(1+(sigphi.*dt)./(2*epsNot.*epsphi));

jr   = dt./(epsr.*epsNot);
jphi = dt./(epsphi.*epsNot);
jz   = dt./(epsz.*epsNot);

% D-Field coeff.
dr1 = (m*dt)./(muNot.*mur);
dr2 = dt./(muNot*mur);

dz1 = (m*dt)./(muNot.*muz(2:end,:));
dz2 = dt./(muNot*muz(2:end,:));
dz3 = (4*dt)./(muNot.*muz(1,:));

dphi1 = dt./(muNot*muphi(2:end,:));
dphi2 = dt./(muNot*muphi(1,:));

rVec = dr.*(1:gridR);
riHalf = repmat(rVec(1:end-1)', 1, gridZ );
ri     = repmat(rVec(1:end)', 1, gridZ ) - .5*dr;



%% Compute PML Coefficients
xVec = (1:pmlWidth)-.5;
gridPlotZ = (1:gridZ)-1;
gridPlotStagZ = (1:gridZ-1)-.5;
gridPlotR = (1:gridR)-1;
gridPlotStagR = (1:gridR-1)-.5;

sigVec = (abs(xVec+.5).^sigM./pmlWidth.^sigM).*sigMax;
sigVecZ = [sigVec(end:-1:1) zeros(1,gridZ-2*length(sigVec)) sigVec];
sigVecR = [zeros(1,gridR-length(sigVec)) sigVec];
sigVecStag = (abs(xVec).^sigM./pmlWidth.^sigM).*sigMax;
sigVecStagZ = [sigVecStag(end:-1:1) zeros(1,gridZ-2*length(sigVec)-1) sigVecStag];
sigVecStagR = [zeros(1,gridR-length(sigVec)-1) sigVecStag];

kappaVec = 1+(abs(xVec+.5).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVecZ = [kappaVec(end:-1:1) ones(1,gridZ-2*length(kappaVec)) kappaVec];
kappaVecR = [ones(1,gridR-length(kappaVec)) kappaVec];
kappaVecStag = 1+(abs(xVec).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVecStagZ = [kappaVecStag(end:-1:1) ones(1,gridZ-2*length(kappaVec)-1) kappaVecStag];
kappaVecStagR = [ones(1,gridR-length(kappaVecStag)-1) kappaVecStag];

aVec = (abs(xVec+.5).^aM./pmlWidth.^aM).*aMax;
aVecZ = [aVec(end:-1:1) zeros(1,gridZ-2*length(aVec)) aVec];
aVecR = [zeros(1,gridR-length(aVec)) aVec];
aVecStag = (abs(xVec).^aM./pmlWidth.^aM).*aMax;
aVecStagZ = [aVecStag(end:-1:1) zeros(1,gridZ-2*length(aVec)-1) aVecStag];
aVecStagR = [zeros(1,gridR-length(aVecStag)-1) aVecStag];

% Plots of vectors used in PML
if displayC == 1;
figure('outerposition', [100, 400, 1500, 600])
subplot(1,3,1)
plot(gridPlotZ, sigVecZ, 'ro', gridPlotStagZ, sigVecStagZ, 'b+' , gridPlotR, sigVecR, 'r.', gridPlotStagR, sigVecStagR, 'b*')
title('PML {\sigma}')

subplot(1,3,2)
plot(gridPlotZ, kappaVecZ, 'ro', gridPlotStagZ, kappaVecStagZ, 'b+', gridPlotR, kappaVecR, 'r.', gridPlotStagR, kappaVecStagR, 'b*')
title('PML {\kappa}')

subplot(1,3,3)
plot(gridPlotZ, aVecZ, 'ro', gridPlotStagZ, aVecStagZ, 'b+', gridPlotR, aVecR, 'r.', gridPlotStagR, aVecStagR, 'b*')
title('PML a')
end

% Generate coefficient vectors used for PML during computation
tau = (kappaVecZ.*epsNot)./(kappaVecZ.*aVecZ+sigVecZ);
tauStag = (kappaVecStagZ.*epsNot)./(kappaVecStagZ.*aVecStagZ+sigVecStagZ);
bVec = exp(-dt./tau);
bVecStag = exp(-dt./tauStag);
gVec = (sigVecZ./(kappaVecZ.*(kappaVecZ.*aVecZ+sigVecZ))).*(1-exp(-dt./tau));
gVecStag = (sigVecStagZ./(kappaVecStagZ.*(kappaVecStagZ.*aVecStagZ+sigVecStagZ))).*(1-exp(-dt./tauStag));

for k = 1:length(gVec)
    if isnan(gVec(k))
        gVec(k) = 0;
    end
end

for k = 1:length(gVecStag)
    if isnan(gVecStag(k))
        gVecStag(k) = 0;
    end
end

bz1 = repmat(bVec,gridR,1);
bz2 = repmat(bVecStag,gridR,1);
gz1 = repmat(gVec,gridR,1);
gz2 = repmat(gVecStag,gridR,1);
kz1 = repmat(1./kappaVecZ,gridR,1);
kz2 = repmat(1./kappaVecStagZ,gridR,1);

tauR = (kappaVecR.*epsNot)./(kappaVecR.*aVecR+sigVecR);
tauStagR = (kappaVecStagR.*epsNot)./(kappaVecStagR.*aVecStagR+sigVecStagR);
bVecR = exp(-dt./tauR);
bVecStagR = exp(-dt./tauStagR);

gVecR = (sigVecR./(kappaVecR.*(kappaVecR.*aVecR+sigVecR))).*(1-exp(-dt./tauR));
gVecStagR = (sigVecStagR./(kappaVecStagR.*(kappaVecStagR.*aVecStagR+sigVecStagR))).*(1-exp(-dt./tauStagR));

for k = 1:length(gVecR)
    if isnan(gVecR(k))
        gVecR(k) = 0;
    end
end

for k = 1:length(gVecStagR)
    if isnan(gVecStagR(k))
        gVecStagR(k) = 0;
    end
end

br1 = repmat(bVecR,gridZ,1)';
br2 = repmat(bVecStagR,gridZ,1)';
gr1 = repmat(gVecR,gridZ,1)';
gr2 = repmat(gVecStagR,gridZ,1)';
kr1 = repmat(1./kappaVecR,gridZ,1)';
kr2 = repmat(1./kappaVecStagR,gridZ,1)';

if displayC == 1
figure('outerposition', [1 100 1900 1000])
subplot(2,4,1)
imagesc(bz1')
xlabel('r')
ylabel('z')
title('bz1')
colorbar
axis square

subplot(2,4,2)
imagesc(bz2')
xlabel('r')
ylabel('z')
title('bz2')
colorbar
axis square

subplot(2,4,3)
imagesc(gz1')
xlabel('r')
ylabel('z')
title('gz1')
colorbar
axis square

subplot(2,4,4)
imagesc(gz2')
xlabel('r')
ylabel('z')
title('gz2')
colorbar
axis square

subplot(2,4,5)
imagesc(br1')
xlabel('r')
ylabel('z')
title('br1')
colorbar
axis square

subplot(2,4,6)
imagesc(br2')
xlabel('r')
ylabel('z')
title('br2')
colorbar
axis square

subplot(2,4,7)
imagesc(gr1')
xlabel('r')
ylabel('z')
title('gr1')
colorbar
axis square

subplot(2,4,8)
imagesc(gr2')
xlabel('r')
ylabel('z')
title('gr2')
colorbar
axis square
return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
 figure('outerposition', [100 50 500 800])   
set(1,'color',[.9 .9 .9])

for n = 1:maxIter
    
        
% Source Condition  
    if pulseSetting == 1            % Gaussian Source
        source = amp*exp(-(maxT-n*dt)^2/(width/5)^2);
        Ephi(sourceR(1):sourceR(2), sourceZ(1):sourceZ(2)) = source   ;
        
    elseif pulseSetting == 2        % Modulated Gaussian Source
        source = amp*exp(-(maxT-n*dt)^2/(width/5)^2) * amp*sin(omega*n*dt);
        Ephi(sourceR(1):sourceR(2), sourceZ(1):sourceZ(2)) = source   ;
        
    elseif pulseSetting == 3        % CW source
        source = amp*sin(omega*n*dt);
        Ephi(sourceR(1):sourceR(2), sourceZ(1):sourceZ(2)) = source   ;
        
    elseif pulseSetting == 4
        source = amp*sech((tp*f/2)*(dt*n-tp/2)/(tp/2)) *sin(omega*n*dt);
        Ephi(sourceR(1):sourceR(2), sourceZ(1):sourceZ(2)) = source   ;
        
    elseif pulseSetting == 5
        Bvec = 0:round(totalW);
        widthBeam = widthBeamFocus.*(sqrt(1+(z./zR).^2));
        if pulsedBeam == 0
            gauss = 1;
        elseif pulsedBeam == 1
            gauss = amp*exp(-(maxT-n*dt)^2/(width/5)^2);
        end
        source = amp.*widthBeamFocus/widthBeam.*exp(- Bvec.^2 /widthBeam.^2)*sin(omega*n*dt).*gauss;
        for k = 1:length(source)
            Hphi(positionsB(k,1),positionsB(k,2)) = source(k).*2.105e-3;
            Hr(positionsB(k,1),positionsB(k,2))   = - source(k).*2.105e-3;
            Hz(positionsB(k,1),positionsB(k,2))   = - source(k).*2.105e-3;
        end
    end
    
    
% Field Update

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  H-field Updates  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%% Hr FIELD %%%%%%%%%%
    QEzphi = bz2.*QEzphi - gz2.*(Ephi(:,2:end) - Ephi(:,1:end-1))./dz ;  
    Hr = Hr - ...
         dr1.*Ez./ri(:,1:end-1) +...
         dr2.*(kz2.*(Ephi(:,2:end) - Ephi(:,1:end-1))./dz + QEzphi);
     
    %%%%%%%%%% Hz FIELD %%%%%%%%%%
    QErphi(2:end,:) = br2.*QErphi(2:end,:) - gr2.*(ri(2:end,:).*Ephi(2:end,:) - ri(1:end-1,:).*Ephi(1:end-1,:))./dr;
    Hz(2:end,:) = Hz(2:end,:) + ...
         dz1.*Er(2:end,:)./riHalf - ...
         dz2.*(kr2.*(ri(2:end,:).*Ephi(2:end,:) - ri(1:end-1,:).*Ephi(1:end-1,:))./dr + QErphi(2:end,:))./riHalf;
     
    % r=0 Singularity Update for m = 0
    if m == 0
    Hz(1,:) = Hz(1,:) - dz3.*Ephi(1,:)./dr; 
    end
    
    % r=0 Singularity Update for m >= 1
    if m > 0
    Hz(1,:) = zeros(size(Hz(1,:)));    
    end
    
    %%%%%%%%%% Hphi FIELD %%%%%%%%%%
    QEzr(2:end,:) = bz2(2:end,:).*QEzr(2:end,:) - gz2(2:end,:).*(Er(2:end,2:end) - Er(2:end,1:end-1))./dz;
    QErz(2:end,:) = br2(:,2:end).*QErz(2:end,:) - gr2(:,2:end).*(Ez(2:end,:) - Ez(1:end-1,:))./dr;
    Hphi(2:end,:) = Hphi(2:end,:) - ...
                    dphi1.*( kz2(2:end,:).*(Er(2:end,2:end) - Er(2:end,1:end-1))./dz - kr2(:,2:end).*(Ez(2:end,:) - Ez(1:end-1,:))./dr + QEzr(2:end,:)  - QErz(2:end,:) );    
    % r=0 Singularity Update   
    Hphi(1,:) = zeros(size(Hphi(1,:)));             
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  E-field Updates  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%% Er FIELD %%%%%%%%%%
    QHzphi(2:end,2:end-1) = bz1(2:end, 2:end-1).*QHzphi(2:end, 2:end-1) - gz1(2:end, 2:end-1).*(Hphi(2:end,2:end) - Hphi(2:end,1:end-1))./dz;
    
    Er(2:end,2:end-1) = cr1(2:end,2:end-1).*Er(2:end,2:end-1) - ...
                        cr2(2:end,2:end-1).*(kz1(2:end, 2:end-1).*(Hphi(2:end,2:end) - Hphi(2:end,1:end-1))./dz + QHzphi(2:end,2:end-1) )- ...
                        cr3(2:end,2:end-1).*Hz(2:end,2:end-1)./riHalf(:,2:end-1) + ...
                        jr(2:end,2:end-1).*Jr(2:end,2:end-1);

    %%%%%%%%%% Ez FIELD %%%%%%%%%%
    QHrphi(1:end-1,:) = br1(1:end-1,2:end).*QHrphi(1:end-1,:) - gr1(1:end-1,2:end).*(riHalf(1:end,2:end).*Hphi(2:end,:) - [zeros(1, gridZ-1); riHalf(1:end-1,2:end)].*Hphi(1:end-1,:))./dr;
    Ez(1:end-1,:) = cz1(1:end-1,:).*Ez(1:end-1,:) + ...
                    cz2(1:end-1,:).*Hr(1:end-1,:)./ri(1:end-1,2:end) + ...
                    cz3(1:end-1,:).*(kr1(1:end-1,2:end).*(riHalf(1:end,2:end).*Hphi(2:end,:) - [zeros(1, gridZ-1); riHalf(1:end-1,2:end)].*Hphi(1:end-1,:))./(dr) + QHrphi(1:end-1,:) )./ri(1:end-1,2:end) + ...
                    jz(1:end-1,:).*Jz(1:end-1,:);

    %%%%%%%%%% Ephi FIELD %%%%%%%%%%
    QHzr(1:end-1,2:end-1) = bz1(1:end-1,2:end-1).*QHzr(1:end-1,2:end-1) - gz1(1:end-1,2:end-1).*(Hr(1:end-1,2:end) - Hr(1:end-1,1:end-1))./dz; 
    QHrz(1:end-1,2:end-1) = br1(1:end-1,2:end-1).*QHrz(1:end-1,2:end-1) - gr1(1:end-1,2:end-1).*(Hz(2:end,2:end-1) - Hz(1:end-1,2:end-1))./dr;
    Ephi(1:end-1,2:end-1) = cphi1(1:end-1,2:end-1).*Ephi(1:end-1,2:end-1) + ...
                            cphi2(1:end-1,2:end-1).*(kz1(1:end-1,2:end-1).*(Hr(1:end-1,2:end) - Hr(1:end-1,1:end-1))./dz + QHzr(1:end-1,2:end-1)) - ...
                            cphi3(1:end-1,2:end-1).*(kr1(1:end-1,2:end-1).*(Hz(2:end,2:end-1) - Hz(1:end-1,2:end-1))./dr + QHrz(1:end-1,2:end-1)) + ...
                            jphi(1:end-1,2:end-1).*Jphi(1:end-1,2:end-1);            

    %%%%%%%%%% Boundary Update %%%%%%%%%%
    Er(:,1)     = 0;
    Er(:,end)   = 0;
    Ez(end,:)   = 0;
    Ephi(:,1)   = 0;
    Ephi(:,end) = 0;
    Ephi(end,:) = 0;
    

%% Frequency Domain Update

    ErInit = ErInit + Er(detR, detZ(1))*exp(-1i.*2.*pi.*dt*n.*testFreq);
    ErFinal = ErFinal + Er(detR, detZ(2))*exp(-1i.*2.*pi.*dt*n.*testFreq);

    EphiInit = EphiInit + Ephi(detR, detZ(1))*exp(-1i.*2.*pi.*dt*n.*testFreq);
    EphiFinal = EphiFinal + Ephi(detR, detZ(2))*exp(-1i.*2.*pi.*dt*n.*testFreq);

    EzInit = EzInit + Ez(detR, detZ(1))*exp(-1i.*2.*pi.*dt*n.*testFreq);
    EzFinal = EzFinal + Ez(detR, detZ(2))*exp(-1i.*2.*pi.*dt*n.*testFreq);


%% Continuous plotting

    if mod(n, skipframes) == 0 && n>startFrame
    display(strcat('n = ', num2str(n)))

    if surfaceView == 0
        subplot(2,3,1)
        imagesc(flipud(real(Er)'))
        title('Er')
        caxis([cMin cMax])
        colormap(kMap)

        subplot(2,3,2)
        imagesc(flipud(real(Ez)'))
        title('Ez')
        caxis([cMin cMax])
        colormap(kMap)

        subplot(2,3,3)
        imagesc(flipud(real(Ephi)'))
        title('Ephi')
        caxis([cMin cMax])
        colormap(kMap)

        subplot(2,3,4)
        imagesc(flipud(real(Hr)'))
        title('Hr')
        caxis([cMin cMax].*3e-3) 
        colormap(kMap)

        subplot(2,3,5)
        imagesc(flipud(real(Hz)'))
        title('Hz')
        caxis([cMin cMax].*3e-3)
        colormap(kMap)

        subplot(2,3,6)
        imagesc(flipud(real(Hphi)'))
        title('Hphi')
        caxis([cMin cMax].*3e-3) 
        colormap(kMap)
        hold on

    elseif surfaceView == 1

        imagesc(flipud([fliplr(Ephi(2:end,:)') , Ephi(2:end,:)' ]))
        xlabel('r', 'fontsize', 16)
        ylabel('z', 'fontsize', 16)
        title('E_{\phi}', 'fontsize', 20)
        caxis([cMin cMax])
        colormap(kMap)
        axis equal  

    elseif surfaceView == 2

        imagesc(flipud([fliplr(Er(2:end,:)') , Er(2:end,:)' ]))
        xlabel('r', 'fontsize', 16)
        ylabel('z', 'fontsize', 16)
        title('E_r', 'fontsize', 20)
        caxis([cMin cMax])
        colormap(kMap)
        axis equal  
        
    elseif surfaceView == 3
        
        imagesc(flipud([fliplr(Ez') , Ez' ]))
        xlabel('r', 'fontsize', 16)
        ylabel('z', 'fontsize', 16)
        title('E_z', 'fontsize', 20)
        caxis([cMin cMax])
        colormap(kMap)
        axis equal     

    elseif surfaceView == 4
        Enew = Ephi(2:end,2:end)+Er(2:end,2:end)+Ez(2:end,:) ;      
        subplot(1,3,1)
        imagesc(flipud([fliplr((Enew(2:end,:))') , (Enew(2:end,:))' ]))
        colormap(kMap)
        xlabel('r', 'fontsize', 16)
        ylabel('z', 'fontsize', 16)
        title('E_r', 'fontsize', 20)
        caxis([cMin cMax])
        axis equal      
        
        subplot(1,3,2)
        surf(X,Y, repmat( Enew(:,planePos), 1 , length(Y) )', 'edgealpha', .05)
        axis([min(min(X)) max(max(X)) min(min(Y)) max(max(Y)) -.5 .5])
        caxis([-.5 .5])
        colormap(kMap)

        subplot(1,3,3)
        I = abs([Enew(gridR-10:-1:2, planePos); Enew(2:gridR-10, planePos) ]).^2;
        plot(I)
        axis([0 198 -2 2])

    elseif surfaceView == 5
        Enew = Er(:,2:end) + Ephi(:,2:end) + Ez ;
        imagesc(flipud([fliplr((Enew(2:end,:))') , (Enew(2:end,:))' ]))
        colormap(kMap)
        xlabel('r', 'fontsize', 16)
        ylabel('z', 'fontsize', 16)
        title('E_T', 'fontsize', 20)
        caxis([cMin cMax])
        axis equal      

    end
    
        s = getframe;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generation and plotting of frequency domain results 

magIr = abs(ErInit).^2;
magFr = abs(ErFinal).^2;
rescale = max(max(magIr));
magIr = magIr/rescale;
magFr = magFr/rescale;

figure('outerposition', [100, 700, 900, 500])
subplot(1,2,1)
surf( testFreq/1e12, detR*dr*1e2, magIr, 'edgealpha', .1)
title('Initial')
xlabel('{\omega} [THz]')
ylabel('r [cm]')

subplot(1,2,2)
surf( testFreq/1e12, detR*dr*1e2, magFr, 'edgealpha', .1)
title('Final')
xlabel('{\omega} [THz]')
ylabel('r [cm]')




    