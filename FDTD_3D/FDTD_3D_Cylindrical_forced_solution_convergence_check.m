clear all;
close all;

run colorSchemes.m

% Physical Constants
epsNot = 8.854187817620e-12;        % Permittivity of Free Space
muNot  = 1.25663706e-6;             % Magnetic Constant
cNot   = 2.99792458e8;              % Speed of light in vacuum
etaNot = sqrt(muNot/epsNot);        % Characteristic impedence of free space
planck = 6.62606957e-34;            % Plancks constant
dirac  = planck/(2*pi);             % Dirac constant (Reduced Plancks constant)
e = 1.60217646e-19;                 % Charge of an electron (Coulombs)
me = 9.10938188e-31;                % Electron Mass (kg)

epsNot = 1;
muNot = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Info:                                                           %
% Tests the cylindrical coordinate system finite diff. solver against by  %
% forcing an analytic solution and comparing the result to the true sltn  %
%                                                                         %
%                                                                         %
%  This Source Code Form is subject to the terms of the Mozilla Public    %
%  License, v. 2.0. If a copy of the MPL was not distributed with this    %
%  file, You can obtain one at http://mozilla.org/MPL/2.0/.               %
%                                                                         %
%  Copyright (C) 2012-2013 Brian D. Hong                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


start = 5;
steps = 6;
points = [];
totalPoints = [];
error = [];
maxs = [];
counter = 0;

for MULT = start:steps+start
    counter = counter +1;
%% Parameters and Settings
%%%%%%%%%%%%%%%%%%%%%%% Define FDTD solution space %%%%%%%%%%%%%%%%%%%%%%%% 
rangeR      = 8e-6;             % r Space size (1 = 1meter)                [USER DEFINED]
rangeZ      = 12e-6;             % z Space size (1 = 1meter)                [USER DEFINED]
timeDom     = 650e-15;          % Time domain  (1 = 1 sec)                 [USER DEFINED]

rangeR = 5;
rangeZ = 10;
timeDom = 3;
L = 10;

gridR       = 8*MULT;               % Number of cells along r                  [USER DEFINED]       
gridZ       = 8*MULT;              % Number of cells along z                  [USER DEFINED]
gridPHI     = 4*2*MULT;              % Number of cells about axis (must be 4*k) [USER DEFINED]

points(counter) = gridR;
totalPoints(counter) = gridR*gridZ*gridPHI;

dr          = rangeR/gridR;       
dz          = rangeZ/gridZ;    
dphi        = 2*pi/gridPHI;
% dt          = .5*min([dr, dz])/(cNot);   

dt          = .5/( 1/(sqrt(epsNot*muNot)) .* sqrt((1/dr.^2) + (1/(dr/2*dphi).^2) + (dz)^2 ) );  

maxIter     = round(timeDom/dt);
maxIter = 30;

% dt = .05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
skipframes  = 1;               % Number of iters to skip while graphing   [USER DEFINED]
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

pulseSetting = 1;               % Choose which source condition to use     [USER DEFINED]
sourcePosR   = [3e-6 3e-6];     % r Location of source in meters           [USER DEFINED]
sourcePosZ   = [3e-6 3e-6];     % z Location of source in meters           [USER DEFINED]    
sourcePosPhi = [pi   pi];     % PHI Location of source in radians        [USER DEFINED] 
f            = 200e12;          % Frequency                                [USER DEFINED]
amp          = 1;               % Amplitude of Wave                        [USER DEFINED]        
sourceR      = round(sourcePosR/dr);
sourceZ      = round(sourcePosR/dz);
sourcePhi    = round(sourcePosPhi/dphi);

sourceR = [20 40]
sourceZ = [20 20]
sourcePhi = [1 gridPHI]
%%%%%%% Gaussian Pulse parameters %%%%%%%
width        = 5e-15;         % Approximate width of pulse in time (s)   [USER DEFINED]
maxT         = width; 
%%%%%%%%% Sine Wave Parameters %%%%%%%%%%
lambda       = cNot/f;          % Wavelength
omega        = 2*pi*f;          % Angular Frequency
%%%%%%%%% SECH pulse Parameters %%%%%%%%%
tp           = 20/f;            % Approximate width of pulse in time (s)   [USER DEFINED]
%%%%%%% Gaussian Beam parameters %%%%%%%%
widthBeamFocus2 = 1e-6;          % Width of Gaussian                       [USER DEFINED]
z = -5e-6;                       % Starting position along beam axis       [USER DEFINED]
pulsedBeam = 0;                  % 1 = pulsed beam, 0 = CW beam            [USER DEFINED]

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
pmlWidth     = 0;              % Width of PML in grid points              [USER DEFINED]
sigM         = 3;               % Power of sigma in PML                    [USER DEFINED]
sigMax       = (sigM+1)*...   % Max value of sigma in PML
              .8/(etaNot*dr);
kappaMax     = 1;               % Max value of kappa in PML                [USER DEFINED]
kappaM       = 3;               % Power of kappa in PML                    [USER DEFINED]
aMax         = 1;               % Max value of kappa in PML                [USER DEFINED]
aM           = 3;               % Power of a in PML                        [USER DEFINED]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Other Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
planeNumber = 5;               % Number of viewing planes

planePos = 8e-6;               % Position of cross sectional viewing plane
planePos = round(planePos/dz);
zMax1 = 1;
zMax2 = 1;
cMin = -.5;
cMax = .5;
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

% positionsB = zeros(gridR-2,2);
% for x = 1:gridR-2
%        positionsB(x,1) = rScaled(x+1);
%        positionsB(x,2) = zScaled(x+1)+sourceZ(1);
%        if imag(positionsB(x,2)) ~= 0
%            positionsB(x,2) = positionsB(x-1,2);
%        end
% end

phiEx = linspace(0, 2*pi, gridPHI);
phi = phiEx(1:end);
X = [];
Y = [];
for k = 1:gridPHI;

    
    [T S] = pol2cart(phi(k).*ones(1,gridR), 0:gridR-1);
    X = [X ; T];
    Y = [Y ; S];
    
end



%% Initialize solution vectors
% General Cylindrical coordinate form: f(r,phi,z)
Hr     = zeros(gridR   , gridPHI, gridZ-1 );
Hz     = zeros(gridR-1 , gridPHI, gridZ   );
Hphi   = zeros(gridR-1 , gridPHI, gridZ-1 );

Er     = zeros(gridR-1 , gridPHI, gridZ  );
Ez     = zeros(gridR   , gridPHI, gridZ-1);
Ephi   = zeros(gridR   , gridPHI, gridZ  );

Jr     = zeros(gridR-1 , gridPHI, gridZ  );
Jz     = zeros(gridR   , gridPHI, gridZ-1);
Jphi   = zeros(gridR   , gridPHI, gridZ  );

sigr   = zeros(gridR-1 , gridPHI, gridZ  );
sigz   = zeros(gridR   , gridPHI, gridZ-1);
sigphi = zeros(gridR   , gridPHI, gridZ  );

epsr   = ones(gridR-1 , gridPHI, gridZ  );
epsz   = ones(gridR   , gridPHI, gridZ-1);
epsphi = ones(gridR   , gridPHI, gridZ  );

mur    = ones(gridR   , gridPHI, gridZ-1);
muz    = ones(gridR-1 , gridPHI, gridZ  );
muphi  = ones(gridR-1 , gridPHI, gridZ-1);

%% Create optical or other physical elements


% Generate a waveguide with higher refractive index then the cladding.
% nIndex = 3;
% epsilons = nIndex.^2;
% widthI = 20;
% 
% epsr(1:widthI, :,:) = epsilons;
% epsz(1:widthI, :,:) = epsilons;
% epsphi(1:widthI, :,:) = epsilons;

%% Compute coefficient vectors
% E-Field coeff.
cr1 = (1-(sigr.*dt)./(2*epsNot.*epsr))./(1+(sigr.*dt)./(2*epsNot.*epsr));
cr2 = (dt./(epsNot.*epsr))./(1+(sigr.*dt)./(2*epsNot.*epsr));

cz1 = (1-(sigz.*dt)./(2*epsNot.*epsz))./(1+(sigz.*dt)./(2*epsNot.*epsz));
cz2 = ((dt)./(epsNot.*epsz))./(1+(sigz.*dt)./(2*epsNot.*epsz));

cphi1 = (1-(sigphi.*dt)./(2*epsNot.*epsphi))./(1+(sigphi.*dt)./(2*epsNot.*epsphi));
cphi2 = (dt./(epsNot.*epsphi))./(1+(sigphi.*dt)./(2*epsNot.*epsphi));

jr   = dt./(epsr.*epsNot);
jphi = dt./(epsphi.*epsNot);
jz   = dt./(epsz.*epsNot);

% D-Field coeff.
dr1 = dt./(muNot.*mur);
dz1 = dt./(muNot.*muz);
dphi1 = dt./(muNot.*muphi);

iNot     = 0:gridR-1;
iNotHalf = iNot(2:end)-1/2;

iN       = zeros(gridR, gridPHI, gridZ);
iNH      = zeros(gridR-1, gridPHI, gridZ);

QEzphi = zeros(size(Hr));
QEzr   = zeros(size(Hphi));
QErz   = zeros(size(Hphi));
QErphi = zeros(size(Hz));
QHzphi = zeros(size(Er));
QHrphi = zeros(size(Ez)); 
QHzr   = zeros(size(Ephi));
QHrz   = zeros(size(Ephi));

for phi = 1:gridPHI
    for k = 1:gridZ
        iN(:,phi,k) = iNot;
        iNH(:,phi,k) = iNotHalf;
    end
end

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

bz1T = repmat(bVec,gridR,1);
bz2T = repmat(bVecStag,gridR,1);
gz1T = repmat(gVec,gridR,1);
gz2T = repmat(gVecStag,gridR,1);
kz1T = repmat(1./kappaVecZ,gridR,1);
kz2T = repmat(1./kappaVecStagZ,gridR,1);

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

br1T = repmat(bVecR,gridZ,1)';
br2T = repmat(bVecStagR,gridZ,1)';
gr1T = repmat(gVecR,gridZ,1)';
gr2T = repmat(gVecStagR,gridZ,1)';
kr1T = repmat(1./kappaVecR,gridZ,1)';
kr2T = repmat(1./kappaVecStagR,gridZ,1)';


bz1 = zeros(gridR  , gridPHI, gridZ  ) ; 
bz2 = zeros(gridR  , gridPHI, gridZ-1) ; 
br1 = zeros(gridR  , gridPHI, gridZ  ) ;
br2 = zeros(gridR-1, gridPHI, gridZ  ) ; 
gz1 = zeros(gridR  , gridPHI, gridZ  ) ; 
gz2 = zeros(gridR  , gridPHI, gridZ-1) ; 
gr1 = zeros(gridR  , gridPHI, gridZ  ) ; 
gr2 = zeros(gridR-1, gridPHI, gridZ  ) ;
kz1 = zeros(gridR  , gridPHI, gridZ  ) ; 
kz2 = zeros(gridR  , gridPHI, gridZ-1) ;
kr1 = zeros(gridR  , gridPHI, gridZ  ) ;
kr2 = zeros(gridR-1, gridPHI, gridZ  ) ;

for k = 1:gridPHI;
    bz1(:,k,:) = bz1T;
    bz2(:,k,:) = bz2T;
    br1(:,k,:) = br1T;
    br2(:,k,:) = br2T;
    gz1(:,k,:) = gz1T;
    gz2(:,k,:) = gz2T;
    gr1(:,k,:) = gr1T;
    gr2(:,k,:) = gr2T;
    kz1(:,k,:) = kz1T;
    kz2(:,k,:) = kz2T;
    kr1(:,k,:) = kr1T;
    kr2(:,k,:) = kr2T;
end

if displayC == 1;
figure('outerposition', [1 100 1900 1000])
subplot(2,4,1)
imagesc(bz1T')
xlabel('r')
ylabel('z')
title('bz1')
colorbar
axis square

subplot(2,4,2)
imagesc(bz2T')
xlabel('r')
ylabel('z')
title('bz2')
colorbar
axis square

subplot(2,4,3)
imagesc(gz1T')
xlabel('r')
ylabel('z')
title('gz1')
colorbar
axis square

subplot(2,4,4)
imagesc(gz2T')
xlabel('r')
ylabel('z')
title('gz2')
colorbar
axis square

subplot(2,4,5)
imagesc(br1T')
xlabel('r')
ylabel('z')
title('br1')
colorbar
axis square

subplot(2,4,6)
imagesc(br2T')
xlabel('r')
ylabel('z')
title('br2')
colorbar
axis square

subplot(2,4,7)
imagesc(gr1T')
xlabel('r')
ylabel('z')
title('gr1')
colorbar
axis square

subplot(2,4,8)
imagesc(gr2T')
xlabel('r')
ylabel('z')
title('gr2')
colorbar
axis square
return
end

%% Exact solution
% ErFun = @(r,phi,z,t,L) r.^2.*exp(-r.^2).*exp(-(z-t-L/2).^2).*cos(phi);
% EphiFun = @(r,phi,z,t,L) r.^6.*exp(-r.^4).*exp(-(z-t-L/2).^2).*sin(phi.^2);
% EzFun = @(r,phi,z,t,L) r.^3.*exp(-r.^4).*exp(-(z-t-L/2).^2).*sin((phi.^2-phi.^3)/10);

ErFun = @(r,phi,z,t,L) r.^2.*exp(-r.^2).*exp(-(z-t-L/2).^2).*cos(phi);
EphiFun = @(r,phi,z,t,L) r.^6.*exp(-r.^4).*exp(-(z-t-L/2).^2).*sin(phi);
EzFun = @(r,phi,z,t,L) (-1/2).*exp( (-r.^2).*(1+r.^2) ).*sqrt(pi).*r.*(exp(r.^2).*r.^4 + exp(r.^4).*(3-2.*r.^2)).*cos(phi).*( erf(L/2 + t) - erf(L/2 + t - z) );


HrFun = @(r,phi,z,t,L, mu) ( 1/(mu) ).* (...
        exp(-r.^4).*( exp(-(1/4).*(L-2*z).^2) - exp(-(1/4).*(L+2*t-2*z).^2) ).* (r.^6) .*  sin(phi) - (1/4).*exp( -r.^2.*(1+r.^2) ).*sqrt(pi).*...
        ( exp(r.^2).*r.^4 + exp(r.^4).*(3-2.*r.^2) ).*( (2/sqrt(pi)).*(-exp(-L.^2/4) + exp(-(1/4).*(L+2.*t).^2) + exp(-(1/4).*(L-2*z).^2) - exp(-(1/4).*(L+2*t-2*z).^2 ) ) -...
        L.*erf(L/2) + (L + 2.*t).* erf(L/2 + t) + (L - 2.*z).*erf(L/2-z) - (L+2.*t -2.*z).*erf(L/2 + t - z) ).*sin(phi) );                                                 
HphiFun = @(r,phi,z,t,L, mu) -(1/(4*mu)).*exp( -r.^2.*(1+r.^2) ).*cos(phi).*(...
          exp(-(1/4).*(L+2*t).^2 - z.^2) .* (2.*exp(L.*z).*( exp(t.*(L+t)) - exp(2.*t.*z) ) .* ( exp(r.^4).*(3 - 10*r.^2 + 4*r.^4) + exp(r.^2).*( 5.*r.^4 - 4.*r.^8 )) - ...
          2.*exp(z.^2).*(-1+exp(t.*(L+t))).*( exp(r.^2) .* ( 5.*r.^4 - 4.*r.^8 ) + exp(r.^4) .* (3 + 4*r.^2.*(-3 + r.^2))  )) + ...
          sqrt(pi).*( exp(r.^2).*(5.*r.^4-4.*r.^8) + exp(r.^4).*(3 + 4*r.^2.*(-3 + r.^2))).* ( L.*erfc(L/2) - (L+2.*t).*erfc(L/2+t) - ...
          (L-2*z) .* erfc(L/2-z) + (L+2*t-2*z).*erfc(L/2 + t - z) ) );
      
HzFun = @(r,phi,z,t,L, mu) (1/(2*mu)).*(...
        exp( -r.^2.*(1+r.^2) ).*sqrt(pi).*r.*( exp(r.^4) + exp(r.^2) .* ( 7.*r.^4 - 4.*r.^8 ) ).*( erf(L/2 - z) - erf(L/2 + t - z) ).*sin(phi) );

JrFun = @(r,phi,z,t,L,mu,eps)  ( 1/(mu) ).*exp(-r.^2.*(1+r.^2)).*cos(phi).* ( ...
        exp(r.^4 - (1/4).*(L+2.*t).^2 + L.*z - z.^2 ) .* r.^2 .* (exp(t.*(L+t)) .* (L-2.*z) + exp(2.*t.*z).*(-1+eps*mu).*(L+2*t-2*z) ) + ...
        sqrt(pi).*(-exp(r.^2).*r.^4 + exp(r.^4).*(1-6*r.^2 + 2.*r.^4)).* ( - erf(L/2 - z) + erf(L/2 + t - z) ) );
JphiFun = @(r,phi,z,t,L,mu,eps)  ( 1/(mu) ).*exp(-r.^2.*(1+r.^2) - (1/4).*(L+2.*t-2.*z).^2 ).* (...
          exp(r.^2).*r.^6.*( exp(t.*(L+t-2.*z)) .* ( L - 2.*z) + (-1 + eps*mu).*(L+2*t-2*z) ) + exp((L/2+t-z).^2) .*sqrt(pi) .* ...
          ( exp(r.^4) - exp(r.^2).*r.^4.*(17-32.*r.^4+8.*r.^8) ).* ( erf(L/2 - z) - erf(L/2 + t - z) ) ).*sin(phi);  
JzFun = @(r,phi,z,t,L,mu,eps)  ( 1/(mu) ).*exp(-r.^2.*(1+r.^2) ).*r.*cos(phi).*(...
         exp( -(1/4).*(L+2.*t).^2 - z.^2 ) .* ( exp(z.^2) .* ( exp(r.^4 + t.*(L+t)) .* (-20 + 22.*r.^2 - 4.*r.^4) + exp(r.^4).*(20 + 3.*eps.*mu - 2.*(11 + eps*mu).*r.^2 + 4.*r.^4) + ...
         exp(r.^2).*r.^2.*(-12 + eps*mu*r.^2 + 28*r.^4 - 8*r.^8) + 4.*exp(r.^2 + t.*(L+t)).*r.^2.*(3-7*r.^4 + 2.*r.^8) ) + ...
         exp(L*z).*( exp(2.*t.*z).* ( exp(r.^4).* (-17-3.*eps.*mu + 2.*(10+eps*mu).*r.^2 - 4.*r.^4 )  + ...
         exp(r.^2).*r.^2.*(12 + (1-eps*mu).*r.^2 - 28.*r.^4 + 8.*r.^8) ) + ...
         exp(t.*(L+t)).*(-exp(r.^2).*r.^2.*(12 + r.^2 - 28.*r.^4 + 8.*r.^8 ) + exp(r.^4).*(17+4*(r.^2).*(-5+r.^2))))) + ...
         sqrt(pi).*( exp(r.^4).*(10-11.*r.^2+2.*r.^4) - 2.*exp(r.^2).*(r.^2).*(3-7.*r.^4+2.*r.^8)).*( L.*erfc(L/2) - (L+2.*t).*erfc(L/2+t) - (L-2*z) .* erfc(L/2-z) + (L+2*t-2*z).*erfc(L/2 + t - z) ) );

 
%% Source Condition  
n = 0;
L2 = L;
% Time positions
t1  = (n)*dt;
t12 = (n-1/2)*dt;

% Radial positions (matrix)
r1  = zeros(gridR, gridPHI, gridZ);
r12 = zeros(gridR-1, gridPHI, gridZ);
for k = 1:gridR
    r1(k,:,:) = dr.*(k-1);
end
for k = 1:gridR-1;

    r12(k,:,:) = dr*(k-1/2);
    
end

% z-positions (matrix)
z1 = zeros(gridR, gridPHI, gridZ); 
z12 = zeros(gridR, gridPHI, gridZ-1); 
for k = 1:gridZ
    z1(:,:,k) = dz.*(k-1);
end
for k = 1:gridZ-1;
    z12(:,:,k) = dz.*(k-1/2);
end

% phi-positions (matrix)
phi1 = zeros(gridR, gridPHI, gridZ); 
phi12 = zeros(gridR, gridPHI, gridZ); 
for k = 1:gridPHI
    phi1(:,k,:) = dphi.*(k);
end
for k = 1:gridPHI;
    phi12(:,k,:) = dphi.*(k+1/2);
end


HrExact = HrFun(r1(:,:,2:end),phi12(:,:,2:end),z12,t12,L,muNot) ;
HzExact = HzFun(r12,phi12(2:end,:,:),z1(2:end,:,:),t12,L,muNot) ;
HphiExact = HphiFun(r12(:,:,2:end),phi1(2:end,:,2:end),z12(2:end,:,:),t12,L,muNot) ;

EzExact = EzFun(r1(:,:,2:end),phi1(:,:,2:end),z12,t1,L) ;
ErExact = ErFun(r12,phi1(2:end,:,:),z1(2:end,:,:),t1,L) ;
EphiExact = EphiFun(r1,phi12,z1,t1,L) ;


JzExact = JzFun(r1(:,:,2:end),phi1(:,:,2:end),z12,t12,L,muNot,epsNot) ;
JrExact = JrFun(r12,phi1(2:end,:,:),z1(2:end,:,:),t12,L,muNot,epsNot) ;
JphiExact = JphiFun(r1,phi12,z1,t12,L,muNot,epsNot) ;

Jz = JzExact;
Jr = JrExact;
Jphi = JphiExact;

Hr = HrExact;
Hz = HzExact;
Hphi = HphiExact;
Ez = EzExact;
Er = ErExact;
Ephi = EphiExact;


     %% SYMBOLIC
     
     syms r phi z t L mu eps

ErS = r.^2.*exp(-r.^2).*exp(-(z-t-L/2).^2).*cos(phi); %% CHECKS OUT


EphiS = r.^6.*exp(-r.^4).*exp(-(z-t-L/2).^2).*sin(phi); %% CHECKS OUT
EzS = (-1/2).*exp( (-r.^2).*(1+r.^2) ).*sqrt(pi).*r.*(exp(r.^2).*r.^4 + exp(r.^4).*(3-2.*r.^2)).*cos(phi).*( erf(L/2 + t) - erf(L/2 + t - z) ); %% CHECKS OUT


HrS = ( 1/(mu) ).* (...
        exp(-r.^4).*( exp(-(1/4).*(L-2*z).^2) - exp(-(1/4).*(L+2*t-2*z).^2) ).* (r.^6) .*  sin(phi) - (1/4).*exp( -r.^2.*(1+r.^2) ).*sqrt(pi).*...
        ( exp(r.^2).*r.^4 + exp(r.^4).*(3-2.*r.^2) ).*( (2/sqrt(pi)).*(-exp(-L.^2/4) + exp(-(1/4).*(L+2.*t).^2) + exp(-(1/4).*(L-2*z).^2) - exp(-(1/4).*(L+2*t-2*z).^2 ) ) -...
        L.*erf(L/2) + (L + 2.*t).* erf(L/2 + t) + (L - 2.*z).*erf(L/2-z) - (L+2.*t -2.*z).*erf(L/2 + t - z) ).*sin(phi) );  %% CHECKS OUT  
    
    
HphiS =  -(1/(4*mu)).*exp( -r.^2.*(1+r.^2) ).*cos(phi).*(...
          exp(-(1/4).*(L+2*t).^2 - z.^2) .* (2.*exp(L.*z).*( exp(t.*(L+t)) - exp(2.*t.*z) ) .* ( exp(r.^4).*(3 - 10*r.^2 + 4*r.^4) + exp(r.^2).*( 5.*r.^4 - 4.*r.^8 )) - ...
          2.*exp(z.^2).*(-1+exp(t.*(L+t))).*( exp(r.^2) .* ( 5.*r.^4 - 4.*r.^8 ) + exp(r.^4) .* (3 + 4*r.^2.*(-3 + r.^2))  )) + ...
          sqrt(pi).*( exp(r.^2).*(5.*r.^4-4.*r.^8) + exp(r.^4).*(3 + 4*r.^2.*(-3 + r.^2))).* ( L.*erfc(L/2) - (L+2.*t).*erfc(L/2+t) - ...
          (L-2*z) .* erfc(L/2-z) + (L+2*t-2*z).*erfc(L/2 + t - z) ) );  %% CHECKS OUT

HzS = (1/(2*mu)).*(...
        exp( -r.^2.*(1+r.^2) ).*sqrt(pi).*r.*( exp(r.^4) + exp(r.^2) .* ( 7.*r.^4 - 4.*r.^8 ) ).*( erf(L/2 - z) - erf(L/2 + t - z) ).*sin(phi) ); %% CHECKS OUT

JrS = ( 1/(mu) ).*exp(-r.^2.*(1+r.^2)).*cos(phi).* ( ...
        exp(r.^4 - (1/4).*(L+2.*t).^2 + L.*z - z.^2 ) .* r.^2 .* (exp(t.*(L+t)) .* (L-2.*z) + exp(2.*t.*z).*(-1+eps*mu).*(L+2*t-2*z) ) + ...
        sqrt(pi).*(-exp(r.^2).*r.^4 + exp(r.^4).*(1-6*r.^2 + 2.*r.^4)).* ( - erf(L/2 - z) + erf(L/2 + t - z) ) );%% CHECKS OUT
    
JphiS = ( 1/(mu) ).*exp(-r.^2.*(1+r.^2) - (1/4).*(L+2.*t-2.*z).^2 ).* (...
          exp(r.^2).*r.^6.*( exp(t.*(L+t-2.*z)) .* ( L - 2.*z) + (-1 + eps*mu).*(L+2*t-2*z) ) + exp((L/2+t-z).^2) .*sqrt(pi) .* ...
          ( exp(r.^4) - exp(r.^2).*r.^4.*(17-32.*r.^4+8.*r.^8) ).* ( erf(L/2 - z) - erf(L/2 + t - z) ) ).*sin(phi); %% CHECKS OUT  

%CHECKS OUT     
JzS =  ( 1/(mu) ).*exp(-r.^2.*(1+r.^2) ).*r.*cos(phi).*(...
         exp( -(1/4).*(L+2.*t).^2 - z.^2 ) .* ( exp(z.^2) .* ( exp(r.^4 + t.*(L+t)) .* (-20 + 22.*r.^2 - 4.*r.^4) + exp(r.^4).*(20 + 3.*eps.*mu - 2.*(11 + eps*mu).*r.^2 + 4.*r.^4) + ...
         exp(r.^2).*r.^2.*(-12 + eps*mu*r.^2 + 28*r.^4 - 8*r.^8) + 4.*exp(r.^2 + t.*(L+t)).*r.^2.*(3-7*r.^4 + 2.*r.^8) ) + ...
         exp(L*z).*( exp(2.*t.*z).* ( exp(r.^4).* (-17-3.*eps.*mu + 2.*(10+eps*mu).*r.^2 - 4.*r.^4 )  + ...
         exp(r.^2).*r.^2.*(12 + (1-eps*mu).*r.^2 - 28.*r.^4 + 8.*r.^8) ) + ...
         exp(t.*(L+t)).*(-exp(r.^2).*r.^2.*(12 + r.^2 - 28.*r.^4 + 8.*r.^8 ) + exp(r.^4).*(17+4*(r.^2).*(-5+r.^2))))) + ...
         sqrt(pi).*( exp(r.^4).*(10-11.*r.^2+2.*r.^4) - 2.*exp(r.^2).*(r.^2).*(3-7.*r.^4+2.*r.^8)).*( L.*erfc(L/2) - (L+2.*t).*erfc(L/2+t) - (L-2*z) .* erfc(L/2-z) + (L+2*t-2*z).*erfc(L/2 + t - z) ) );
%CHECKS OUT
tic

display('Calculating derivatives....')
dErSphi = (1/r)*diff(ErS, phi);  
dErSz   = diff(ErS, z);
   
dEphiSr = (1/r)*diff(EphiS*r, r);
dEphiSz   = diff(EphiS, z);

dEzSphi = (1/r)*diff(EzS, phi);
dEzSr = diff(EzS, r);

toc

tic 
display ('Creating matlab functions...')
dErphi = matlabFunction(dErSphi);
dErz   = matlabFunction(dErSz);

dEphir = matlabFunction(dEphiSr);
dEphiz   = matlabFunction(dEphiSz);

dEzphi = matlabFunction(dEzSphi);
dEzr   = matlabFunction(dEzSr);


toc
clear r phi z t L mu eps

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
close all
figure('outerposition', [1 50 1920 1150]) 
h = axes;
set(h, 'cameraposition', [0.5 0.5 9.16025], ...
       'cameraviewangle', 8, ...
       'cameratarget', [gridR/2 gridPHI/2, gridZ/2],...
       'DataAspectRatio', [1 1 1],...
       'View', [-30 20], ...
       'XLim', [-10 gridR], ...
       'Ylim', [-2.28008 gridPHI], ...
       'Zlim', [0 gridZ]);
set(h, 'clim', cVar)
axis equal

for n = 1:maxIter
    
        
%% Source Condition  

% Time positions
t1  = (n)*dt;
t12 = (n-1/2)*dt;


HrExact = HrFun(r1(:,:,2:end),phi12(:,:,2:end),z12,t12,L2,muNot) ;
HzExact = HzFun(r12,phi12(2:end,:,:),z1(2:end,:,:),t12,L2,muNot) ;
HphiExact = HphiFun(r12(:,:,2:end),phi1(2:end,:,2:end),z12(2:end,:,:),t12,L2,muNot) ;

EzExact = EzFun(r1(:,:,2:end),phi1(:,:,2:end),z12,t1,L2) ;
ErExact = ErFun(r12,phi1(2:end,:,:),z1(2:end,:,:),t1,L2) ;
EphiExact = EphiFun(r1,phi12,z1,t1,L2) ;


JzExact = JzFun(r1(:,:,2:end),phi1(:,:,2:end),z12,t12,L2,muNot,epsNot) ;
JrExact = JrFun(r12,phi1(2:end,:,:),z1(2:end,:,:),t12,L2,muNot,epsNot) ;
JphiExact = JphiFun(r1,phi12,z1,t12,L2,muNot,epsNot) ;

Jz = JzExact;
Jr = JrExact;
Jphi = JphiExact;


display(strcat('Max Jphi = ', num2str( max(max(max(Jphi)) ) ) ) )
%% Field Update

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  H-field Updates  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(strcat('Computing step...', num2str(n)))
    %%%%%%%%%% Hr FIELD %%%%%%%%%%
    QEzphi = bz2.*QEzphi - gz2.*( Ephi(:,:,2:end) - Ephi(:,:,1:end-1) )./dz   ;
    Hr(2:end,:,:) = Hr(2:end,:,:) + ...
        dr1(2:end,:,:).*( (       Ephi(2:end,:,2:end)         - Ephi(2:end,:,1:end-1) )./dz + QEzphi(2:end,:,:) )   - ...
        dr1(2:end,:,:).*(   [Ez(2:end,2:end,:) Ez(2:end,1,:)] -     Ez(2:end,:,:)     )./(iN(2:end,:,2:end)*dr*dphi) ;
    
%     r=0 Singularity Update 
    Hr(1,:,:) = Hr(1,:,:) + .5* (...
                .5*dr1(1,:,:).*( (  [Ez(2,3*gridPHI/4+2:end,  :    )   Ez(2,1:3*gridPHI/4+1,   :   ) ] - [Ez(2,gridPHI/4+2:end,:), Ez(2,1:gridPHI/4+1,:)                ] )./dr + QEzphi(1,:,:) )+ ...
                .5*dr1(1,:,:).*(  - [Er(1,3*gridPHI/4+2:end,2:end  )   Er(1,1:3*gridPHI/4+1, 2:end ) ] + [Er(1,gridPHI/4+2:end,2:end)   Er(1,1:gridPHI/4+1,2:end)       ] - ...
                                    [Er(1,  gridPHI/4+2:end,1:end-1)   Er(1,1:  gridPHI/4+1,1:end-1) ] + [Er(1,3*gridPHI/4+2:end,1:end-1) Er(1,1:3*gridPHI/4+1,1:end-1) ] )./dz +...
                .5*dr1(1,:,:).*( (  [Ez(2,3*gridPHI/4+1:end,  :    )   Ez(2,1:3*gridPHI/4,   :   ) ] - [Ez(2,gridPHI/4+1:end,:), Ez(2,1:gridPHI/4,:)                ] )./dr + QEzphi(1,:,:) )+ ...
                .5*dr1(1,:,:).*(  - [Er(1,3*gridPHI/4+1:end,2:end  )   Er(1,1:3*gridPHI/4, 2:end ) ] + [Er(1,gridPHI/4+1:end,2:end)   Er(1,1:gridPHI/4,2:end)       ] - ...
                                    [Er(1,  gridPHI/4+1:end,1:end-1)   Er(1,1:  gridPHI/4,1:end-1) ] + [Er(1,3*gridPHI/4+1:end,1:end-1) Er(1,1:3*gridPHI/4,1:end-1) ] )./dz );
    
                                %%%%%%%%%% Hz FIELD %%%%%%%%%%
    QErphi(2:end,:,:) = br2(2:end,:,:).*QErphi(2:end,:,:) - gr2(2:end,:,:).*( iN(3:end,:,:).*Ephi(3:end,:,:)   - iN(2:end-1,:,:).*Ephi(2:end-1,:,:) )./( dr.*iNH(2:end,:,:) );

   Hz(1:end,:,:) = Hz(1:end,:,:) + ...
     dz1(1:end,:,:).*(   [Er(1:end,2:end,:) Er(1:end,1,:)]  -             Er(1:end,:,:)          )./( iNH(1:end,:,:)*dr*dphi )  - ...
     dz1(1:end,:,:).*( (   iN(2:end,:,:).*Ephi(2:end,:,:)   - iN(1:end-1,:,:).*Ephi(1:end-1,:,:) )./(   dr.*iNH(1:end,:,:)   ) + QErphi(1:end,:,:) );


    %%%%%%%%%% Hphi FIELD %%%%%%%%%%
    QEzr = bz2(2:end,:,:).*QEzr - gz2(2:end,:,:).*( Er(:,:,2:end) - Er(:,:,1:end-1) )./dz;
    QErz = br2(:,:,2:end).*QErz - gr2(:,:,2:end).*( Ez(2:end,:,:) - Ez(1:end-1,:,:) )./dr;
    Hphi = Hphi + ...
           dphi1.*( ( Ez(2:end,:,:) - Ez(1:end-1,:,:) )./dr + QErz )- ...
           dphi1.*( ( Er(:,:,2:end) - Er(:,:,1:end-1) )./dz + QEzr );
                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  E-field Updates  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%%%%%%%%% Er FIELD %%%%%%%%%%
    QHzphi(:,:,2:end-1) = bz1(2:end,:,2:end-1).*QHzphi(:,:,2:end-1) - gz1(2:end,:,2:end-1).*(Hphi(:,:,2:end) - Hphi(:,:,1:end-1))./dz;
    Er(:,:,2:end-1) = cr1(:,:,2:end-1).*Er(:,:,2:end-1) - ...
                      cr2(:,:,2:end-1).*( ( Hphi(:,:,2:end) -              Hphi(:,:,1:end-1)            )./dz  + QHzphi(:,:,2:end-1) ) + ...
                      cr2(:,:,2:end-1).*( Hz(:,:,2:end-1)   - [Hz(:,end,2:end-1) Hz(:,1:end-1,2:end-1)] )./( iNH(:,:,2:end-1)*dr*dphi ) - ...
                      jr(:,:,2:end-1).*Jr(:,:,2:end-1);
                  
    %%%%%%%%% Ez FIELD %%%%%%%%%%
    QHrphi(2:end-1,:,:) = br1(2:end-1,:,2:end).*QHrphi(2:end-1,:,:)- gr1(2:end-1,:,2:end).*( iNH(2:end,:,2:end).*Hphi(2:end,:,:) - iNH(1:end-1,:,2:end).*Hphi(1:end-1,:,:)  )./(dr*iN(2:end-1,:,2:end));     
    Ez(2:end-1,:,:) = cz1(2:end-1,:,:).*Ez(2:end-1,:,:) - ...
                      cz2(2:end-1,:,:).*( (           Hr(2:end-1,:,:)         -  [Hr(2:end-1,end,:) Hr(2:end-1,1:end-1,:)] )./(iN(2:end-1,:,2:end)*dr*dphi) + QHrphi(2:end-1,:,:) )+ ...
                      cz2(2:end-1,:,:).*( iNH(2:end,:,2:end).*Hphi(2:end,:,:) -   iNH(1:end-1,:,2:end).*Hphi(1:end-1,:,:)  )./(    dr*iN(2:end-1,:,2:end)     ) - ...
                      jz(2:end-1,:,:).*Jz(2:end-1,:,:);
                  
    % r=0 Singularity Update 
    HphiC = zeros(1,1,gridZ-1);
    HphiCU = zeros(1,gridPHI,gridZ-1);
    for k = 1:gridPHI;
        HphiC(1,1,:) = HphiC(1,1,:)+Hphi(1,k,:);
    end
    for k = 1:gridPHI;
        HphiCU(1,k,:) = HphiC;
    end
    
    Ez(1,:,:) = cz1(1,:,:).*Ez(1,:,:) + ...
                cz2(1,:,:).*( 4/(dr*gridPHI) ).*HphiCU - ...
                jz(1,:,:).*Jz(1,:,:);
              
    %%%%%%%%%% Ephi FIELD %%%%%%%%%%
    QHzr(2:end-1,:,2:end-1) = bz1(2:end-1,:,2:end-1).*QHzr(2:end-1,:,2:end-1) - gz1(2:end-1,:,2:end-1).*( Hr(2:end-1,:,2:end) - Hr(2:end-1,:,1:end-1) )./dz ; 
    
    QHrz(2:end-1,:,2:end-1) = br1(2:end-1,:,2:end-1).*QHrz(2:end-1,:,2:end-1) - gr1(2:end-1,:,2:end-1).*( Hz(2:end,:,2:end-1) - Hz(1:end-1,:,2:end-1) )./dr;
    Ephi(2:end-1,:,2:end-1) = cphi1(2:end-1,:,2:end-1).*Ephi(2:end-1,:,2:end-1) + ...
                              cphi2(2:end-1,:,2:end-1).*( ( Hr(2:end-1,:,2:end) - Hr(2:end-1,:,1:end-1) )./dz + QHzr(2:end-1,:,2:end-1) ) - ...
                              cphi2(2:end-1,:,2:end-1).*( ( Hz(2:end,:,2:end-1) - Hz(1:end-1,:,2:end-1) )./dr + QHrz(2:end-1,:,2:end-1) ) - ...
                              jphi(2:end-1,:,2:end-1).*Jphi(2:end-1,:,2:end-1);            
    % r=0 Singularity Update                  
    Ephi(1,:,2:end-1) = cphi1(1,:,2:end-1).*Ephi(1,:,2:end-1) + ...
                        cphi2(1,:,2:end-1).*( Hr(1,:,2:end) - Hr(1,:,1:end-1) )./dz - ...
                        cphi2(1,:,2:end-1).*( Hz(1,:,2:end-1) -[ Hz(1,gridPHI/2:end,2:end-1) Hz(1,1:gridPHI/2-1,2:end-1)] )./dr - ...
                        jphi(1,:,2:end-1).*Jphi(1,:,2:end-1);    
                    
    %%%%%%%%%% Boundary Update %%%%%%%%%%
    Er  ( :  ,: ,1   ) = 0;
    Er  ( :  ,: ,end ) = 0;
    Ez  ( end,: ,:   ) = 0;
    Ephi( :  ,: ,1   ) = 0;
    Ephi( :  ,: ,end ) = 0;
    Ephi( end,: ,:   ) = 0;
    

%% Continuous plotting
    planePos = 25;
    planePos2 = 70;
    
    phiPlane =8;
    phiPlane2 = 18;

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comparison = Hz(:,:,round(gridZ/2))-HzExact(:,:,round(gridZ/2));
error(counter) =  norm(comparison);
maxs(counter) = max(max(max(comparison)));
averages(counter) = sum(sum(sum(comparison)))/(gridZ*gridR*gridPHI);
end

%%
close all

figure
plot( log10(totalPoints), log10(error./points.^2), 'g+', log10(totalPoints), log10(maxs), 'b*');
% log10(points), log10(averages), 'ro')

