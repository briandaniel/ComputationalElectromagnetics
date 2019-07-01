clear all;
close all;

%% Program Info  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Simulates FDTD in 3D Cartesian Coordinates                              %
% Utilizes CPML Boundary                                                  %
%                                                                         %
%  This Source Code Form is subject to the terms of the Mozilla Public    %
%  License, v. 2.0. If a copy of the MPL was not distributed with this    %
%  file, You can obtain one at http://mozilla.org/MPL/2.0/.               %
%                                                                         %
%  Copyright (C) 2012-2013 Brian D. Hong                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% Physical Constant
epsNot = 8.854187817620e-12;    % Permittivity of Free Space
muNot = 1.25663706e-6;          % Magnetic Constant
cNot = 2.99792458e8;            % Speed of light
etaNot = sqrt(muNot/epsNot);    % Characteristic impedence of free space

%% Parameters
% Define FDTD physical space
domain = [3e-6 3e-6 3e-6] ;     % Computational Domain [X Y Z]             [USER DEFINED]
gridC  = [50   50    50] ;     % Number of gridpoints [X Y Z]             [USER DEFINED]
dx = domain(1)/gridC(1);        % Cell width  
dy = domain(2)/gridC(2);        % Cell length  
dz = domain(3)/gridC(3);        % Cell height
dmin = min([dx dy dz]) ;        % Minimum cell size
% Define FDTD temporal space
maxIter = max(gridC)*6;             % Maximum number of iterations         [USER DEFINED]
dt = .5*dmin/(cNot*sqrt(3));   % Time step (1 = 1sec)                 [USER DEFINED]
Scr = cNot*dt/dmin;    % Courant Number

% Define positions
planePosZ = .5e-6;                  % Position of cross sectional viewing plane
planePosZ = round(planePosZ/dz);
sourcePos = [30 30, 30 30, 30 30];  % Source position 
                                    % [xMin xMax, yMin yMax, zMin zMax]    [USER DEFINED]
                                    
zPlaneNum = 8;                      % Number of Zplanes                                    
cVar = [-1 1];                      % Color variation in plot
skipFrames = 30;

% Surface view display setting
surfaceView = 3;
% 0 = Total E Field
% 1 = Ex
% 1.5 = Ex Middle Plane
% 2 = Ey
% 3 = Ez
% 3.1 = Ez contour (Not working)
% 3.5 = Ez Middle Plane
% 4 = Hx
% 5 = Hy
% 6 = Hz

% Define PML Parameters
% Display PML parameter graphs and stop before FDTD simulation?
displayC = 0;    % 0 if no, 1 if yes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PML Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmlWidth     = 10;              % Width of PML in grid points              [USER DEFINED]
sigM         = 3;               % Power of sigma in PML                    [USER DEFINED]
sigMax       = (sigM+1)*...   % Max value of sigma in PML
              .8/(etaNot*dmin);
kappaMax     = 1;               % Max value of kappa in PML                [USER DEFINED]
kappaM       = 3;               % Power of kappa in PML                    [USER DEFINED]
aMax         = 1;               % Max value of kappa in PML                [USER DEFINED]
aM           = 3;               % Power of a in PML                        [USER DEFINED]

% Initial Pulse Parameters
% Gaussian Pulse parameters
width = 5;                     % Width of Gaussian                        [USER DEFINED]
maxT = width*6;
magG = 150;

% Blackman-Harris Parameters
a = [0.35875 0.48829 0.14128 0.01168];
widthB = 30*dt;
mag = 1;

% Sine Wave Parameters
f = 2501e12;                     % Frequency                                [USER DEFINED]
amp = 1;                        % Amplitude of Wave                        [USER DEFINED]
lambda =cNot/f;                 % Wavelength
omega = 2*pi*f;                 % Angular Frequency

%% Initialize solution vectors
% General Cartesian coordinates: f(x,y,z)

% Magnetic Field Vectors
Hx = zeros(gridC(1)-1, gridC(2)  , gridC(3)  );
Hy = zeros(gridC(1)  , gridC(2)-1, gridC(3)  );
Hz = zeros(gridC(1)  , gridC(2)  , gridC(3)-1);

% Electric Field Vectors
Ex = zeros(gridC(1)  , gridC(2)+1, gridC(3)+1);
Ey = zeros(gridC(1)+1, gridC(2)  , gridC(3)+1);
Ez = zeros(gridC(1)+1, gridC(2)+1, gridC(3)  );

% Magnetic Source
Mx = zeros(gridC(1)-1, gridC(2)  , gridC(3)  );
My = zeros(gridC(1)  , gridC(2)-1, gridC(3)  );
Mz = zeros(gridC(1)  , gridC(2)  , gridC(3)-1);

% Current Source
Jx = zeros(gridC(1)  , gridC(2)+1, gridC(3)+1);
Jy = zeros(gridC(1)+1, gridC(2)  , gridC(3)+1);
Jz = zeros(gridC(1)+1, gridC(2)+1, gridC(3)  );

% Sigma Star
sigS = zeros(gridC(1)  , gridC(2)  , gridC(3)  );
% sigSx = zeros(gridC(1)-1, gridC(2)  , gridC(3)  );
% sigSy = zeros(gridC(1)  , gridC(2)-1, gridC(3)  );
% sigSz = zeros(gridC(1)  , gridC(2)  , gridC(3)-1);

% Sigma
sig = zeros(gridC(1)+1 , gridC(2)+1, gridC(3)+1);
% sigX = zeros(gridC(1)  , gridC(2)+1, gridC(3)+1);
% sigY = zeros(gridC(1)+1, gridC(2)  , gridC(3)+1);
% sigZ = zeros(gridC(1)+1, gridC(2)+1, gridC(3)  );

% Relative permeability
mu  = ones(gridC(1)  , gridC(2)  , gridC(3)  );
% muX = ones(gridC(1)-1, gridC(2)  , gridC(3)  );
% muY = ones(gridC(1)  , gridC(2)-1, gridC(3)  );
% muZ = ones(gridC(1)  , gridC(2)  , gridC(3)-1);

% Relative permittivity
eps = ones(gridC(1)+1 , gridC(2)+1, gridC(3)+1);
% epsX = ones(gridC(1)  , gridC(2)+1, gridC(3)+1);
% epsY = ones(gridC(1)+1, gridC(2)  , gridC(3)+1);
% epsZ = ones(gridC(1)+1, gridC(2)+1, gridC(3)  );

% Vectors for use in PML updates
% QH for updating Efields
QEyz = zeros(gridC(1)-1, gridC(2)  , gridC(3)  );
QEzy = zeros(gridC(1)-1, gridC(2)  , gridC(3)  );

QEzx = zeros(gridC(1)  , gridC(2)-1, gridC(3)  );
QExz = zeros(gridC(1)  , gridC(2)-1, gridC(3)  ); 

QEyx = zeros(gridC(1)  , gridC(2)  , gridC(3)-1);
QExy = zeros(gridC(1)  , gridC(2)  , gridC(3)-1);

% QE for updating Hfields
QHyz = zeros(gridC(1)  , gridC(2)+1, gridC(3)+1);
QHzy = zeros(gridC(1)  , gridC(2)+1, gridC(3)+1);

QHzx = zeros(gridC(1)+1, gridC(2)  , gridC(3)+1);
QHxz = zeros(gridC(1)+1, gridC(2)  , gridC(3)+1);

QHyx = zeros(gridC(1)+1, gridC(2)+1, gridC(3)  );
QHxy = zeros(gridC(1)+1, gridC(2)+1, gridC(3)  );

% Meshgrid for plots
[XX YY] = meshgrid(1:gridC(1), 1:gridC(2));

%% Compute coefficient vectors
% E-Field coeff.
c1 = (1-(sig.*dt)./(2*epsNot.*eps))./(1+(sig.*dt)./(2*epsNot.*eps));
c2 = (dt./(epsNot.*eps))./(1+(sig.*dt)./(2.*epsNot.*eps));

% D-Field coeff.
d1 = (1-(sigS.*dt)./(2*muNot.*mu))./(1+(sigS.*dt)./(2*muNot.*mu));
d2 = (dt./(muNot.*mu))./(1+(sigS.*dt)./(2*muNot.*mu));


%% Compute PML Coefficients

% Compute X direction
xVec = (1:pmlWidth)-.5;
gridPlot = (1:gridC(1)+1)-1;
gridPlotStag = (1:gridC(1))-.5;

sigVec = (abs(xVec+.5).^sigM./pmlWidth.^sigM).*sigMax;
sigVecStag = (abs(xVec).^sigM./pmlWidth.^sigM).*sigMax;

sigVec = [sigVec(end:-1:1) zeros(1,gridC(1)-2*length(sigVec)+1) sigVec];
sigVecStag = [sigVecStag(end:-1:1) zeros(1,gridC(1)-2*length(sigVecStag)) sigVecStag];

kappaVec = 1+(abs(xVec+.5).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVec = [kappaVec(end:-1:1) ones(1,gridC(1)-2*length(kappaVec)+1) kappaVec];
kappaVecStag = 1+(abs(xVec).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVecStag = [kappaVecStag(end:-1:1) ones(1,gridC(1)-2*length(kappaVecStag)) kappaVecStag];

aVec = (abs(xVec+.5).^aM./pmlWidth.^aM).*aMax;
aVec = [aVec(end:-1:1) zeros(1,gridC(1)-2*length(aVec)+1) aVec];
aVecStag = (abs(xVec).^aM./pmlWidth.^aM).*aMax;
aVecStag = [aVecStag(end:-1:1) zeros(1,gridC(1)-2*length(aVecStag)) aVecStag];

% Plots of vectors used in PML
if displayC == 1;
figure('outerposition', [100, 100, 1500, 900])
subplot(2,3,1)
plot(gridPlot, sigVec, 'ro', gridPlotStag, sigVecStag, 'b+')
title('PML {\sigma}')

subplot(2,3,2)
plot(gridPlot, kappaVec, 'ro', gridPlotStag, kappaVecStag, 'b+')
title('PML {\kappa}')

subplot(2,3,3)
plot(gridPlot, aVec, 'ro', gridPlotStag, aVecStag, 'b+')
title('PML a')
end

% Generate coefficient vectors used for PML during computation
tau = (kappaVec.*epsNot)./(kappaVec.*aVec+sigVec);
tauStag = (kappaVecStag.*epsNot)./(kappaVecStag.*aVecStag+sigVecStag);
bVec = exp(-dt./tau);
bVecStag = exp(-dt./tauStag);
gVec = (sigVec./(kappaVec.*(kappaVec.*aVec+sigVec))).*(1-exp(-dt./tau));
gVecStag = (sigVecStag./(kappaVecStag.*(kappaVecStag.*aVecStag+sigVecStag))).*(1-exp(-dt./tauStag));

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

if displayC == 1
subplot(2,3,4)
plot(gridPlot, bVec, 'ro', gridPlotStag, bVecStag, 'b+')
title('coef b')

subplot(2,3,5)
plot(gridPlot, gVec, 'ro', gridPlotStag, gVecStag, 'b+')
title('coef g')
end
bxIn  = repmat(bVecStag',[1, gridC(2),gridC(3)]);
bxOut = repmat(bVec',[1, gridC(2)+1,gridC(3)+1]);
gxIn  = repmat(gVecStag',[1, gridC(2),gridC(3)]);
gxOut = repmat(gVec',[1, gridC(2)+1,gridC(3)+1]);
kxIn = repmat(1./kappaVecStag',[1, gridC(2),gridC(3)]);
kxOut = repmat(1./kappaVec',[1,gridC(2)+1,gridC(3)+1]);


if displayC == 1;
    figure('outerposition', [100 100 1400 900])
    h = axes;
    hold on
    
    [X Y] = meshgrid(gridPlot, 1:gridC(2)+1);
    [X2 Y2] = meshgrid(gridPlotStag, 1:gridC(2));
    for k = 1:round(gridC(3)/zPlaneNum):gridC(3);
        surf( X', Y', k.*ones(size(X')), bxOut(1:gridC(1), 1:gridC(2)+1, k), 'edgealpha', (1-(k./gridC(3)).^(1/2))/2, 'facealpha', (1-(k./gridC(3)).^(1/2)))
    display(num2str((1-(k./gridC(3)).^2)))
    end
    
    for k = 1:round(gridC(3)/zPlaneNum):gridC(3);
        surf( X2', Y2', k.*ones(size(X2'))-1/2, bxIn(1:gridC(1)-1, 1:gridC(2), k), 'edgealpha',(1-(k./gridC(3)).^(1/2))/2, 'facealpha', (1-(k./gridC(3)).^(1/2)))
    end
    
    set(h, 'cameraposition', [0.5 0.5 9.16025], ...
           'cameraviewangle', 5, ...
           'cameratarget', [gridC(1)/2 gridC(2)/2, gridC(3)/2],...
           'DataAspectRatio', [1 1 1],...
           'View', [-15 20], ...
           'XLim', [-11.5899 48.9016], ...
           'Ylim', [-2.28008 37.7199], ...
           'Zlim', [0 20]);
           
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    hold off
    axis equal
title('bxIn')
end

%% Compute Y direction
xVec = (1:pmlWidth)-.5;
gridPlot = (1:gridC(2)+1)-1;
gridPlotStag = (1:gridC(2))-.5;

sigVec = (abs(xVec+.5).^sigM./pmlWidth.^sigM).*sigMax;
sigVecStag = (abs(xVec).^sigM./pmlWidth.^sigM).*sigMax;

sigVec = [sigVec(end:-1:1) zeros(1,gridC(2)-2*length(sigVec)+1) sigVec];
sigVecStag = [sigVecStag(end:-1:1) zeros(1,gridC(2)-2*length(sigVecStag)) sigVecStag];

kappaVec = 1+(abs(xVec+.5).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVec = [kappaVec(end:-1:1) ones(1,gridC(2)-2*length(kappaVec)+1) kappaVec];
kappaVecStag = 1+(abs(xVec).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVecStag = [kappaVecStag(end:-1:1) ones(1,gridC(2)-2*length(kappaVecStag)) kappaVecStag];

aVec = (abs(xVec+.5).^aM./pmlWidth.^aM).*aMax;
aVec = [aVec(end:-1:1) zeros(1,gridC(2)-2*length(aVec)+1) aVec];
aVecStag = (abs(xVec).^aM./pmlWidth.^aM).*aMax;
aVecStag = [aVecStag(end:-1:1) zeros(1,gridC(2)-2*length(aVecStag)) aVecStag];

% Generate coefficient vectors used for PML during computation
tau = (kappaVec.*epsNot)./(kappaVec.*aVec+sigVec);
tauStag = (kappaVecStag.*epsNot)./(kappaVecStag.*aVecStag+sigVecStag);
bVec = exp(-dt./tau);
bVecStag = exp(-dt./tauStag);
gVec = (sigVec./(kappaVec.*(kappaVec.*aVec+sigVec))).*(1-exp(-dt./tau));
gVecStag = (sigVecStag./(kappaVecStag.*(kappaVecStag.*aVecStag+sigVecStag))).*(1-exp(-dt./tauStag));

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

byIn  = repmat(bVecStag,       [gridC(1)  , 1, gridC(3)  ]);
byOut = repmat(bVec,           [gridC(1)+1, 1 ,gridC(3)+1]);
gyIn  = repmat(gVecStag,       [gridC(1)  , 1, gridC(3)  ]);
gyOut = repmat(gVec,           [gridC(1)+1, 1 ,gridC(3)+1]);
kyIn  = repmat(1./kappaVecStag,[gridC(1)  , 1, gridC(3)  ]);
kyOut = repmat(1./kappaVec,    [gridC(1)+1, 1 ,gridC(3)+1]);

%% Compute Z direction
xVec = (1:pmlWidth)-.5;
gridPlot = (1:gridC(3)+1)-1;
gridPlotStag = (1:gridC(3))-.5;

sigVec = (abs(xVec+.5).^sigM./pmlWidth.^sigM).*sigMax;
sigVecStag = (abs(xVec).^sigM./pmlWidth.^sigM).*sigMax;

sigVec = [sigVec(end:-1:1) zeros(1,gridC(3)-2*length(sigVec)+1) sigVec];
sigVecStag = [sigVecStag(end:-1:1) zeros(1,gridC(3)-2*length(sigVecStag)) sigVecStag];

kappaVec = 1+(abs(xVec+.5).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVec = [kappaVec(end:-1:1) ones(1,gridC(3)-2*length(kappaVec)+1) kappaVec];
kappaVecStag = 1+(abs(xVec).^kappaM./pmlWidth.^kappaM).*(kappaMax-1);
kappaVecStag = [kappaVecStag(end:-1:1) ones(1,gridC(3)-2*length(kappaVecStag)) kappaVecStag];

aVec = (abs(xVec+.5).^aM./pmlWidth.^aM).*aMax;
aVec = [aVec(end:-1:1) zeros(1,gridC(3)-2*length(aVec)+1) aVec];
aVecStag = (abs(xVec).^aM./pmlWidth.^aM).*aMax;
aVecStag = [aVecStag(end:-1:1) zeros(1,gridC(3)-2*length(aVecStag)) aVecStag];

% Generate coefficient vectors used for PML during computation
tau = (kappaVec.*epsNot)./(kappaVec.*aVec+sigVec);
tauStag = (kappaVecStag.*epsNot)./(kappaVecStag.*aVecStag+sigVecStag);
bVec = exp(-dt./tau);
bVecStag = exp(-dt./tauStag);
gVec = (sigVec./(kappaVec.*(kappaVec.*aVec+sigVec))).*(1-exp(-dt./tau));
gVecStag = (sigVecStag./(kappaVecStag.*(kappaVecStag.*aVecStag+sigVecStag))).*(1-exp(-dt./tauStag));

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

bVecStagZ = ones(1,1,length(bVecStag));
bVecStagZ(1,1,:) = bVecStag;

bVecZ = ones(1,1,length(bVec));
bVecZ(1,1,:) = bVec;

gVecStagZ = ones(1,1,length(gVecStag));
gVecStagZ(1,1,:) = gVecStag;

gVecZ = ones(1,1,length(gVec));
gVecZ(1,1,:) = gVec;

kappaVecStagZ = ones(1,1,length(kappaVecStag));
kappaVecStagZ(1,1,:) = kappaVecStag;

kappaVecZ = ones(1,1,length(kappaVec));
kappaVecZ(1,1,:) = kappaVec;

bzIn  = repmat(bVecStagZ,       [gridC(1)  , gridC(2)  , 1]);
bzOut = repmat(bVecZ,           [gridC(1)+1, gridC(2)+1, 1]);
gzIn  = repmat(gVecStagZ,       [gridC(1)  , gridC(2)  , 1]);
gzOut = repmat(gVecZ,           [gridC(1)+1, gridC(2)+1, 1]);
kzIn  = repmat(1./kappaVecStagZ,[gridC(1)  , gridC(2)  , 1]);
kzOut = repmat(1./kappaVecZ,    [gridC(1)+1, gridC(2)+1, 1]);

if displayC == 1;
    return
end

%%
display(' ')
display('Memory status after creation of solution vectors')


close all
figure('outerposition', [100 50 800 800]) 
h = axes;
set(h,'DataAspectRatio', [1 1 1],...
       'View', [-30 4], ...
       'XLim', [-10 gridC(1)], ...
       'Ylim', [-2.28008 gridC(2)], ...
       'Zlim', [0 gridC(3)]);
set(h, 'clim', cVar)
axis equal
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN SOLUTION LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for n = 1:maxIter
    
% Source update
%     source = exp(-.5*((maxT*dmin-n*dt*cNot/Scr)/(width*dmin))^2) *amp*sin(omega*n*dt)
      source = magG*exp(-.5*((maxT*dmin-n*dt*cNot/Scr)/(width*dmin))^2);
%     source = amp*sin(omega*n*dt);

%     Jx(sourcePos(1):sourcePos(2), sourcePos(3):sourcePos(4), sourcePos(5):sourcePos(6)) = source;
%     Jy(sourcePos(1):sourcePos(2), sourcePos(3):sourcePos(4), sourcePos(5):sourcePos(6)) = source;
%     Jz(sourcePos(1):sourcePos(2), sourcePos(3):sourcePos(4), sourcePos(5):sourcePos(6)) = source;
%     Ez(sourcePos(1):sourcePos(2), sourcePos(3):sourcePos(4), sourcePos(5):sourcePos(6)) = source;
      Ez(sourcePos(1):sourcePos(2), sourcePos(3):sourcePos(4), sourcePos(5):sourcePos(6)) = source;
%       Ey(sourcePos(1):sourcePos(2), sourcePos(3):sourcePos(4), sourcePos(5):sourcePos(6)) = source;
%       Ex(sourcePos(1):sourcePos(2), sourcePos(3):sourcePos(4), sourcePos(5):sourcePos(6)) = source;
%%%%%%%%%%%%%%%%%%%%% Cartesian Coordinates: f(x,y,z) %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  H-field Updates  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%% Hx FIELD %%%%%%%%%%
    QEzy = bzIn(2:end,:,:).*QEzy - gzIn(2:end,:,:).*(Ey(2:end-1,:    ,2:end) - Ey(2:end-1, :      , 1:end-1) )./dz; 
    QEyz = byIn(2:end,:,:).*QEyz - gyIn(2:end,:,:).*(Ez(2:end-1,2:end,:    ) - Ez(2:end-1, 1:end-1, :      ) )./dy; 
    Hx = d1(2:end,:,:).*Hx + ...
         d2(2:end,:,:).*( kzIn(2:end,:,:).*(Ey(2:end-1,:    ,2:end) - Ey(2:end-1, :      , 1:end-1) )./dz - ...
                          kyIn(2:end,:,:).*(Ez(2:end-1,2:end,:    ) - Ez(2:end-1, 1:end-1, :      ) )./dy - ...
                          Mx + QEzy - QEyz);

    %%%%%%%%%% Hy FIELD %%%%%%%%%%
    QExz = bxIn(:,2:end,:).*QExz - gxIn(:,2:end,:).*(Ez(2:end,2:end-1, :   ) - Ez(1:end-1,2:end-1, :)      )./dx;
    QEzx = bzIn(:,2:end,:).*QEzx - gzIn(:,2:end,:).*(Ex(:    ,2:end-1,2:end) - Ex(:      ,2:end-1,1:end-1) )./dz; 

    Hy = d1(:,2:end,:).*Hy + ...
         d2(:,2:end,:).*( kxIn(:,2:end,:).*(Ez(2:end,2:end-1, :   ) - Ez(1:end-1,2:end-1, :     ) )./dx - ...
                          kzIn(:,2:end,:).*(Ex(:    ,2:end-1,2:end) - Ex(:      ,2:end-1,1:end-1) )./dz - ...
                          My + QExz - QEzx);                
                      
    %%%%%%%%%% Hz FIELD %%%%%%%%%%
    QEyx = byIn(:,:,2:end).*QEyx - gyIn(:,:,2:end).*(Ex(:    ,2:end  ,2:end-1) - Ex(:      ,2:end  ,2:end-1) )./dy;
    QExy = bxIn(:,:,2:end).*QExy - gxIn(:,:,2:end).*(Ey(2:end, :     ,2:end-1) - Ey(1:end-1, :     ,2:end-1) )./dx; 

    Hz = d1(:,:,2:end).*Hz + ...
         d2(:,:,2:end).*( kyIn(:,:,2:end).*(Ex(:    ,2:end  ,2:end-1) - Ex(:      ,2:end  ,2:end-1) )./dy - ...
                          kxIn(:,:,2:end).*(Ey(2:end, :     ,2:end-1) - Ey(1:end-1, :     ,2:end-1) )./dx - ...
                          Mz + QEyx - QExy);       

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  E-field Updates  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%% Ex FIELD %%%%%%%%%%
    QHyz(:,2:end-1,2:end-1) = byOut(2:end,2:end-1,2:end-1).*QHyz(:,2:end-1,2:end-1) - gyOut(2:end,2:end-1,2:end-1).*( Hz(:,2:end, :   ) - Hz(:,1:end-1, :     ) )./dy;     
    QHzy(:,2:end-1,2:end-1) = bzOut(2:end,2:end-1,2:end-1).*QHzy(:,2:end-1,2:end-1) - gzOut(2:end,2:end-1,2:end-1).*( Hy(:, :   ,2:end) - Hy(:, :     ,1:end-1) )./dz;
    
    Ex(:,2:end-1,2:end-1) = c1(2:end,2:end-1,2:end-1).*Ex(:,2:end-1,2:end-1) + ...
                            c2(2:end,2:end-1,2:end-1).*( kyOut(2:end,2:end-1,2:end-1).*( Hz(:,2:end, :   ) - Hz(:,1:end-1, :     ) )./dy - ...
                                                         kzOut(2:end,2:end-1,2:end-1).*( Hy(:, :   ,2:end) - Hy(:, :     ,1:end-1) )./dz - ...
                                                         Jx(:,2:end-1,2:end-1) + QHyz(:,2:end-1,2:end-1) - QHzy(:,2:end-1,2:end-1));
  
    %%%%%%%%%% Ey FIELD %%%%%%%%%%
    QHzx(2:end-1,:,2:end-1) = bzOut(2:end-1,2:end,2:end-1).*QHzx(2:end-1,:,2:end-1) - gzOut(2:end-1,2:end,2:end-1).*(Hx(:    ,:,2:end) - Hx(:      ,:,1:end-1) )./dz;     
    QHxz(2:end-1,:,2:end-1) = bxOut(2:end-1,2:end,2:end-1).*QHxz(2:end-1,:,2:end-1) - gxOut(2:end-1,2:end,2:end-1).*(Hz(2:end,:, :   ) - Hz(1:end-1,:,:      ) )./dx; 
    
    Ey(2:end-1,:,2:end-1) = c1(2:end-1,2:end,2:end-1).*Ey(2:end-1,:,2:end-1) + ...
                            c2(2:end-1,2:end,2:end-1).*( kzOut(2:end-1,2:end,2:end-1).*( Hx(:    ,:,2:end) - Hx(:      ,:,1:end-1) )./dz - ...                            
                                                         kxOut(2:end-1,2:end,2:end-1).*( Hz(2:end,:, :   ) - Hz(1:end-1,:,:      ) )./dx - ...
                                                         Jy(2:end-1,:,2:end-1) + QHzx(2:end-1,:,2:end-1) - QHxz(2:end-1,:,2:end-1));

    %%%%%%%%%% Ez FIELD %%%%%%%%%%
    QHxy(2:end-1,2:end-1,:) = bxOut(2:end-1,2:end-1,2:end).*QHxy(2:end-1,2:end-1,:) - gxOut(2:end-1,2:end-1,2:end).*( Hy(2:end, :   ,:) - Hy(1:end-1, :     ,:) )./dx;    
    QHyx(2:end-1,2:end-1,:) = byOut(2:end-1,2:end-1,2:end).*QHyx(2:end-1,2:end-1,:) - gyOut(2:end-1,2:end-1,2:end).*( Hx( :   ,2:end,:) - Hx( :     ,1:end-1,:) )./dy;    
    
    Ez(2:end-1,2:end-1,:) = c1(2:end-1,2:end-1,2:end).*Ez(2:end-1,2:end-1,:) + ...        
                            c2(2:end-1,2:end-1,2:end).*( kxOut(2:end-1,2:end-1,2:end).*( Hy(2:end, :   ,:) - Hy(1:end-1, :     ,:) )./dx - ...                            
                                                         kyOut(2:end-1,2:end-1,2:end).*( Hx( :   ,2:end,:) - Hx( :     ,1:end-1,:) )./dy - ...
                                                         Jz(2:end-1,2:end-1,:) + QHxy(2:end-1,2:end-1,:) - QHyx(2:end-1,2:end-1,:));        
                                                     
%%%%%%%%%% Boundary Update %%%%%%%%%%
    Ez(1  ,:  ,:  ) = 0;
    Ez(end,:  ,:  ) = 0;
    Ez(:  ,1  ,:  ) = 0;
    Ez(:  ,end,:  ) = 0;
    
    Ey(1  ,:  ,:  ) = 0;
    Ey(end,:  ,:  ) = 0;
    Ey(:  ,:  ,1  ) = 0;
    Ey(:  ,:  ,end) = 0;
    
	Ex(:  ,end,:  ) = 0;
    Ex(:  ,1  ,:  ) = 0;
    Ex(:  ,:  ,1  ) = 0;
    Ex(:  ,:  ,end) = 0;
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT  UTILITIES  %%%%%%%%%%%%%%%%%%%%%%%% %%


display(strcat('n = ', num2str(n)))

if mod(n,skipFrames) == 0
if surfaceView == 0
    
    % Plot Total Efield
    hold on
    for k = pmlWidth:round(gridC(3)/zPlaneNum):gridC(3)-pmlWidth;
        surf( XX', YY', k.*ones(size(XX')) + Ex(XX(1,:), YY(:,1), k) + Ey(XX(1,:), YY(:,1), k) + Ez(XX(1,:), YY(:,1), k), Ex(XX(1,:), YY(:,1), k) + Ey(XX(1,:), YY(:,1), k) + Ez(XX(1,:), YY(:,1), k), 'edgealpha', (1-(k./gridC(3)).^(1/2))/2, 'facealpha', (1-(k./gridC(3)).^(1/2)))
%         surf( XX', YY', k.*ones(size(XX')) + Ex(XX(1,:), YY(:,1), k) + Ey(XX(1,:), YY(:,1), k) + Ez(XX(1,:), YY(:,1), k), Ex(XX(1,:), YY(:,1), k) + Ey(XX(1,:), YY(:,1), k) + Ez(XX(1,:), YY(:,1), k), 'edgealpha', .1, 'facealpha', .05)
    end
    hold off 
    set(h, 'clim', cVar)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    colorbar('East')
    title('Total E Field', 'fontsize', 18) 
    
end

if surfaceView == 1;
    % Plot Ex field

    hold on
    for k = pmlWidth:round(gridC(3)/zPlaneNum):gridC(3)-pmlWidth;
        surf( XX', YY', k.*ones(size(XX'))+ Ex(XX(1,:), YY(:,1), k), Ex(XX(1,:), YY(:,1), k), 'edgealpha', (1-(k./gridC(3)).^(1/2))/2, 'facealpha', (1-(k./gridC(3)).^(1/2)))
    end
    hold off 
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    colorbar('East')
    title('Ex', 'fontsize', 18)
    
elseif surfaceView == 1.5;
    surf(fliplr(Ex(:,:,round(gridC(3)/2))), 'edgealpha', .1)
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    title('Ex', 'fontsize', 18)
    set(h, 'clim', cVar)
    axis([0 gridC(1) 0 gridC(2) -1 1])

elseif surfaceView == 2;
    % Plot Ez field

    hold on
    for k = pmlWidth:round(gridC(3)/zPlaneNum):gridC(3)-pmlWidth;
        surf( XX', YY', k.*ones(size(XX'))+ Ey(XX(1,:), YY(:,1), k), Ey(XX(1,:), YY(:,1), k), 'edgealpha', (1-(k./gridC(3)).^(1/2))/2, 'facealpha', (1-(k./gridC(3)).^(1/2)))
    end
    hold off 
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    colorbar('East')
    title('Ey', 'fontsize', 18)
       
elseif surfaceView == 3;
    % Plot Ez field

    hold on
    for k = pmlWidth:round(gridC(3)/zPlaneNum):gridC(3)-pmlWidth;
        surf( XX', YY', k.*ones(size(XX'))+ Ez(XX(1,:), YY(:,1), k), Ez(XX(1,:), YY(:,1), k), 'edgealpha', (1-(k./gridC(3)).^(1/2))/2, 'facealpha', (1-(k./gridC(3)).^(1/2)))
    end
    hold off 
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    colorbar('East')
    title('Ez', 'fontsize', 18)
    
elseif surfaceView == 3.1;
    % Plot Ez field

    hold on
    for k = pmlWidth:round(gridC(3)/zPlaneNum):gridC(3)-pmlWidth;
        contour3( XX', YY', Ez(XX(1,:), YY(:,1), k), Ez(XX(1,:), YY(:,1), k))
    end
    hold off 
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    colorbar('East')
    title('Ez', 'fontsize', 18)
        
elseif surfaceView == 3.5;
    surf(fliplr(Ez(:,:,round(gridC(3)/2))), 'edgealpha', .1)
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    title('Ez', 'fontsize', 18)
    set(h, 'clim', cVar)
    axis([0 gridC(1) 0 gridC(2) -1 1])
end

    s = getframe;
    if n < maxIter
    cla
    end
end
end












    
  
    