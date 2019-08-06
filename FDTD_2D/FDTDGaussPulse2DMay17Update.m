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

%% Parameters
% Define FDTD Space
dx = .01   ;                    % Cell Size (1 = 1meter)                   [USER DEFINED]
dt = 1.67e-11 ;                 % Time step (1 = 1sec)                    [USER DEFINED]
Cn = cNot*dt/dx;                % Courant Number (only works in Vacuum)
planeSize = 250;                % Total number of cells = planeSize^2      [USER DEFINED]


% Define vectors
xDom = 1:planeSize;
yDom = xDom;
[xDom yDom] = meshgrid(xDom, yDom);

Dz = zeros(planeSize,planeSize);
Ez = zeros(planeSize,planeSize);
Hy = zeros(planeSize,planeSize-1);
Hx = zeros(planeSize-1,planeSize);
EzStore = Ez;
EzStoreOld = EzStore;

%% Hard source parameters
sourcePos = [.5*planeSize, .5*planeSize];   % Location of source           [USER DEFINED]

% Gaussian Pulse parameters
width = 10;                                  % Width of Gaussian            [USER DEFINED]
maxT = 40.0;                                % Decay of Gaussian            [USER DEFINED]

% Sine Wave Parameters
f = 5e8;                                    % Frequency                    [USER DEFINED]
amp = 1;                                    % Amplitude of Wave            [USER DEFINED]
lambda =cNot/f;                             % Wavelength
omega = 2*pi*f;                             % Angular Frequency

% Ricker Wavelet Parameters
fp = 5e9;                                   % Peak Frequency               [USER DEFINED]
md = 1;                                     % Temporal Delay Multiple      [USER DEFINED]
dr = md/fp;                                 % Temporal Delay

%% Parameters associated with the dielectric 
% posMed = 200;                   % Position of mediunm                      [USER DEFINED]
% sigma = .01 ;                   % Conductivity of dielectric               [USER DEFINED]
% dielC = 4;                      % Dielectric Constant (rel. Permittivity)  [USER DEFINED]
%       
% epsDiel = dielC*epsNot;         % Permittivity of dielectric
% 
% factor = dt*sigma/(2*epsNot*dielC);
% alpha = (1-factor)/(1+factor);
% beta = .5/(dielC*(1+factor));

epsInv = 1;



%%%%%% MAIN SOLUTION LOOP %%%%%%
figure('outerposition', [710   500   780   622])


for n = 1:2000
    
% Generate Gaussian Pulse
source = exp(-.5*((maxT-n)/width)^2);

% Generate wave packet
% source = exp(-.5*((maxT-n)/width)^2) * amp*sin(omega*n*dt);
    
% Generate sine hard source
% source = amp*sin(omega*n*dt);

% Generate Ricker Wavelet
% source = (1-2*(pi*fp*(n*dt-dr))^2)*exp(-(pi*fp*(n*dt-dr))^2);

%% Without dielectric

% Boundary Update
% Dz(1,:)   = 0 ;
% Dz(end,:) = 0 ;
% Dz(:,1)   = 0 ;
% Dz(:,end) = 0 ;

% Calculate Dz Field and convert to Ez
Dz(2:end-1, 2:end-1) = Dz(2:end-1, 2:end-1) + ...
                       Cn*( Hy(2:end-1,2:end) - Hy(2:end-1,1:end-1) ) - ...
                       Cn*( Hx(2:end,2:end-1) - Hx(1:end-1,2:end-1) );
      
Ez = epsInv*Dz;
Ez(sourcePos(1), sourcePos(2)) = source;


% MUR boundary update
% Ez(end,:) = - EzStoreOld(end-1,:) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(end-1,:)+EzStoreOld(end,:)) + 2*dx/(cNot*dt+dx)*(EzStore(end,:)+EzStore(end-1,:))+...
%           (cNot*dt)^2/(2*dx*(cNot*dt+dx))*([EzStore(end,2:end) 0] - 2*EzStore(end,:) + [0 EzStore(end,1:end-1)] +...
%                                            [EzStore(end-1,2:end) 0] - 2*EzStore(end-1,:) + [0 EzStore(end-1,1:end-1)]);
%                                        
% Ez(1,:) = - EzStoreOld(2,:) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(2,:)+EzStoreOld(1,:)) + 2*dx/(cNot*dt+dx)*(EzStore(1,:)+EzStore(2,:))+...
%           (cNot*dt)^2/(2*dx*(cNot*dt+dx))*([EzStore(1,2:end) 0] - 2*EzStore(1,:) + [0 EzStore(1,1:end-1)] +...
%                                            [EzStore(2,2:end) 0] - 2*EzStore(2,:) + [0 EzStore(2,1:end-1)]);
%                                        
% Ez(:,1) = - EzStoreOld(:,2) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(:,2)+EzStoreOld(:,1)) + 2*dx/(cNot*dt+dx)*(EzStore(:,1)+EzStore(:,2))+...
%           (cNot*dt)^2/(2*dx*(cNot*dt+dx))*(vertcat(EzStore(2:end,1), 0) - 2*EzStore(:,1) + vertcat(0, EzStore(1:end-1,1)) +...
%                                            vertcat(EzStore(2:end,2), 0) - 2*EzStore(:,2) + vertcat(0, EzStore(1:end-1,2)));                                       
% 
% Ez(:,end) = - EzStoreOld(:,end-1) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(:,end-1)+EzStoreOld(:,end)) + 2*dx/(cNot*dt+dx)*(EzStore(:,end)+EzStore(:,end-1))+...
%           (cNot*dt)^2/(2*dx*(cNot*dt+dx))*(vertcat(EzStore(2:end,end), 0) - 2*EzStore(:,end) + vertcat(0, EzStore(1:end-1,end)) +...
%                                            vertcat(EzStore(2:end,end-1), 0) - 2*EzStore(:,end-1) + vertcat(0, EzStore(1:end-1,end-1)));
%                                        
% % Corner Update                                       
% 
% Ez(1,1) = .5*( - EzStoreOld(1,2) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(1,2)+EzStoreOld(1,1)) + 2*dx/(cNot*dt+dx)*(EzStore(1,1)+EzStore(1,2)) + ...
%                - EzStoreOld(2,1) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(2,1)+EzStoreOld(1,1)) + 2*dx/(cNot*dt+dx)*(EzStore(1,1)+EzStore(2,1)));
% 
% 
% Ez(end,end) =  .5*(- EzStoreOld(end,end-1) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(end,end-1)+EzStoreOld(end,end)) + 2*dx/(cNot*dt+dx)*(EzStore(end,end)+EzStore(end,end-1))+ ...
%                    - EzStoreOld(end-1,end) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(end-1,end)+EzStoreOld(end,end)) + 2*dx/(cNot*dt+dx)*(EzStore(end,end)+EzStore(end-1,end)));
% 
% 
% Ez(1,end) = .5*(- EzStoreOld(1,end-1) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(1,end-1)+EzStoreOld(1,end)) + 2*dx/(cNot*dt+dx)*(EzStore(1,end)+EzStore(1,end-1))+ ...
%                 - EzStoreOld(2,end)    + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(2,end)+EzStoreOld(1,end))  + 2*dx/(cNot*dt+dx)*(EzStore(1,end)+EzStore(2,end)));
% 
% Ez(end,1) = .5*(- EzStoreOld(end,2) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(end,2)+EzStoreOld(end,1)) + 2*dx/(cNot*dt+dx)*(EzStore(end,1)+EzStore(end,2))+ ...
%                 - EzStoreOld(end-1,1) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(end-1,1)+EzStoreOld(end,1)) + 2*dx/(cNot*dt+dx)*(EzStore(end,1)+EzStore(end-1,1)));
%             

% Ez(:,1) = -EzStoreOld(:,2) + (cNot*dt-dx)/(cNot*dt+dx)*(Ez(:,2)+EzStoreOld(:,1)) + 2*dx/(cNot*dt+dx)*(EzStore(:,2)+EzStore(:,1));
%     
% Calculate Hy/Hx fields
Hy = Hy + Cn*(Ez(:,2:end) - Ez(:,1:end-1));
Hx = Hx - Cn*(Ez(2:end,:) - Ez(1:end-1,:));

EzStoreOld = EzStore;
EzStore = Ez;



%%
if mod(n,2) == 0;
% Generate plots/movie
grid = 2;
surf(xDom(1:grid:end, 1:grid:end), yDom(1:grid:end, 1:grid:end), ...
       Ez(1:grid:end, 1:grid:end))
caxis([-.05,.1])
axis([0 planeSize 0 planeSize 0 .5])
title('Gauss Pulse 2D Propagation')
xlabel('x')
ylabel('y')
zlabel('Ez')

timeStep = num2str(n);
time = strcat('Time Step: ', timeStep);
text(sourcePos(1), .5*sourcePos(2), .7, time);

% subplot(2,1,2)
% plot(xDom,Hy)
% axis([0 200 -1.2 1.2])
% xlabel('z')
% ylabel('Hy')
M(n/2) = getframe(gcf);
end

end
    

% movie2avi(M, 'Gauss Pulse 1D', 'fps', 40);

