%% Header
close all;
clear all;

%% Information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This program computes solutions to the time dependent nonlinear         %
% schrodinger equation of the form                                        %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%           dPsi      d^2Psi                                              %
%           ---- = c1*------ + c2*V(x)*Psi + c3*|Psi|^2*Psi               %
%            dt        dx^2                                               %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Properties:                                                             %
% 1. Fourth order (in time) using the fourth order Runge-Kutta time       %  
%    stepping.                                                            %
% 2. The linear part is computed using a fft and if periodic boundary     %
%    conditions are satisfied the spatial dimension has spectral          %
%    convergence                                                          %
% 3. The nonlinear part is computed exactly at each time step, which is   %
%    used by the Runge-Kutta solver.                                      %
% 4. The code accesses an external function 'nls.m' to perform the        %
%    calculations for the PDE specfied above.                             %
% 5. The code is currently set to simulate the harmonic oscillator with   %
%    the parabolic trapping potential.                                    %
% 6. The code is run in imaginary time prior to running in real time in   %
%    order to produce the ground state                                    %
%                                                                         %
%                                                                         %
% Brian Hong, ACMS, 2012                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
N = 200;                % Grid Points
tSteps = 1e3;           % Time steps -- Sim 1: Evolution to ground state
tLength = 1.2;          % Max time   -- Sim 1: Evolution to ground state

tSteps2 = 6e3;          % Time steps -- Sim 2: Evolution in real time
tLength2 = 10;          % Max time   -- Sim 2: Evolution in real time

xMin = -2*pi;           % x domain min
xMax = 2*pi;            % x domain max
xScale = 1;             % x scale

Nonlin = 0;             % Nonlinearity multiplier
f = 1/(pi);             % Oscillation frequency
omega = 2*pi*f;         % Angular frequency
dirac = 1;              % Normalized reduced plancks constant
me = 1;                 % Normalized e- mass
dt = tLength/tSteps;    % Time steps

c1 = (1i*dirac)/(2*me); % 1. d^2/dx^2 Multiplier
c2 = -1i/dirac;         % 2. V(x)     Multiplier
c3 = -1i*Nonlin;        % 3. |Psi|^2  Multiplier

%% Vector Initilization
% Domain computations
L = xMax-xMin;
dx = (xMax-xMin)/(N);
xDom = dx*(-N/2:N/2-1);
k = (2*pi/L)*[0:N/2-1, -N/2:-1]';

% Initial condition -- select from a variety of initial conditions:
Uinit = sech(xDom/xScale)';    
% Uinit = sech((xDom+1)/xScale)' - sech((xDom-1)/xScale)';
% Uinit = sech((xDom+1.5))'+ sech((xDom-1.5))' - sech((xDom))' ;    
% Uinit =(2*xDom.^2 - 1).*exp(-xDom.^2/2);
% Uinit =(8*xDom.^3 - 12*xDom).*exp(-xDom.^2/2);

% Normalize the initial solution Integral(|Psi|^2) = 1
U = Uinit/sqrt(trapz(xDom,abs(Uinit).^2));

% Harmonic oscillator parabolic potential:
Vx = me*(omega.^2)*xDom.^2'/2;

% Call nls function as specified in 'nls.m'
nls_solver = @(t,x) nls(t,x, k, Vx, c1, c2, c3);

% Converting dt --> -i*dt
dt = -1i*dt;
t = 0;

figure('outerposition', [200 400 1200 700])
%% 1st Sltn Loop
for k = 1:tSteps;   
    
    % 4th order Runge-Kutta Solver -- calls 'nls.m' to calculate values at 
    % various time steps
    k1 = dt * nls_solver(t, U);
	k2 = dt * nls_solver(t + dt/2, U + k1/2 );
	k3 = dt * nls_solver(t + dt/2, U + k2/2 );
	k4 = dt * nls_solver(t + dt, U + k3 );
    U = U + ( k1 + 2*k2 + 2*k3 + k4 ) / 6;
    U = U/sqrt(trapz(xDom,abs(U).^2));
    t = t + dt;
    
    if mod(k,30) == 0
        plot(xDom, real(U))
        axis([xMin xMax -1 1])
        text(-4, -.1, 'dt -> -i*dt ', 'fontsize', 18, 'fontname', 'cordia new')
        text(-4, -.5, strcat( 't =  ', num2str(t)), 'fontsize', 18, 'fontname', 'cordia new')
        getframe;
    end
    
    
end

%% Second loop to run the solution in real time.
tSteps = tSteps2;
tLength = tLength2;
dt = tLength/tSteps;
t = 0;

for k = 1:tSteps;   
     
    % 4th order Runge-Kutta Solver -- calls 'nls.m' to calculate values at 
    % various time steps
    k1 = dt * nls_solver(t, U);
	k2 = dt * nls_solver(t + dt/2, U + k1/2 );
	k3 = dt * nls_solver(t + dt/2, U + k2/2 );
	k4 = dt * nls_solver(t + dt, U + k3 );
    U = U + ( k1 + 2*k2 + 2*k3 + k4 ) / 6;
    U = U/sqrt(trapz(xDom,abs(U).^2));
    t = t + dt;
    
    if mod(k,30) == 0
        plot(xDom, real(U))
        axis([xMin xMax -1 1])
        text(-4, -.1, 'dt -> dt ', 'fontsize', 18, 'fontname', 'cordia new')
        text(-4, -.5, strcat( 't =  ', num2str(t)), 'fontsize', 18, 'fontname', 'cordia new')
        getframe;
    end
    
    
end
