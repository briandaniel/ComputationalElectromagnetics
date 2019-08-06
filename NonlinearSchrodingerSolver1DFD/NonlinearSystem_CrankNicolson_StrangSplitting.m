%% Header
close all;
clear all;

%Physical Constants
epsNot = 8.854187817620e-12;        % Permittivity of Free Space
muNot  = 1.25663706e-6;             % Magnetic Constant
cNot   = 2.99792458e8;              % Speed of light in vacuum
etaNot = sqrt(muNot/epsNot);        % Characteristic impedence of free space
planck = 6.62606957e-34;            % Plancks constant
dirac  = planck/(2*pi);             % Dirac constant (Reduced Plancks constant)
e = 1.60217646e-19;                 % Charge of an electron (Coulombs)
me = 9.10938188e-31;                % Electron Mass (kg)

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
% 1. Second order (in both space and time) Crank-Nicholsen Finite         %   
%    Difference method is utilized to computer the linear solution.       %
% 2. The nonlinear part of the solution is combined using second order    %
%    Strang Splitting.                                                    %
% 3. The code is currently set to simulate the harmonic oscillator with   %
%    the parabolic trapping potential.                                    %
% 4. The code is run in imaginary time prior to running in real time in   %
%    order to produce the ground state                                    %
%                                                                         %
%                                                                         %
% Brian Hong, ACMS, 2012                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
n = 200;                % Grid Points
tSteps = 200;           % Time steps -- Sim 1: Evolution to ground state
tLength = 2;            % Max time   -- Sim 1: Evolution to ground state

tSteps2 = 500;          % Time steps -- Sim 2: Evolution in real time
tLength2 = 20;          % Max time   -- Sim 2: Evolution in real time

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
xDom = linspace(xMin,xMax,n);      
dx = xDom(2)-xDom(1);

% Harmonic oscillator parabolic potential:
Vx = me*(omega.^2)*xDom.^2/2;

% Initial condition -- select from a variety of initial conditions:
Uinit = sech(xDom/xScale)';    
% Uinit = sech((xDom+1)/xScale)' - sech((xDom-1)/xScale)';
% Uinit = sech((xDom+1.5))'+ sech((xDom-1.5))' - sech((xDom))' ;    
% Uinit =(2*xDom.^2 - 1).*exp(-xDom.^2/2);
% Uinit =(8*xDom.^3 - 12*xDom).*exp(-xDom.^2/2);


% Normalize the initial solution Integral(|Psi|^2) = 1
U = Uinit/sqrt(trapz(xDom,abs(Uinit).^2));
Usltn = zeros(tSteps,n);

% Generate exact linear solution for graph
E = dirac*omega;
xSQ = (xDom).^2;
A_o = 1;
psi_0 = A_o*exp(-sqrt(dirac*c2/-1i)*omega/2/dirac*xSQ);
psi_0 = psi_0 / sqrt(trapz(xDom, abs(psi_0).^2));

% Converting dt --> -i*dt
dt = -1i*dt;

% Generate (n+1) time step matrix
u1 = ( -dt*c1/(2*dx^2)*ones(1,n-1) )                ;
u2 = ( (1+dt*c1/dx^2)*ones(1,n) - (c2*dt/2)*Vx ) ;
u3 = ( -dt*c1/(2*dx^2)*ones(1,n-1) )                ;

% Generate (n) time step matrix
v1 = ( dt*c1/(2*dx^2)*ones(1,n-1) )                     ;
v2 = ( (1-dt*c1/(dx^2))*ones(1,n) + (c2*dt/2)*Vx ) ;
v3 = ( dt*c1/(2*dx^2)*ones(1,n-1) )                     ;
B = diag(v1,-1) + diag(v2) + diag(v3, 1);

%% 1st Sltn Loop
figure('outerposition', [200 400 1200 700])


U = exp(dt*(c3.*(abs(U)).^2)/2).*U;
for r = 1:tSteps;
    
    
    %%%%%%%%%%%%%%%%% TRIDIAGONAL SOLUTION TO LINEAR PART %%%%%%%%%%%%%%%%%
    % 1. Compute RHS
    U = B*U;
    a = u1; b = u2; c = u3;
    
    for k = 1:n-1;
        mu = a(k)/b(k);
        b(k+1) = b(k+1) - mu*c(k);
        U(k+1) = U(k+1) - mu*U(k);
    end
    
    U(n) = U(n)/b(n);
    
    for k = n-1:-1:1
        U(k) = ( U(k) - c(k)*U(k+1) )/(b(k));
    end
    %%%%%%%%%%%%%%% END TRIDIAGONAL SOLUTION TO LINEAR PART %%%%%%%%%%%%%%%
     
   
    %%%%%%%%%%%%% STRANG SPLITTING TO COMPUTE NONLINEAR PART %%%%%%%%%%%%%%
    if r<tSteps
        U = exp(dt*(c3.*(abs(U)).^2)).*U;
    else
        U = exp(dt*(c3.*(abs(U)).^2)/2).*U; 
    end
    %%%%%%%%%%% END STRANG SPLITTING TO COMPUTE NONLINEAR PART %%%%%%%%%%%%
    
    % Normalize the solution during evolution to ground state.
    U = U/sqrt(trapz(xDom,abs(U).^2));
    
    % Continuous plotting to display solution. Run time much faster if this
    % is removed.
    Usltn(r,:) = U;
    plot(xDom, real(U), 'b-', xDom, psi_0, 'b:', xDom, Vx, 'r--');
    text(-2, -.1, 'dt -> -i*dt ', 'fontsize', 18, 'fontname', 'cordia new')
    text(-2, -.5, strcat( 'time step: ', num2str(r)), 'fontsize', 18, 'fontname', 'cordia new')
    axis( [ xMin, xMax, -1, 1] )
    getframe;

end


%% Second loop to run the solution in real time.
tSteps = tSteps2;
tLength = tLength2;
dt = tLength/tSteps;

% Generate (n+1) time step matrix
u1 = ( -dt*c1/(2*dx^2)*ones(1,n-1) )                ;
u2 = ( (1+dt*c1/dx^2)*ones(1,n) - (c2*dt/2)*Vx ) ;
u3 = ( -dt*c1/(2*dx^2)*ones(1,n-1) )                ;

% Generate (n) time step matrix
v1 = ( dt*c1/(2*dx^2)*ones(1,n-1) )                     ;
v2 = ( (1-dt*c1/(dx^2))*ones(1,n) + (c2*dt/2)*Vx ) ;
v3 = ( dt*c1/(2*dx^2)*ones(1,n-1) )                     ;
B = diag(v1,-1) + diag(v2) + diag(v3, 1);


U = exp(dt*(c3.*(abs(U)).^2)/2).*U;
for r = 1:tSteps;
    
    %%%%%%%%%%%%%%%%% TRIDIAGONAL SOLUTION TO LINEAR PART %%%%%%%%%%%%%%%%%
    % 1. Compute RHS
    U = B*U;
    a = u1; b = u2; c = u3;
    
    for k = 1:n-1;
        mu = a(k)/b(k);
        b(k+1) = b(k+1) - mu*c(k);
        U(k+1) = U(k+1) - mu*U(k);
    end
    
    U(n) = U(n)/b(n);
    
    for k = n-1:-1:1
        U(k) = ( U(k) - c(k)*U(k+1) )/(b(k));
    end
    %%%%%%%%%%%%%%% END TRIDIAGONAL SOLUTION TO LINEAR PART %%%%%%%%%%%%%%%
     
   
    %%%%%%%%%%%%% STRANG SPLITTING TO COMPUTE NONLINEAR PART %%%%%%%%%%%%%%
    if r<tSteps
        U = exp(dt*(c3.*(abs(U)).^2)).*U;
    else
        U = exp(dt*(c3.*(abs(U)).^2)/2).*U; 
    end
    %%%%%%%%%%% END STRANG SPLITTING TO COMPUTE NONLINEAR PART %%%%%%%%%%%%
    
    % Continuous plotting to display solution. Run time much faster if this
    % is removed.
    Usltn(r,:) = U;
    plot(xDom, real(U), 'b-', xDom, psi_0, 'b:', xDom, Vx, 'r--');
    text(-2, -.1, 'dt -> -i*dt ', 'fontsize', 18, 'fontname', 'cordia new')
    text(-2, -.5, strcat( 'time step: ', num2str(r)), 'fontsize', 18, 'fontname', 'cordia new')
    axis( [ xMin, xMax, -1, 1] )
    getframe;
    
end


