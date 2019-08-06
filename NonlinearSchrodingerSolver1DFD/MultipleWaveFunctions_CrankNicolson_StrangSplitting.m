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
Npsi = 8;               % Number of seperate wavefunctions in domain
n = 200;                % Grid Points
tSteps = 2000;           % Time steps -- Sim 1: Evolution to ground state
tLength = 50;            % Max time   -- Sim 1: Evolution to ground state

tSteps2 = 500;          % Time steps -- Sim 2: Evolution in real time
tLength2 = 5000;          % Max time   -- Sim 2: Evolution in real time

xMin = -2*pi;           % x domain min
xMax = 2*pi;            % x domain max
xScale = 1;             % x scale

Nonlin = .0;             % Nonlinearity multiplier
f = 1/(pi);             % Oscillation frequency
omega = 2*pi*f;         % Angular frequency
dirac = 1;              % Normalized reduced plancks constant
me = 1;                 % Normalized e- mass
VxConst = .05;          % Multiplier for potential V(x)
dt = tLength/tSteps;    % Time steps

c1 = (1i*dirac)/(2*me); % 1. d^2/dx^2 Multiplier
c2 = -1i/dirac;         % 2. V(x)     Multiplier
c3 = -1i*Nonlin;        % 3. |Psi|^2  Multiplier

%% Vector Initilization
% Domain computations
dx = (xMax-xMin)/(n);
xDom = (-n/2:n/2-1)*dx;
% xDom = linspace(xMin,xMax,n);      
% dx = xDom(2)-xDom(1);
Vcoul = zeros(size(xDom));
% coul = -(1./abs(xDom));
% Lorentzian shape for potential
gamma = .1;
coul = 1./(pi*gamma*(1+(xDom/gamma).^2));



% Harmonic oscillator parabolic potential:
Vx = VxConst*me*(omega.^2)*xDom.^2/2;

% Initial condition -- select from a variety of initial conditions:
% Uinit = sech((xDom+1)/xScale)' - sech((xDom-1)/xScale)';
% Uinit = sech((xDom+1.5))'+ sech((xDom-1.5))' - sech((xDom))' ;    
% Uinit =(2*xDom.^2 - 1).*exp(-xDom.^2/2);
% Uinit =(8*xDom.^3 - 12*xDom).*exp(-xDom.^2/2);

Uinit1 = sech((xDom)/xScale)';   
Uinit2 = sech((xDom)/xScale)';   

% Normalize the initial solution Integral(|Psi|^2) = 1
Uinit1 = Uinit1/sqrt(trapz(xDom,abs(Uinit1).^2));
Uinit2 = Uinit2/sqrt(trapz(xDom,abs(Uinit2).^2));

UT = [Uinit1, Uinit2, Uinit1, Uinit2];
rho = zeros(length(UT),1);
R = zeros(size(UT));
Q = zeros(size(UT));

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
% Normalize the solution during evolution to ground state.
[Q R] = qr(UT);

UT = Q(:,1:Npsi);
 plot(xDom,UT(:,1), xDom, UT(:,2))
 legend('\Psi_1', '\Psi_2')

% Calculate electron density for current time step
for s = 1:Npsi
   rho = rho + abs(UT(:,s)).^2;
end
% Vcoul = fftshift( ifft( (fft([ rho(2:end); rho(1) ]) .* fftcoul) ) );
Vcoul = conv(rho,coul, 'same');
% plot(xDom, Vcoul, 'ro', xDom, Vcoul2, 'b*')



for s = 1:Npsi
    U = UT(:,s);  
    U = exp(dt*(c3.*Vcoul)/2).*U; 
    UT(:,s) = U;
end


 plot(xDom, Vcoul/160, xDom, rho, xDom, UT(:,1), xDom, UT(:,2))
 legend('V_{coulmb}', '\rho', '\Psi_1', '\Psi_2')



for r = 1:tSteps;
   
    for s = 1:Npsi
        
        % Calculate individual solutions based on e- density
        U = UT(:,s);    

        %%%%%%%%%%%%%%% TRIDIAGONAL SOLUTION TO LINEAR PART %%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%% END TRIDIAGONAL SOLUTION TO LINEAR PART %%%%%%%%%%%%%
        UT(:,s) = U;
    end
    
    % Calculate electron density for current time step
    for s = 1:Npsi
       rho = rho + abs(UT(:,s)).^2;
    end
%     Vcoul = fftshift( ifft( (fft([ rho(2:end); rho(1) ]) .* fftcoul) ) );
    Vcoul = conv(rho,coul, 'same');
    %%%%%%%%%%% STRANG SPLITTING TO COMPUTE NONLINEAR PART %%%%%%%%%%%%
    for s = 1:Npsi
        U = UT(:,s);  
        if r<tSteps
            U = exp(dt*(c3.*Vcoul)).*U;
        else
            U = exp(dt*(c3.*Vcoul)/2).*U; 
        end
         UT(:,s) = U;
    end
    %%%%%%%%% END STRANG SPLITTING TO COMPUTE NONLINEAR PART %%%%%%%%%%
    
    % Normalize the solution during evolution to ground state.
    [Q R] = qr(UT);

    UT = Q(:,1:Npsi);

    % Continuous plotting to display solution. Run time much faster if this
    % is removed.
    % plot(xDom, real(UT), 'b-', xDom, psi_0, 'b:', xDom, Vx, 'r--');
    if mod(r,1) == 0
    plot(xDom, real(UT), xDom, Vx, 'r--');
    text(-2, -.1, 'dt -> -i*dt ', 'fontsize', 18, 'fontname', 'cordia new')
    text(-2, -.5, strcat( 'time step: ', num2str(r)), 'fontsize', 18, 'fontname', 'cordia new')
    axis( [ xMin, xMax, -.15, .15] )
    getframe;
    end
    
end


    
return

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


