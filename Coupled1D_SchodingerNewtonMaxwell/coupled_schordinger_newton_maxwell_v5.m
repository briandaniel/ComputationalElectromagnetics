%% Header
clear all
close all

%Physical Constants
epsNot = 8.854187817620e-12;        % Permittivity of Free Space
muNot  = 1.25663706e-6;             % Magnetic Constant
cNot   = 2.99792458e8;              % Speed of light in vacuum
etaNot = sqrt(muNot/epsNot);        % Characteristic impedence of free space
planck = 6.62606957e-34;            % Plancks constant
hbar  = planck/(2*pi);             % Dirac constant (Reduced Plancks constant)
e = 1.60217646e-19;                 % Charge of an electron (Coulombs)
me = 9.10938188e-31;                % Electron Mass (kg)
jeV = 6.24150974e18;


%% Information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This program computes solutions to three models coupled into a single   %
% solver in a 1D toy setting. The models are:                             %
%                                                                         %
%  I. Classical Maxwells' Equations for light propagation                 %
%  II. Schrodinger Equation for quantum electron evolution                %
%  III. Classical Newton's equation for ion motion                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  I. FDTD solver for maxwells equations                                  %
%-------------------------------------------------------------------------%
%                                                                         %
%      dHy        dEx                                                     % 
%  1. ----- = a1*-----                                                    % 
%      dt         dz                                                      % 
%                                                                         %
%      dEx        dHy                                                     % 
%  2. ----- = a2*----- + a3*Je + a4*Jion                                  % 
%      dt         dz                                                      % 
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%      where                                                              %
%                          dPsi_j      dPsi*_j                            % 
%      Je =  Sum [ Psi*_j -------- +  --------- Psi_j ]                   % 
%                            dz          dz                               % 
%                                                                         %
%      Jion = Sum [ Zn * Vn * delta( z - Rn ) ]                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  II. Crank Nicholson scheme with Strang-Splitting to solve the          %
%      Nonlinear Schrodinger equation for quantum electrons               %
%-------------------------------------------------------------------------%
%                                                                         %  
%  dPsi      d^2Psi                                                       % 
%  ---- = b1*------ + b2*(z*Ex) + b3*F(rho(z)) + b4*G(Rn)                 % 
%   dt        dz^2                                                        % 
%                                                                         % 
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%      where                                                              %
%      rho(x) = Sum |Psi_i(x)|^2                                          %
%                                                                         %
%                    rho( z')                                             %
%      F(rho) = Int ---------- dz'                                        %
%                   | z - z' |                                            %
%                i.e. F is the convolution of rho(x) and u(x)             %
%                                                                         %
%                      Zn                                                 %
%      G(Rn) = Sum ---------- dz'                                         %
%                  | z - Rn |                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  III. Finite difference solution to Newton's equation                   %
%-------------------------------------------------------------------------%
%   The second order differential equation is split into a system of two  %
% first order differential equations.                                     %
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%                                                                         %
%          dRn                                                            %
%     1.  ----- = Vn                                                   %
%           dt                                                            %
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%     where                                                               %
%     Rn = position of nth ion (Rion)                                     %
%     Vn = velocity of nth ion (Vion)                                     %
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%                                                                         %
%          dVn                                                            %
%     2.  ----- = c1*F_Ion + c2*F_e + c3*F_Lorentz                        %
%          dt                                                             %
%                                                                         % 
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%     where                                                               %
%                            Rn - Rm                                      %
%     F_Ion = Zn*Zm* Sum  --------------                                  % 
%                    m~=n | Rn  - Rm |^3                                  %   
%                                                                         %
%                  d         rho(z)                                       %
%     F_Ion = Zn* ---- Int ---------- dz                                  % 
%                  dR      | z - Rn |                                     %   
%                                                                         %
%     F_Lorentz = Zn* ( Ex + Vn x By)                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
% Brian Hong, ACMS, January 2013                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DOMAIN SETTINGS %%%%%%%%%%1%%%%%%%%%%%%%%%%%%%
Npsi = 2;                   % Number of electron wavefunctions in domain
Nion = 2;                   % Number of ions
N = 401;                    % Grid Points
EMpoints = 15;              % Additional points in EM field domain
NHy = N + EMpoints*2;       % H-field grid points 
NEx = N + EMpoints*2 + 1;   % H-field grid points


mult = 1;
tSteps = 1000*mult;         % Time steps -- Sim 1: Evolution to ground state
tLength = 20e-16*mult;       % [s] Max time -- Sim 1: Evolution to ground state

zMin = 0e-9;                % [m] z domain min 
zMax = 6e-9;                % [m] z domain max 
surfacePos = 2e-9;         % [m] position of surface


dz = (zMax-zMin)/(N-1);     % z step
dt = tLength/tSteps;        % time step
%=========================================================================%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mz = me*2.4e4;                % Ion mass
Zn = 1;                     % Ion charge
guess = .5;                 % Guess initial spacing for optimization
spin = 1;
                            %    / Interaction Normalization 
gamma =4e-11;              %<--| (gamma --> 0 implies 
                            %    \ 1/sqrt(gamma^2+z^2) --> 1/|z|.
                            
epsFactor = 1/(4*pi*epsNot);
%=========================================================================%


%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDITIONAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%
EMpml = EMpoints - 5;
PSIpml = 10;
%=========================================================================%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E.M. Field (Ex, Hy, Je, Jion)

slowLight = 1/10;
epsNot = 1/slowLight* 8.854187817620e-12;        % Permittivity of Free Space
muNot  = 1/slowLight* 1.25663706e-6;             % Magnetic Constant
cNot   = 1/sqrt(muNot*epsNot);      % Speed of light in vacuum

muy  = muNot*ones(N+2*EMpoints,1)  ;
epsx = epsNot*ones(N+2*EMpoints+1,1);

a1 = 1 ./( muy );                           % Coef. for dEx/dz in Hy update
a2 = 1 ./( epsx );                          % Coef. for dHy/dz in Ex update
a3 = ( 1i*hbar*e )/( 2*me );                % Coef. for electron current
a4 = -e;                                    % Coef. for ion current

% Electrons (Psi, rho)
b1 = ( 1i * hbar )/( 2 * me );              % Coef. for dPsi^2/dz^2 term
b2 = ( -1i .* e )/( hbar );                 % Coef. for E-field contribution
b3 = epsFactor * ( -1i .* e^2 )/( hbar );   % Coef. for electron contribution
b4 = epsFactor * ( -1i .* e^2 )/( hbar );   % Coef. for ion contribution

% Ion (Vion, Rion)
c1 = epsFactor* e^2;                        % Coef. of ionic force
c2 = epsFactor* e^2;                        % Coef. of electronic force
c3 = e;                                     % Coef. of Lorentz force
%=========================================================================%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PML Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigM         = 3;               % Power of sigma in PML                    [USER DEFINED]
sigMax       = 3*(sigM+1)*... % Max value of sigma in PML
              .8/(etaNot*dz);
kappaMax     = 1;               % Max value of kappa in PML                [USER DEFINED]
kappaM       = 3;               % Power of kappa in PML                    [USER DEFINED]
aMax         = 1;               % Max value of kappa in PML                [USER DEFINED]
aM           = 3;               % Power of a in PML                        [USER DEFINED]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Domain computations & Vector Initilization
% Domain Computations
L = zMax-zMin;
zDomHalf = dz*( 0:N-1 );
zGraph = ( zDomHalf - zMin )*1e9;
zD = zDomHalf;
newPts = (NHy-N)/2;
zDomH = [dz*(-newPts:-1) , zDomHalf , ( zDomHalf(end) + (1:newPts)*dz )];
zDom = [ zDomH - .5*dz, zDomH(end)+.5*dz];
% kz = (2*pi/L)*[0:N/2-1, -N/2:-1]';

% Electron vectors
PSI = zeros(N,Npsi);
uZ = zeros(N,Npsi);

% Ion vectors
Vion = zeros(Nion,1);
Rion = zeros(Nion,1);
Zion = zeros(Nion,1);
Mion = zeros(Nion,1);

% Electromagnetic Radiation vectors
Ex = zeros(N+2*EMpoints+1,1);
Hy = zeros(N+2*EMpoints,1);


%% Pre-loop computations
% Initial condition
Uinit1 = sech((zDomHalf))';   
 
% Normalize the initial solution Integral(|Psi|^2) = 1
Uinit1 = Uinit1/sqrt(trapz(zDomHalf,abs(Uinit1).^2));
PSI(:,1) = Uinit1;

% Calculate interaction terms for use in fft convolution
mid = (zDomHalf(end) - zDomHalf(1))/2 + zDomHalf(1);
coul = 1./sqrt(gamma^2 + ( zDomHalf - mid ).^2); 
fftCoul = fft([coul'; ones(N-1,1)*coul(end)]);

% Normalize the solution before evolution to ground state
[Q R] = qr(PSI);
PSI = Q(:,1:Npsi);

% Fill terms with mass and charge of ions
for s = 1:Nion
    Zion(s) = Zn;
    Mion(s) = mz;
end


%% PML Calculations
% Compute PML Coefficients
zVec = (1:EMpml)-.5;
zVecShift = zVec+.5;

% Basic Vectors
sigVec = (abs(zVec+.5).^sigM./EMpml.^sigM).*sigMax;
sigVecStag = (abs(zVec).^sigM./EMpml.^sigM).*sigMax;

kappaVec = 1+(abs(zVec+.5).^kappaM./EMpml.^kappaM).*(kappaMax-1);
kappaVecStag = 1+(abs(zVec).^kappaM./EMpml.^kappaM).*(kappaMax-1);

aVec = (abs(zVec+.5).^aM./EMpml.^aM).*aMax;
aVecStag = (abs(zVec).^aM./EMpml.^aM).*aMax;

% Full vector set spanning grid space
sigVecD = [sigVec(end:-1:1) zeros(1,NEx-2*length(sigVec)) sigVec];
sigVecStagD = [sigVecStag(end:-1:1) zeros(1,NEx-2*length(sigVec)-1) sigVecStag];

kappaVecD = [kappaVec(end:-1:1) ones(1,NEx-2*length(kappaVec)) kappaVec];
kappaVecStagD = [kappaVecStag(end:-1:1) ones(1,NEx-2*length(kappaVec)-1) kappaVecStag];

aVecD = [aVec(end:-1:1) zeros(1,NEx-2*length(aVec)) aVec];
aVecStagD = [aVecStag(end:-1:1) zeros(1,NEx-2*length(aVec)-1) aVecStag];

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

bx = bx';
by = by';
gx = gx';
gy = gy';
kx = kx';
ky = ky';


%% Optimize ground state distance
figure('outerposition', [100 400 1000 740])

% Function for finding the energy of various positions
energyFun = @(x) systemEnergy(x, N, Npsi, Nion, zDomHalf, dz, dt, ...
                              tSteps, PSI, spin,  Zion,  b1, b3, b4, ...
                              gamma, e, hbar, me, epsFactor, fftCoul);    

% Options for minimizer                          
options = optimset ( 'tolx', 1e-4, 'maxfunevals', 30);

% Find minimum energy ion spacing, and put into "Rion"
% minParam = fminsearch( energyFun, guess, options);

minParam = guess;
mid = (zDomHalf(end) - zDomHalf(1))/2 + zDomHalf(1);
count = linspace(-1,1,Nion);
Rground = mid + minParam*(zMax - zMin)./3*count; 

%%
Rion = Rground - Rground(1) + surfacePos;
groundSteps = tSteps*2;

% Evaluate final ground state for optimal ion distance
[ PSIground, Etot, Ue, Te, Ubond, Uion] = ...
    ground( Rion, N, Npsi, Nion, zDomHalf, dz, dt, groundSteps, PSI, ... 
            spin,  Zion,  b1, b3, b4, gamma, e, hbar, ...
            me, epsFactor, fftCoul);
Vground = Vion;
Rground = Rion;
% Display energy results
display('')
display(strcat('Optimal distance:  :', num2str( (Rground(2) - Rground(1)) *1e9), 'nm'))
display(strcat('ke:  :', num2str(Te*jeV)))
display(strcat('bond energy:  :', num2str(Ubond*jeV)))
display(strcat('e- potential:  :', num2str(Ue*jeV)))
display(strcat('ion potential:  :', num2str(Uion*jeV)))
display(strcat('Total energy:  :', num2str( Etot*jeV )))
display('')


%% 2nd Solution loop (real time)
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DOMAIN SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tSteps = 20000;          % Time steps -- Sim 2:
tLength = 8e-15;      % Max time   -- Sim 2: 
dt = tLength/tSteps;        % time step

normDt = dt;
dt = .5*dz/(cNot);  

dt/normDt
%=========================================================================%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on = 1, off = 0
emSwitch = 1;               % Turn EM field on or off
motionSwitch = 1;           % Turn ion motion on or off
psiSwitch = 1;              % Turn electrons on or off 
%=========================================================================%

%%%%%%%%%%%%%%%%%%%%%%%% Hard source parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian = 1,  Modulated Gaussian = 2,  Sine Wave = 3, 
% Ricker Wavelet = 4;
pulseSetting = 2;               % Choose which source condition to use     [USER DEFINED]
sourcePos    = EMpml+1;
% Gaussian Pulse parameters
width        = 2e-15;          % Approximate width of pulse in time (sec) [USER DEFINED]
maxT         = width*.5;
% Sine Wave Parameters
f            = 5e15;              % Frequency                                [USER DEFINED]
amp          = 2e11;               % Amplitude of Wave                        [USER DEFINED]
lambda       = cNot/f;          % Wavelength
omega        = 2*pi*f;          % Angular Frequency
% Ricker Wavelet Parameters
fp           = 5e8;             % Peak Frequency                           [USER DEFINED]
md           = 1;               % Temporal Delay Multiple                  [USER DEFINED]
dr           = md/fp;           % Temporal Delay
%=========================================================================%

%%%%%%%%%%%%%%%%%%%%%%% Frequency domain parameters %%%%%%%%%%%%%%%%%%%%%%%
freqWindow = [ .001*f, 3.5*f];
stepsF       = 200;
detInitial   = sourcePos+1 ;        % Position of initial detector             [USER DEFINED]
detPost      = NEx - 20;            % Position of latter detector              [USER DEFINED]
detRef       = sourcePos+1;         % Position of reflection detector
%=========================================================================%


%%%%%%%%%%%%%%%%%%%%%%% Frequency domain parameters %%%%%%%%%%%%%%%%%%%%%%%
% Frequency domain calculations
freqMin = freqWindow(1);
freqMax = freqWindow(2);
testFreq = linspace(freqMin, freqMax, stepsF);
testWL = cNot./testFreq;
EfInit  = zeros(size(testFreq));
EfRef  = zeros(size(testFreq));
EfFinal = zeros(size(testFreq));
changeTime = width*1.5/dt;

PSI = PSIground;
Vion = Vground;
Rion = Rground';

rho = zeros(length(PSI),1);
for s = 1:Npsi 
    rho = rho + (abs(PSI(:,s)).^2);
end
maxPSI = max(max(abs(real(PSI))));
maxRho = max(max(abs(real(rho))));

% Electromagnetic Radiation vectors
Ex = zeros(N+2*EMpoints+1,1);
Hy = zeros(N+2*EMpoints,1);
By = zeros(size(Hy));
HyOld = zeros(size(Hy));
ByOld = zeros(size(Hy));
Jx = zeros(size(Ex));
muy = ones(size(Hy));
epsx = ones(size(Ex));
source = 0;

% PML Auxilliary variables
QEx = zeros(size(Hy));
QHy = zeros(size(Ex));

% First derivative of U using fourth order F.D. approximation.
s1 = 8*ones(1,N-1);
s2 = ones(1,N-2);
D1 = diag(s1,1) - diag(s1, -1) - diag(s2,2) + diag(s2,-2);
% Applying the boundary conditions
D1(end,1) = 8;
D1(end-1,1) = -1;
D1(end,2) = -1;
D1(1,end) = -8;
D1(1,end-1) = 1;
D1(2,end) = 1;
% Divide each element by constant 12*dx
D1 = D1./(12*dz);


% Generate (n+1) time step matrix
u1 = ( -dt*b1/(2*dz^2)*ones(1,N-1) ) ;
u2 = ( (1+dt*b1/dz^2)*ones(1,N) )    ;
u3 = ( -dt*b1/(2*dz^2)*ones(1,N-1) ) ;

% Generate (n) time step matrix
v1 = ( dt*b1/(2*dz^2)*ones(1,N-1) )  ;
v2 = ( (1-dt*b1/(dz^2))*ones(1,N) )  ;
v3 = ( dt*b1/(2*dz^2)*ones(1,N-1) )  ;

B = diag(v1,-1) + diag(v2) + diag(v3, 1);

stationary = Rion(1);
psiTime = 0;
ionTime = 0;
fieldTime = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLUTION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('outerposition', [200 50 1500 1140])
for r = 1:tSteps
  
    tic
    % step 1 of 5 %
    %%%%%%%%%%%%%%%%%%%%%% ION POSITION (R_n) UPDATE %%%%%%%%%%%%%%%%%%%%%%
    Rold = Rion;
    if motionSwitch == 1
        Rion = Rion + dt*Vion;
    end
    
    Ravg = ( Rold + Rion )/2;
    %=====================================================================%
    ionTime = ionTime + toc;
    
    tic
    if emSwitch == 1;
    % step 2 of 4 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ex UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Source Condition  
    if pulseSetting == 1            % Gaussian Source
        source = exp(-(maxT-r*dt)^2/(width/5)^2);
    elseif pulseSetting == 2        % Modulated Gaussian Source
        source = exp(-(maxT-r*dt)^2/(width/5)^2) * amp*sin(omega*r*dt);
    elseif pulseSetting == 3        % CW source
        source = amp*sin(omega*r*dt);
    elseif pulseSetting == 4        % Ricker Wavelet source
        source = (1-2*(pi*fp*(r*dt-dr))^2)*exp(-(pi*fp*(r*dt-dr))^2);
    end
    Ex(sourcePos) = source/2 +  Ex(sourcePos);

    % Current vectors
    Jx = zeros(NEx,1);
    Jpsi = zeros(N-1,1);
    Jion = zeros(NEx,1);

    % Current contribution from electrons (Jpsi)
    for s = 1:Npsi;
        Jpsi = Jpsi + ...
              ( ( conj(PSI(2:end,s)) + conj(PSI(1:end-1,s)) )/2 .*     ...
                (       PSI(2:end,s) - PSI(1:end-1,s)       )/dz   ) - ...
              ( (       PSI(2:end,s) + PSI(1:end-1,s)       )/2 .*     ...
                ( conj(PSI(2:end,s)) - conj(PSI(1:end-1,s)) )/dz   );
    end

    % Current contribution from ions (Jion)
    for s = 1:Nion
        [C I1] = min(abs(zDom - Ravg(s)));
        if ( zDom(I1) - Ravg(s) < 0)
            I2 = I1+1;
        else
            I2 = I1;
            I1 = I1-1;
        end

        rightWeight = (Ravg(s) - zDom(I1))/dz;
        leftWeight = (zDom(I2) - Ravg(s))/dz;

        Jion([I1, I2]) = Jion([I1, I2]) + ...
                         Zion(s)*Vion(s)*[leftWeight ; rightWeight];
    end

    % Total current (Jx)
    Jx(EMpoints + 2: NEx-EMpoints-1) = a3*Jpsi;
    Jx = Jx + a4*Jion;

    % E-field Update
    QHy(2:end-1) = bx(2:end-1).*QHy(2:end-1) ... 
                  - gx(2:end-1).*(Hy(2:end) - Hy(1:end-1))/dz;

    Ex(2:end-1) = Ex(2:end-1) + a2(2:end-1).* dt .* ...
                 ( ( Hy(2:NEx-1) - Hy(1:NEx-2) )/dz + QHy(2:end-1)) ...
                 - dt.*Jx(2:end-1) ;   

    % Frequency domain update
    EfFinal = EfFinal + Ex(detPost).*exp(-1i.*2.*pi.*dt*r.*testFreq);
    if r < changeTime
        EfInit = EfInit + Ex(detInitial).*exp(-1i.*2.*pi.*dt*r.*testFreq);           
    elseif r > changeTime
        EfRef = EfRef + Ex(detInitial).*exp(-1i.*2.*pi.*dt*r.*testFreq);
    end
    
    % Hy UPDATE 
    Hyold = Hy;

    QEx = by.*QEx - gy.*(Ex(2:end) - Ex(1:end-1))/dz;
    Hy = Hy + a1.*dt.*( (Ex(2:end) - Ex(1:end-1))/dz + QEx );
    %=====================================================================%
    end
    fieldTime  = fieldTime + toc;

    tic
    % step 3 of 4 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PSI UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Strang split requires 3 steps, i, ii, iii. 
    if psiSwitch == 1

    uZ = zeros(1,N);
    for s = 1:Nion
    uZ = uZ - Zion(s)...
         *abs( 1./sqrt(gamma^2 + ( zDomHalf - Rion(s) ).^2) );      
    end
    uZ = uZ';


    % Electric field potential contribution 
    uE = ( b2* ( dt/4 ) ).*...
         ( zDom(EMpoints + 1 : NEx-EMpoints-1).*Ex(EMpoints + 1 : NEx-EMpoints-1)'...
         + zDom(EMpoints + 2 : NEx-EMpoints  ).*Ex(EMpoints + 2 : NEx-EMpoints  )'    );
    uE = uE';

    % Calculate electron density for current time step
    rho = zeros(length(PSI),1);
    for s = 1:Npsi
       rho = rho + abs(PSI(:,s)).^2;
    end

    fftRho = fft( [rho; zeros(N-1,1)]);      
    Vcoul2  =  ifft(fftRho.* fftCoul); 
    Vcoul  = Vcoul2(floor(N/2+1):floor(3*N/2));

    %%%%%%%%%%% STRANG SPLITTING TO COMPPSIE NONLINEAR PART %%%%%%%%%%%%
    for s = 1:Npsi
        U = PSI(:,s);  
        U = exp(dt*(b3.*Vcoul)/2 + b4*dt*uZ/2 + uE/2).*U; 
        PSI(:,s) = U;
    end
    %%%%%%%%% END STRANG SPLITTING TO COMPPSIE NONLINEAR PART %%%%%%%%%%

    for s = 1:Npsi
        PSI(:,s) = B*PSI(:,s);
        a = u1; b = u2; c = u3;
        PSI(:,s) = tridisolve(a, b, c, PSI(:,s));
    end

    % Calculate electron density for current time step
    rho = zeros(length(PSI),1);
    for s = 1:Npsi
       rho = rho + abs(PSI(:,s)).^2;
    end

    fftRho = fft( [rho; zeros(N-1,1)]);      
    Vcoul2  =  ifft(fftRho.* fftCoul); 
    Vcoul  = Vcoul2(floor(N/2+1):floor(3*N/2));


    %%%%%%%%%%% STRANG SPLITTING TO COMPPSIE NONLINEAR PART %%%%%%%%%%%%
    for s = 1:Npsi
        U = PSI(:,s);  
        U = exp(dt*(b3.*Vcoul)/2 + b4*dt*uZ/2 + uE/2).*U; 
        PSI(:,s) = U;
    end
    %%%%%%%%% END STRANG SPLITTING TO COMPPSIE NONLINEAR PART %%%%%%%%%%
    end
    %=====================================================================%
    psiTime = toc + psiTime;
    
    tic
    % step 5 of 5 %
    %%%%%%%%%%%%%%%%%%%%%%%%% ION VELOCITY UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%
    if motionSwitch == 1
% %         Calculate Bfield/Efield at ion locations
%         By = muy.*Hy;
%         ByOld = muy.*HyOld;
%         ByR = spline(zDomH, By, Rion);
%         ByOldR = spline(zDomH, ByOld, Rion);
%         ExR = spline(zDom, Ex, Rion);

%         Calculate inverse Bfield Multiplier
%         Binv = 1./(1 - ( (e*Zion*dt)./(2*Mion) ).*ByR);
        Binv = 1;
        % Calculate ion-ion interaction (Coulomb Force)
        Fion = zeros(Nion,1);
        for s = 1:Nion   
            for ss = 1:Nion
                if s ~= ss;
                    zPos = Rion(s) - Rion(ss);
                    const = Zion(s)*Zion(ss)*e^2.*epsFactor/Mion(s);
                    Fion(s) = const* ( zPos / ( gamma^2 + zPos^2 )^(3/2) );
                end
            end
        end

        % Calculate electron-ion interaction (Coulomb Force)
        rho = zeros(length(PSI),1);
        for s = 1:Npsi 
            rho = rho + (abs(PSI(:,s)).^2);
        end

        fftRho = fft([rho; zeros(N-1,1)]);        
        convRho = (ifft(fftRho.*fftCoul));
        convRho = convRho(floor(N/2+1):floor(3*N/2));

%         d1rho = splineDerivative( zDomHalf', convRho);
        d1rho  = D1* convRho;  
%         d1rhoIon = spline(zDomHalf, d1rho, Rion);
        d1rhoIon = interp1(zDomHalf, d1rho, Rion);
        Fpsi = ( (e^2*Zion*epsFactor)./Mion ) .*d1rhoIon;

        FpsiT = c2*dt*Fpsi;
        FionT = c1*dt*Fion;

        % Calculate E-field interaction (Lorentz Force)
%         Florz = ( (e*Zion.*dt)./Mion ).*(ExR + Vion.*ByOldR/2);
        Florz = 0;
        Vion = Binv.*(Vion + FionT + FpsiT + c3*Florz);
        
    end
    %=====================================================================%
    ionTime = toc + ionTime;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% DATA VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Continuous plotting to display solution.
    if mod(r,100) == 0;
        hold on
        
       colors = hsv(round(3*Npsi));
        for s = 1:Npsi
%             plot( zDomHalf, .9*real(PSI(:,s))/(maxPSI) , 'color', colors(s+2*Npsi,:))
        end
        
        plot(zDomHalf, -.9*abs(uZ)/max(abs(uZ)), 'k--')
        plot(zDomHalf, .4*real(rho)/maxRho, 'r.')
        plot(zDom, Ex/amp, ':', 'linewidth', 2) 
        plot(zDom, .4*Jx/(max(abs(Jx))) )
        plot(Rion, zeros(size(Rion)), 'bo', 'markersize', 12, 'linewidth', 1 );
        plot(Rion, zeros(size(Rion)), 'b*', 'markersize', 12, 'linewidth', 1 );
        for s = 1:Nion-1
            
            plot( [Rion(s+1) Rion(s)] , [ -.02 -.02 ], 'k-')
            text( .97*( Rion(s+1) + Rion(s) )/2, -.05 , strcat( num2str( ...
                  round ((Rion(s+1) - Rion(s) ).*1e14 )/1e5 ) , ' [nm]'));
        end
        %% Calculate energy of system

        % Second derivative of U using fourth order F.D. approximation.
        %
        %          -U(x+2h) + 16U(x+h) - 30U(x) + 16U(x-h) - U(x-2h)
        % U_xx =  ---------------------------------------------------  + O(h^4)
        %                              12h^2

        % The corresponding matrix, D2 (Second derivative matrix)
        e1 = -30*ones(1,N);
        e2 = 16*ones(1,N-1);
        e3 = -1*ones(1,N-2);
        D2 = diag(e1) + diag(e2,-1) + diag(e2, 1) + diag(e3,-2) + diag(e3,2);
        D2(end,1) = 16;
        D2(end-1,1) = -1;
        D2(end,2) = -1;
        D2(1,end) = 16;
        D2(1,end-1) = -1;
        D2(2,end) = -1;
        % Divide each element by constant 12*dx^2
        D2 = D2./(12*dz^2);


        % Electron-Electron repulsion potential
        rho = zeros(length(PSI),1);
        for s = 1:Npsi 
            rho = rho + (abs(PSI(:,s)).^2);
        end
        fftRho = fft( [rho; zeros(N-1,1)]);      
        Vcoul2 =  ifft(fftRho.* fftCoul); 
        Vcoul  = Vcoul2(floor(N/2+1):floor(3*N/2));    
        Ue =.5*e^2*epsFactor*trapz( rho .* Vcoul ).*dz ;  
        if b3 == 0 
            Ue = 0;
        end

        % Electron kinetic energy
        Te = 0;
        for s = 1:Npsi 
            integrand = conj(PSI(:,s)).* ( D2* PSI(:,s) );
            Te = Te - spin*( hbar^2/(2*me) )*trapz(zDomHalf, integrand);
        end

        % Electron-ion bonding potential
        Ubond = 0;
        for s = 1:Nion;
            integrand = rho'.*1./sqrt(gamma^2 + ( zDomHalf - Rion(s) ).^2);  
            Ubond = Ubond - spin*Zion(s)*e^2*epsFactor*trapz(integrand)*dz;
        end


        % Ion kinetic energy
        Tion = 0;
        for s = 1:Nion;
            Tion = Tion + .5*Vion(s).^2*Mion(s);
        end
        
        % Ion/Ion repulsion potential
        Uion = 0;
        for s = 1:Nion   
            for ss = 1:Nion
                if s ~= ss;
                    const = .5*Zion(s)*Zion(ss)*e^2.*epsFactor*dz;
                    Enm   =  const./sqrt( gamma^2 + ( Rion(s) - Rion(ss) )^2 );
                    Uion = Uion + Enm;
                end
            end
        end

        Etot = (Te + Tion + Ue + Ubond + Uion)*jeV;

          text(zDomHalf(round(2.5*N/4)), .7, strcat( 'Total Energy = ', num2str(Etot), ' [s]'), 'fontsize', 18, 'fontname', 'cordia new')
          text(zDomHalf(round(2.5*N/4)), .8, strcat( 'Rho = ', num2str(rho(round(N/2))), ' [s]'), 'fontsize', 18, 'fontname', 'cordia new')
%         
%         
        text(zDomHalf(round(3*N/4)), .55, strcat( 'E.M. calculation = ', num2str(fieldTime), ' [s]'), 'fontsize', 18, 'fontname', 'cordia new')
        text(zDomHalf(round(3*N/4)), .5, strcat( 'ion calculation = ', num2str(ionTime), ' [s]'), 'fontsize', 18, 'fontname', 'cordia new')
        text(zDomHalf(round(3*N/4)), .45, strcat( 'psi calculation = ', num2str(psiTime), ' [s]'), 'fontsize', 18, 'fontname', 'cordia new')
%         
        text(zDomHalf(round(3*N/4)), .4, strcat( 'Density integral =', num2str(trapz(rho))), 'fontsize', 18, 'fontname', 'cordia new')
        text(zDomHalf(round(3*N/4)), .35, strcat( 'time step = ', num2str( r)), 'fontsize', 18, 'fontname', 'cordia new')
        text(zDomHalf(round(3*N/4)), .3, strcat( 'source = ', num2str( source/amp)), 'fontsize', 18, 'fontname', 'cordia new')
        text(zDomHalf(round(3*N/4)), .25, strcat( 'time = ', num2str((dt*r*1e15)), ' [fs]'), 'fontsize', 18, 'fontname', 'cordia new')
%         
% 
%         xlabel(strcat('z [nm]'), 'fontsize', 20, 'fontname', 'cordia new')
%         ylabel('Energy [eV] + RE{\Psi_i} [a.u.]', 'fontsize', 20, 'fontname', 'cordia new')

        axis( [ zDomHalf(1), zDomHalf(end),-1, 1] )
        getframe;
        
        clf
        hold off
    end
    %=====================================================================%

end


%%
figure('outerposition', [10   550   1400   550])

    % Amplitude in frequency domain
    magI = abs(EfInit).^2;
    magF = abs(EfFinal).^2;
    magR = abs(EfRef).^2;
    hold on
    
    plot(planck*testFreq/e, magR, 'r-')
    plot(planck*testFreq/e, magF, 'color', [.1 .1 .99], 'linewidth', 1)
    plot(planck*testFreq/e, magI, 'k-')
    
    legend('E_R', 'E_T', 'E_I')
    hold off
    ylabel('Power') 
    xlabel('{E = h\nu} [eV]')
    title('Spectrum')
    axis([planck*min(testFreq)/e, planck*max(testFreq)/e,  0 1.2*max([max(magR), max(magF), max(magI)])]) 
    
    
