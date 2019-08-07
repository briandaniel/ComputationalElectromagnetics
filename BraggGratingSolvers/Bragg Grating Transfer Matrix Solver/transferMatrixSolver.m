close all;
clear all;

cNot = 2.99793e8;                    % Speed of Light
lambdaMax = 500e-9;                 % Central wavelength (reflectance max here)
lowerLambda = 0.60*lambdaMax;        % Lower limit of wavelengths
upperLambda = 1.80*lambdaMax;        % Upper limit of wavelengths
lambdaSteps = 500;                   % Number of steps to use in wavelength vector
lambdas = linspace(lowerLambda, upperLambda, lambdaSteps);

% Parameters
nH = 2.32;
nL = 1.38;
LH = .25;
LL = .25;
lH = lambdaMax*LH/(nH);
lL = lambdaMax*LL/(nL);
rho = (nH - nL)/(nH + nL);
la2 = pi*(LL+LH)*1/acos(rho)*lambdaMax;
la1 = pi*(LL+LH)*1/acos(-rho)*lambdaMax;


grateNum = 13;                        % Number of grating sets
grateWidth = lH*grateNum + lL*grateNum;
L = grateWidth;                      % Total Grating length (meters)
dz = lH;                    % Individual section length
grateNum = 2*grateNum+1;

dzo = - L + 10.69e-3;    %????

del = lH+lL;
neff = lambdaMax/(2*del);
% neff = 1.32;
% neff = 2.00;
% del = lambdaMax/(2*neff);             % Grating period
qL = 10;                               % Normalization constant

qo = pi/(2*neff*del)*.55;

% Vector initilization
TM = eye(2,2);
T = ones(2,2);
gdc = 1;
r = zeros(1,lambdaSteps);
t = zeros(1,lambdaSteps);
R = r;
TP = t;

for k = 1:lambdaSteps;
    TM = eye(2,2);
    lambda = lambdas(k);
    sigma = 2*pi*neff*(1/lambda - 1/lambdaMax);
    
    for m = 1:grateNum

    %  Transfer matrix for the particular section
        if mod(grateNum,2) == 0;
            dz = lL;
            q = qo*(nL-nH);
        else
            dz = lH;
            q = qo*(nH-nL);
        end
        gamma = sqrt((abs(q))^2 - gdc^2);
        delta = sigma;

        T= [cosh(gamma*dz)-i*delta/gamma*sinh(gamma*dz),        -i*q/gamma*sinh(gamma*dz); ...
               i*q/gamma*sinh(gamma*dz),            cosh(gamma*dz)+i*delta/gamma*sinh(gamma*dz)];
    
    %   Update overall transfer
        TM = T*TM;

    end

    r(k) = TM(2,1)/TM(1,1);
    R(k) = (abs(r(k)))^2 ;        %Power Reflectivity

    %????
    PHI = 2*pi*neff*dzo/k;
    %phase difference between gratings
    Fp = [exp(-i*PHI) 0; 0 exp(i*PHI)];
    Ffp = TM*Fp*TM;
    t(k) = 1/Ffp(1,1);
    TP(k) = (abs(t(k)))^2;
    
end


plot(lambdas.*10^9, R)
% plot(lambdas.*10^9, TP)
legend('Reflection', 'Transmission')
ylabel('Power')
xlabel('Wavelength(nm)')
text(700,.7,strcat('Grate size = { }',  num2str(round(L*1e7)*1e-1), ' {\mu}m'))
text(700,.64,strcat('Number of grates = { }',  num2str(grateNum)))
text(700,.56,strcat('{\lambda}_o= { }',  num2str(round(lambdaMax*1e10)*1e-1), ' nm'))
text(700,.5,strcat('n_e_f_f = ',  num2str(round(neff*1e3)*1e-3)))

axis([min(lambdas)*10^9 max(lambdas)*10^9 0 1])
