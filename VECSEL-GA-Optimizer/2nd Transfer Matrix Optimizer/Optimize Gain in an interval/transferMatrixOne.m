function [R, r, T, t] = transferMatrixOne(lambda, n, L, sigma)
%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

%Initilization
grateLength = sum(L);
cellNum = length(L);

Mi = zeros(2,2);
r = zeros(length(lambda),1);
t = zeros(length(lambda),1);

eps = n.^2;

ko = 2*pi*n(1)/lambda;
kf = 2*pi*n(end)/lambda;
M = eye(2,2);
freq = cNot/lambda;
eps = n.^2 + i.*sigma./(freq.*2*pi*epsNot);
nW =  sqrt(eps);

for s = 1:cellNum;

    kl = 2*pi*nW(s+1)*L(s)/lambda;
    k = kl/(L(s));
    Mi = [cos(kl), 1/k*sin(kl); -k*sin(kl), cos(kl)];
    M =  Mi*M;

end

r = ((M(2,1)+ko*kf*M(1,2))+1i*(ko*M(2,2)-kf*M(1,1)))...
       /((-M(2,1)+ko*kf*M(1,2))+1i*(ko*M(2,2)+kf*M(1,1)));


t = 2*1i*ko*exp(-1i*kf*grateLength)*(M(1,1)*M(2,2)-M(1,2)*M(2,1))...
       /(-M(2,1)+ko*kf*M(1,2)+1i*(kf*M(1,1)+ko*M(2,2)));

R = abs(r).^2;
T = abs(t).^2;
