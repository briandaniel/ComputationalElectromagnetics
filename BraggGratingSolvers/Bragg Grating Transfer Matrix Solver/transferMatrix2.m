function [R,T] = transferMatrix(lambdaNot, lambdaMin, lambdaMax, n, L)


lambdas = linspace(lambdaMin, lambdaMax, 501);
n = [1 2.32 1.38 2.32 1.38 2.32 1.52];
L = [.1078 .1812 .1078 .1812 .1078].*500e-9;
grateLength = sum(L);
cellNum = length(L);

Mi = zeros(2,2);
r = zeros(length(lambdas),1);
t = zeros(length(lambdas),1);


for w = 1:length(lambdas);
    
    lambda = lambdas(w);
    ko = 2*pi*n(1)/lambda;
    kf = 2*pi*n(end)/lambda;
    M = eye(2,2);
    
    for s = 1:cellNum;
        
        kl = 2*pi*n(s+1)*L(s)/lambda;
        k = kl/(L(s));
        Mi = [cos(kl), 1/k*sin(kl); -k*sin(kl), cos(kl)];
        M =  Mi*M;
        
    end
        
    r(w) = ((M(2,1)+ko*kf*M(1,2))+1i*(ko*M(2,2)-kf*M(1,1)))...
           /((-M(2,1)+ko*kf*M(1,2))+1i*(ko*M(2,2)+kf*M(1,1)));
       
    t(w) = 2*1i*ko*exp(-1i*kf*grateLength)*(M(1,1)*M(2,2)-M(1,2)*M(2,1))...
           /(-M(2,1)+ko*kf*M(1,2)+1i*(kf*M(1,1)+ko*M(2,2)));
    
end

R = abs(r).^2;
T = abs(t).^2;

plot( lambdas*1e9, R, lambdas*1e9, T, lambdas*1e9)


