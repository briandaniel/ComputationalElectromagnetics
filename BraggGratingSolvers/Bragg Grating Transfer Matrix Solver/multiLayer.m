function [ la, Gla, r ] = multiLayer(la0, lambdaMin, lambdaMax, lambdaSteps, nH, nL, N)
%MULTILAYER Summary of this function goes here
%   Detailed explanation goes here


na = 1; 
nb = 1;                                 % refractive indices
LH = 0.25; LL = 0.25;                           % optical thicknesses in units of λ0

la0;                                            % λ0 in units of nm
rho = (nH-nL)/(nH+nL);                          % reflection coefficient ρ
la2 = pi*(LL+LH)*1/acos(rho) * la0;             % Right bandedge
la1 = pi*(LL+LH)*1/acos(-rho) * la0;            % left bandedge
Dla = la2-la1;                                  % bandwidth

N;                                          
n = [na, nH, repmat([nL,nH], 1, N), nb];        % indices for the layers A|H(LH)N |G 
L = [LH, repmat([LL,LH], 1, N)];                % Lengths of layers
la = linspace(lambdaMin,lambdaMax,lambdaSteps); % plotting range is 300 ≤ λ ≤ 800 nm
r = multidiel(n,L,la/la0);
Gla = 100*abs(r).^2; 

% 
% figure; plot(la,Gla); 
% f = linspace(0,6,1201); 
% Gf = 100*abs(multidiel(n,L,1./f)).^2; 
% figure; plot(f,Gf); 

end

