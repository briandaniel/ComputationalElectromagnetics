function [sums] = TROptimSum(lambda, n, L, sigma)

%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

% solution with gain
[R, r, T, t] = transferMatrixOne(lambda, n, L, sigma);
Lrev = L(end:-1:1);
nrev = n(end:-1:1);
[Rrev, rrev, Trev, trev] = transferMatrixOne(lambda, nrev, Lrev, sigma);

%%
TInv = 1/(T*100);
RrevInv = 1/(Rrev*100);
sums = 1/((T*3.4+Rrev)*100);
end