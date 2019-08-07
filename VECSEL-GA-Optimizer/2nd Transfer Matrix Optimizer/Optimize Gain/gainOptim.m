function [gainInv] = gainOptim(lambda, n, L, sigma, scale)

L = L/scale;

%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

% solution with gain
[R, r, T, t] = transferMatrixOne(lambda, n, L, sigma);

%%
gain = (R).*100;
gainInv = -gain;

end