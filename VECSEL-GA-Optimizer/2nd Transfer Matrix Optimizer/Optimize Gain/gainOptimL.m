function [negGain] = gainOptimL(lambda, n, x, lWell, numWell, sigma, scaleFactor)
global storage
global storageInit
global pop

% Grating in active region

lH = x(2)/scaleFactor;                    % Widths of surrounding cladding (High refractive index)
lL = x(1)/scaleFactor;                    % Widths of surrounding cladding (Low refractive index)
LActive = repmat([ lH , lWell, lL], 1, numWell);
L = LActive;

%Constants
cNot = 2.99792458e8;            % Speed of light
epsNot = 8.854187817620e-12;    % Permittivity of Free Space

% solution with gain
[R, r, T, t] = transferMatrixOne(lambda, n, L, sigma);

%%
gain = (R).*100;
negGain= -gain;

if length(storageInit) < pop;
    storageInit = [x(1)/scaleFactor, x(2)/scaleFactor; storageInit];
else
    storage = [x(1)/scaleFactor, x(2)/scaleFactor; storage];
end

% 
% storage
% storageInit

end