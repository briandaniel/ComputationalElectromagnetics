function [ dPsi ] = nls( t, Psi, k, Vx, c1, c2, c3) %#ok<INUSL>
% Solves the NLS using a pseudospectral fft method
% input requires parameters (k, Vx, c1, c2, c3)
% i.e. time dependent schrodinger equation with linear potential and 
% nonlinear square term...
%
%  dPsi      d^2Psi
%  ---- = c1*------ + c2*V(x) + c3*|Psi|^2*Psi
%   dt        dx^2
% 
    uxx = (ifft( (k*1i).^2 .* fft(Psi) ) );
    uNL = abs(Psi).^2.*Psi;    
    dPsi = c1.*uxx + c2*Vx.*Psi + c3*uNL;

end

