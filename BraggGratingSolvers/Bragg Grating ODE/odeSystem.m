function [ dS ] = odeSystem(t, S, kappa, sigma )
% Function of two simple ODE's.

%  S(1) = R          S(2) = S
% dS(1) = dR/dz     dS(2) = dS/dz


dS = zeros(2,1);
dS(1) = sigma*i*S(1) + kappa*i*S(2);
dS(2) = -sigma*i*S(2) - kappa*i*S(1);

% dS = -i*S;

end

%%
