close all;
clear all;

N = 10000;
L = 2*10^-3;
kappa = 8/L;
sigma = (20:-.2:-20)/L;
rho = zeros(1,length(sigma));
options = odeset('RelTol',1e-4,'AbsTol',1e-4);

% sigmahat*L is between -20 and 20;
% kappa*L = 8;


for n = 1:length(sigma);

    meshSize = L/100;
    domain = bvpinit(0:meshSize:L,[1, .5]);
    sol = bvp4c(@(T,Y)odeSystem(T,Y,kappa,sigma(n)), @bounds, domain, options);

    y1Initial = sol.y(1,1);
    y2Initial = sol.y(2,1);

    y1Sltn = y1Initial;
    y2Sltn = y2Initial;
    
    rho(n) = y2Sltn/y1Sltn;
    
end

r = abs(rho).^2;
normalWave = 1./(1+(sigma.*L)/(pi*N));

plot(normalWave, r)

%%% Plots %%%
% figure
% set(1,'outerposition', [1500 200 900 800])
% 
% subplot(2,1,1)
% plot(sol.x, sol.y(1,:), sol.x, sol.y(2,:))
% legend('R','S')
% 
% subplot(2,1,2)
% plot(sol.x, sol.yp(1,:), sol.x, sol.yp(2,:))
% legend('Rprime','Sprime')
