close all;
clear all;

options = odeset('RelTol',1e-6,'AbsTol',1e-8);

L = 4*10^-3;


Rinit = 1;                      % Initial R value
Sinit = -0.1030 + 0.9193i;      % Initial S value
L = 4*10^-3;                    % Final z value

[T,Y] = ode45(@odeSystem,[0 L],[Rinit Sinit],options);



%%% Plots %%%
figure
set(1,'outerposition', [1500 200 900 800])

subplot(2,1,1)
plot(T, Y(:,1))
title('R')

subplot(2,1,2)
plot(T, Y(:,2))
title('S')

% figure
% subplot(2,1,1)
% plot(T, imag(Y(:,1)))
% title('Im R')
% 
% subplot(2,1,2)
% plot(T, imag(Y(:,2)))
% title('Im S')