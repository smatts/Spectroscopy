
% Lutz Langgut PRL acoustic Drexhage experiment

close all
clear all

c_sound = 343 ; % m/sec

f = 470; % Hz
lambda = c_sound / f;

d = (0:0.01:1);
x = 4 .* pi .* d / lambda;

eta = 0.5;

gamma_perp = 1 + eta .* (-3 .* sin(x) ./ x - 6 .* cos(x) ./ x.^2 + 6 .* sin(x) ./ x.^3);
gamma_par = 1 + eta .* (-3 .* cos(x) ./ x.^2 + 3 .* sin(x) ./ x.^3);

plot(d,gamma_perp,'r')
hold on
plot(d, gamma_par,'b')
hold off
