close all
clear all
clc


% everything in column vectors N x 1
lambda = transpose(linspace(400,900,300)) ;    % wavelength (nm)
r_particle = 50;  % particle radius
nenv = 1.33; % effictive refractive inde1x


%------------------------
% calculate Mie absorption cross section (Nenv is somewhere between glas and air)
% the "miegoldabs" takes and gives everything in units of nm (nm^2 ...)

nqu=load('jc_gold.dat');
energ=nqu(:,1);     %Energie in eV
n=nqu(:,2);    %Brechungsindex Gold
k=nqu(:,3);    %Absorption Gold
wel=4.13566733*2.99792458e-1 ./ energ*1000;  



% Interpolation
ior1 = interp1(wel,n, lambda,'spline');  
ior2 = interp1(wel,k, lambda,'spline');

ior = ior1 + 1i .* ior2;


[sigma_ext, sigma_sca] = mie(ior1 + 1i .* ior2, lambda , nenv, r_particle);

 
%------

eps_in = ior .^2;
eps_out = nenv .^2;

alpha = 4 .* pi .* r_particle.^3 .* (eps_in - eps_out) ./ ( eps_in + 2 .* eps_out);

k = 2 .* pi .* nenv ./ lambda;

sigma_ray = k .* imag(alpha);

hold on
plot(lambda, sigma_ext)
plot(lambda, sigma_sca)
plot(lambda, sigma_ext - sigma_sca)

plot(lambda, sigma_ray ,'r')
hold off

%---------------
