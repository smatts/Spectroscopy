close all
clear variables

data = load('../handout/Ewire1400_x_pol_new.dat');

% coordinates and electrical field
x = data(:,1);
y = data(:,2);
Ex = data(:,3) +1i * data(:,4);
Ey = data(:,5) +1i * data(:,6);
Ez = data(:,7) +1i * data(:,8); % is zero up to numerical accuracy

% choose a suitable time t / phase angle phi
% YOUR WORK here

phi = pi/2 ;
Ex = Ex .* exp(1i .*phi);
Ey = Ey .* exp(1i .*phi);
Ez = Ez .* exp(1i .*phi);

% calculate the direction (tx, ty) of the surface tangential 
% at all points (x,y)
% normalize the length of the tangential to 1
% hint: see circshift
% YOUR WORK here

tx = circshift(x,1) - circshift(x,-1);
ty = circshift(y,1) - circshift(y,-1);
L = sqrt(tx.^2 + ty.^2);
tx = tx ./ L;
ty = ty ./ L;

% calculate the direction (nx, ny) of the surface normal 
% at all points (x,y)
% YOUR WORK here

nx = -ty;
ny = tx;

%calculate the projection of the E-field on (n,t)
% YOUR WORK here

En = nx .* Ex + ny .* Ey;
Et = tx .* Ex + ty .* Ey;

% calculate the amplitude of p^(2) using only the (nnn) tensor element
% YOUR WORK here

chi2_nnn = 250;
p2nnn = chi2_nnn .* En .* En;

% calculate the (x,y) components of p^(2)
% YOUR WORK here

p2nnn_x = real(nx .* p2nnn);
p2nnn_y = real(ny .* p2nnn );


%plot the field
%adjust time t above so that the figure looks nice
% YOUR WORK here

quiver(x,y,p2nnn_x,p2nnn_y)


axis equal  % so that a circle is displayed as circle
title('Ewire1400 x pol new.dat')