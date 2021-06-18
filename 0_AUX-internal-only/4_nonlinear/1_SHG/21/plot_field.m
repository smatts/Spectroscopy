close all
clear variables

data = load('../handout/Ewire1400_x_pol_new.dat');

% coordinates and electrical field
x = data(:,1);
y = data(:,2);
Ex = data(:,3) +1i * data(:,4);
Ey = data(:,5) +1i * data(:,6);
Ez = data(:,7) +1i * data(:,8); % is zero up to numerical accuracy

% calculate field vector in 2D (ignore z component)
% choose a suitable time t / phase angle phi
%YOUR WORK here 

phi = pi/4;
vx = real(Ex .* exp(1i .* phi));
vy = real(Ey .* exp(1i .* phi));

%plot the E-field as vector field
%YOUR WORK here 
quiver(x,y,vx,vy)

axis equal  % so that a circle is displayed as circle
title('Ewire1400 x pol new.dat')