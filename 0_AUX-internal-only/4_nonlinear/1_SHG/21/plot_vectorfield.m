close all
clear variables

% plot a vector field, here in 2D

phi = (0:5:360)./ 180 .* pi;

% coordinates of points at which we know the field
x = cos(phi);
y = sin(phi);

% field at these points
vx = cos(phi);
vy = sin(phi);

% plot the field
% YOUR WORK here 
 quiver(x,y,vx,vy)

axis equal  % so that a circle is displayed as circle