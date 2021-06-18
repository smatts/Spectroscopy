close all, clear all

%% select data set
%   x_pol -> along long axis of structure (resonant)
%   y_pol -> along short axis of structure (off-resonant)


%data = load('Esplitring1750_x_pol_new.txt');
% data = load('Esplitring1750_y_pol_new.txt');
 data = load('Ewire1400_x_pol_new.txt');
% data = load('Ewire1400_y_pol.txt');

%% process data      

% get coordinates and local surface vectors
x = data(:,1);
y = data(:,2);
z = data(:,3);

nx = data(:,4);
ny = data(:,5);
nz = data(:,6);

tx = data(:,7);
ty = data(:,8);
tz = data(:,9);

sx = data(:,10);
sy = data(:,11);
sz = data(:,12);

% construct E-fields and add adjustable global phase (-> P2 pointing out of structure for nice plots)

phi = 90;
phase_factor = exp(1i*deg2rad(phi));

Ex = (data(:,13) + 1i* data(:,14)) * phase_factor;
Ey = (data(:,15) + 1i* data(:,16)) * phase_factor;
Ez = (data(:,17) + 1i* data(:,18)) * phase_factor;

normE = sqrt(Ex.*conj(Ex) + Ey.*conj(Ey) + Ez.*conj(Ez));

%% calculate P^(2)

% chi2-components
X1 = 250;
X2 = 0.* 1;
X3=0 .* 3.6;

% create nonlinear polarization in surface local coordinate system
Et = tx.*Ex + ty.*Ey + tz.*Ez;
En = nx.*Ex + ny.*Ey + nz.*Ez;
Es = sx.*Ex + sy.*Ey + sz.*Ez;
    
P2t = 2* X3*Et.*En;
P2n = X1*En.^2 + X2*(Es.^2 + Et.^2);
P2s = 2* X3*Es.*En;
    
% back transformation
P2x = nx.*P2n + tx.*P2t + sx.*P2s;
P2y = ny.*P2n + ty.*P2t + sy.*P2s;
P2z = nz.*P2n + tz.*P2t + sz.*P2s;

%% create figures - 3D data set

figure
quiver3(x,y,z,Ex,Ey,Ez,'LineWidth',1,'AutoScaleFactor',1,'Color',[0.6350 0.0780 0.1840])
axis equal
legend('E')
xlabel('x-coordinate (nm)')
ylabel('y-coordinate (nm)')
zlabel('z-coordinate (nm)')

figure
quiver3(x,y,z,P2x,P2y,P2z,'LineWidth',2,'AutoScaleFactor',2)
axis equal
legend('P^{(2)}')
xlabel('x-coordinate (nm)')
ylabel('y-coordinate (nm)')
zlabel('z-coordinate (nm)')

%% create figures - 2D data set
% reduce dataset to 2D and use circumference at the middle of the structure

% get z = 0 indices
ind = find(~z);

figure
quiver(x(ind),y(ind),Ex(ind),Ey(ind),'LineWidth',1,'AutoScaleFactor',1,'Color',[0.6350 0.0780 0.1840])
axis equal
legend('E')
xlabel('x-coordinate (nm)')
ylabel('y-coordinate (nm)')

figure
quiver(x(ind),y(ind),P2x(ind),P2y(ind),'LineWidth',2,'AutoScaleFactor',2)
axis equal
legend('P^{(2)}')
xlabel('x-coordinate (nm)')
ylabel('y-coordinate (nm)')


%% emission calculation for d << lambda

P2total = [sum(P2x) sum(P2y) sum(P2z)]';
Ishg = sum(P2total.*conj(P2total));

disp(['estimated far field emission ' num2str(Ishg)])














