% plot trajectory of Bloch vector on sphere

close all     % close all windows
clear variables     % clear all variables

% free parameters
Omega_Rabi = 10;       % Rabi frequency in inverse time units
delta_omega = 0;       % (Atom - laser) frequency in inverse time units
pulse_area = pi;       % radians
dt = 1e-3;             % time steps in time units

M = [Omega_Rabi, 0,  delta_omega];  % torque vector
S0 = [0,0,1];          % starting Bloch vector

% calculate number of integration steps until pulse area is reached
% pulse_area = N_steps * dt * Omega_Rabi
N_steps = floor(pulse_area ./ (dt .* Omega_Rabi)) ;

Strace = zeros(N_steps,3);   % here we collect the trace of the Bloch vector

S = S0;  % we start
for id=1:N_steps
    S = S + dt .* cross(M, S) ;  % next Bloch vector
    Strace(id,:) = S;   % store result
end

% check length of S
% |S|^2 at each time step
Slength = vecnorm(Strace,2,2);
plot(Slength);


% plot Bloch sphere
figure
plot3(Strace(:,1),Strace(:,2),Strace(:,3), 'linewidth',3)
hold on   % more to come

% draw direction of torque vector
Mdir = M ./ norm(M) .* 1.2;   % make vector stick out over sphere
Mvec = [ 0  0 0 ; Mdir];
plot3(Mvec(:,1),Mvec(:,2),Mvec(:,3), 'linewidth',2)

% draw lon / lat circles
phi = linspace(0, 2* pi, 100);  % list of angles
for theta = ((1:4).* pi/4)
    [x,y,z] = sph2cart(theta,phi,1); % convert spherical to cartesian
    plot3(x,y,z,'k')   % plot in thin black
    
    [x,y,z] = sph2cart(phi,theta/2 - pi/2, 1); 
    z = y .* 0 + z;   % z is single value, but needs go be vector
    plot3(x,y,z,'k')
    
    [x,y,z] = sph2cart(phi,theta/2 , 1); 
    z = y .* 0 + z;
    plot3(x,y,z,'k')

end

% adjust plot range so that sphere is spherical
mylim = [-1,1] .* 1.1;
xlim(mylim)
ylim(mylim)
zlim(mylim)
axis square
xlabel 'u'
ylabel 'v'
zlabel 'w'


figure
plot(Strace(:,1))
hold on
plot(Strace(:,2))
plot(Strace(:,3))
hold off
legend('u','v','w')
xlabel 'time (arb.u.)'


