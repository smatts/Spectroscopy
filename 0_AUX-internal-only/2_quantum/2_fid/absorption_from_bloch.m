% integrate coherence to get absorbed power

close all     % close all windows
clear variables     % clear all variables

% free parameters
Omega_Rabi = 5;       % Rabi frequency in inverse time units
%delta_omega = 0;       % (Atom - laser) frequency in inverse time units
pulse_area =  pi;       % radians
dt = 1e-3;             % time steps in time units

delta_omega_list = (-40:1:40);
pabs = zeros(size(delta_omega_list));

figure
hold on
xlabel('time (arb.u.)')
ylabel('coherence');

for id2 = 1:length(delta_omega_list)

    delta_omega = delta_omega_list(id2);
    M = [Omega_Rabi, 0,  delta_omega];  % torque vector
    S0 = [0,0,1];          % starting Bloch vector

    % calculate number of integration steps until pulse area is reached
    % here we define pulse are as on peak of resonance
    % we tune the laser, and keep laser power constant
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
%     Slength = vecnorm(Strace,2,2);
%     plot(Slength);

    coherence = 0.5 .* Strace(:,2);
    plot(coherence)

    pabs(id2) = -1 .* sum(coherence);
end

figure
plot(delta_omega_list, pabs)
xlabel('detuning (arb. u.)')
ylabel('absorbed power')

