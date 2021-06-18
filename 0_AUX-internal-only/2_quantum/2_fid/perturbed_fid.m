% perturbed free inductiond decay
close all
clear variables


dt = 0.1;  % time step in time units
tmax = 1000;  % max simualtion time
fsample = 1/ dt;  % sample frequency in inverse time units

time = (0:dt:tmax) ;   % time axis

% list of pump-probe delay values that are used in the final plot
pp_delay_list = (-25:0.1:25);  

beta_pp = 0.5;  % efficiency of the pump-pulse (1 = max, 0=no effect)
phi_qd = 0.4 .* pi; % phase difference between reflected and scattered field
gamma = 0.05;   % decay constant of the coherence, in inverse time units

%--------------------
% calculate reference spectrum, i.e., without pump pulse
reference = exp(-1 .* gamma .* time);  % coherence decays exponentially
ft_ref = fftshift(fft(reference)) ;  % numerical Fourier trasnform 
spec_ref = real(exp(1i .* phi_qd) .*  ft_ref);  % spectrum is real part

%-------------------
% loop over pump-probe delay
data = zeros( length(pp_delay_list),length(time));  % store data

for id = 1:length(pp_delay_list)
    
    pp_delay = pp_delay_list(id);  % current delay
    
    % generate effect of pump pulse: step function at time = pp_delay
    pump = zeros(size(time)) ;
    pump(time > -1 .* pp_delay) = beta_pp;
    
    % signal is reduced by pump pulse
    signal = reference .* (1 - pump);
    ft_sig = fftshift(fft(signal)) ;  % FT
    spec_sig = real(exp(1i .* phi_qd) .* ft_sig); % spectrum
    
    data(id,:) = spec_sig - spec_ref; % we measure difference in spectrum
end

% build frequency axis
n = length(signal);
n_half = fix(n/2); % round to next interger below n/2
n_remainder = n - 2 * n_half;
omega = 2 .* pi .* (-1 * n_half :n_half -1 + n_remainder)./ n .* fsample;

% plot data as matrix / flase color plot
imagesc(pp_delay_list, omega, transpose(data))
xlabel('pp delay')
ylabel('frequency')
ylim([-2,2])
colormap('gray')
c = colorbar;
c.Label.String = 'delta P(omega)';

