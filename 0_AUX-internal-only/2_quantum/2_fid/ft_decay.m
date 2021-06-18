% the FT of an exponential decay is the Lorentz function

close all
clear all


dt = 0.01;
tmax = 10;
fsample = 1/ dt;

time = (0:dt:tmax) ;


gamma = 2;
signal = exp(-1 .* gamma .* time);
signal(time < 0) = 0;


ft = fftshift(fft(signal)) ;

% bastle Frequenzachse
n = length(signal);
n_halbe = fix(n/2); % runde ab auf ganze Zahl
n_rest = n - 2 * n_halbe;
omega = 2 .* pi .* (-1 * n_halbe :n_halbe -1 + n_rest)./ n .* fsample;


lorentz = 1 ./ (gamma + 1i .* omega);

plot(omega, real(ft) .* dt)
hold on
plot(omega, real(lorentz) )
hold off

figure
plot(omega, imag(ft) .* dt)
hold on
plot(omega, imag(lorentz) )
hold off


