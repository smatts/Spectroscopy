close all
clear all

omega_1 = 1;

k_by_m = 0.1;

omega_2 = (0.5:0.01:1.5);

a = (omega_1 .^ 2 + omega_2 .^ 2) ./ 2 + k_by_m;

b = ((omega_1 .^ 2 - omega_2 .^ 2) ./ 2).^2 + k_by_m.^2;

omega_p = sqrt(a + sqrt(b));
omega_m = sqrt(a - sqrt(b));

hold on
plot(omega_2, omega_p)
plot(omega_2, omega_m)
plot(omega_2, omega_2 .* 0 + 1);
%plot(omega_2, omega_2 .* 0 + 1 + k_by_m ./ 2);
plot(omega_2, omega_2 );

hold off
