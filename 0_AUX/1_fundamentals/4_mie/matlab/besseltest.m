
close all
clear all

x = (0:0.01:12);

%    zj = besselj( nhelp + 0.5, xhelp);
%     zh = besselj( nhelp + 0.5, xhelp) + 1i .* bessely( nhelp + 0.5, xhelp);

j0 = sqrt(pi ./ (2 .* x)) .* besselj(0 + 0.5, x);

j1 =  sqrt(pi ./ (2 .* x)) .* besselj(1 + 0.5, x);

hold on
plot(x, j0)
plot(x, j1)
hold off