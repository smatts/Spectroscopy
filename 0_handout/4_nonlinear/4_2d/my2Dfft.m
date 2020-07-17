function [data_out, omega] = my2Dfft(data_in, time_in)
% 2D FFT with zero padding for t < 0
%        and omega axis

% we know that there is no signal (in this Feynman diag) 
% for a different time order = negtaive times
% here we include this information by appending zeros

% time starts at zero for the first pixel and then increases
% t<0 pixel are thus in the second half of the array

    zero = data_in .* 0;
    padded = [data_in, zero; zero, zero];


    data_out =fftshift(fft2(padded));


    % build frequency axis
    n = 2 .* length(time_in);
    n_half = fix(n/2); % round to integer
    n_remain = n - 2 * n_half;
    fs = abs(1./(time_in(2) - time_in(1))); % sampling frequency
    freq = (-1 * n_half :1:n_half -1 + n_remain)./ n .* fs;
    omega = freq .* 2 .* pi;


end
