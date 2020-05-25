close all
clear all

load 'meas_4.3b.mat'

mywl = (600:0.25:850);

excite = interp1(meas.exc.wl, meas.exc.int, mywl);
excite = 100 .* excite ./ max(excite);
excite(isnan(excite)) = 0;

fwm = interp1(meas.fwm.wl, meas.fwm.int, mywl);
fwm = 100 .* fwm ./ max(fwm);
fwm(isnan(fwm)) = 0;

plot(mywl, fwm)
hold on
plot(mywl, excite)
hold off

data_out = [mywl; excite; fwm]';
save('fwm_spectrum.dat','data_out','-ascii')

