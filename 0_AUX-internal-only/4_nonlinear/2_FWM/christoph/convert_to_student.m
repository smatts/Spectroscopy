close all
clear all

load 'meas_4.5abde.mat'

wl = meas{1}.wl;
linphase = meas{1}.linphase;

d = meas{4}.fwm;

save('fwm_wide_2.dat','d','-ascii')

save('wl.dat','wl','-ascii');
save('linphase.dat','linphase','-ascii');

% excite = interp1(meas.exc.wl, meas.exc.int, mywl);
% excite = 100 .* excite ./ max(excite);
% excite(isnan(excite)) = 0;
% 
% fwm = interp1(meas.fwm.wl, meas.fwm.int, mywl);
% fwm = 100 .* fwm ./ max(fwm);
% fwm(isnan(fwm)) = 0;
% 
% plot(mywl, fwm)
% hold on
% plot(mywl, excite)
% hold off
% 
% data_out = [mywl; excite; fwm]'
% 
