close all
clear all
data = import_bodipy('data/BODIPY_650_665.csv');

wl = data.Wavelength;
absorption = data.Excitation;
fluo = data.Emission;


fluo_en = fluo .* wl.^2;
fluo_en = fluo_en ./ max(fluo_en);
absorption = absorption ./ max(absorption);

energy = (1.5:0.01:2.5);
abs_en = interp1(1240 ./ wl, absorption, energy);
fl_en = interp1(1240 ./ wl, fluo_en, energy);

abs_en = abs_en ./ energy;
fl_en = fl_en ./ energy.^3;
abs_en = abs_en ./ max(abs_en);
fl_en = fl_en ./ max(fl_en);


[~, maxid_abs] = max(abs_en);

[~, maxid_fl] = max(fl_en);

pre = 30;
post = 20;

energy_cut = energy(maxid_fl - pre: maxid_fl + post) - energy(maxid_fl);
fl_cut = fl_en(maxid_fl - pre: maxid_fl + post);
abs_cut = fliplr(abs_en(maxid_abs - post: maxid_abs + pre));

plot(energy_cut, abs_cut)
hold on
plot(energy_cut, fl_cut)
hold off

legend('abs','em')