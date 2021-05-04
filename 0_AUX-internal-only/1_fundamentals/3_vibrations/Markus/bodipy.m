close all
clear all

data = importfile('BODIPY_650_665.csv');
wl = data.Wavelength;
absspec = data.Excitation;
absspec = absspec ./ max(absspec);

emspec = data.Emission;

wn = 1e-2 ./ (wl .* 1e-9);
emspec_wn = emspec .* wl.^2;
emspec_wn = emspec_wn ./ max(emspec_wn);


wnx = linspace(1.2e4, 2e4, 1000)';
absspec_x = interp1(wn, absspec, wnx);
emspec_x = interp1(wn, emspec_wn, wnx);



plot(wnx, absspec_x);
hold on
plot(wnx, emspec_x);
hold off

[m,ida ] = max(absspec_x);
[m,idb ] = max(emspec_x);

id_mitte = fix((ida + idb) ./ 2);
wnxm = 2 .* wnx(id_mitte) - wnx;

pa =  absspec_x ./ wnx;
pa = pa ./ max(pa);

pf =  emspec_x ./ wnx.^3;
pf = pf ./ max(pf);


figure
plot(wnx, pa);
hold on
plot(wnxm, pf);
hold off
xlabel('wavenumber absorption (cm^{-1})')
ylabel('dipole moment (normalized)')
legend('absorption','fluorescence')
xlim([1.3e4, 1.9e4])
