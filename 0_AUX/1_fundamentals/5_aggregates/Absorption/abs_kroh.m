close all
clear all

d = importdata('blank_NaOH_abs.csv');
data_blank = d.data;

d = importdata('TDBC_NaOH_c2_1mmK_abs-2.csv');
data_c2 = d.data;

d = importdata('TDBC_NaOH_c3_abs.csv');
data_c3 = d.data;

d = importdata('TDBC_NaOH_c4_abs.csv');
data_c4 = d.data;

d = importdata('TDBC_NaOH_c5_abs.csv');
data_c5 = d.data;


conc = [0, 1e-5, 5e-6, 1e-6, 5e-7];
%conc = conc .* 0 + 1;

wl_fix = 790;
[~ , id_fix] = min(abs(data_c3(:,1) - wl_fix));
[~ , id_fix_c2] = min(abs(data_c2(:,1) - wl_fix));

plot(data_c2(:,1),(data_c2(:,2) - data_c2(id_fix_c2, 2)).* 10 ./ conc(2))

hold on
plot(data_c3(:,1),(data_c3(:,2) - data_c3(id_fix,2)) ./ conc(3))
plot(data_c4(:,1),(data_c4(:,2) - data_c4(id_fix,2))./ conc(4) )
plot(data_c5(:,1),(data_c5(:,2) - data_c5(id_fix,2)) ./ conc(5))

%plot(data_blank(:,1), (data_blank(:,2) - data_blank(id_fix,2)) )

hold off