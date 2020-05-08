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

wl_fix = 650;
[~ , id_fix] = min(abs(data_c3(:,1) - wl_fix));
[~ , id_fix_c2] = min(abs(data_c2(:,1) - wl_fix));

plot(data_c2(:,1),(data_c2(:,2) - data_c2(id_fix_c2, 2)).* 10 ./ conc(2))

hold on
plot(data_c3(:,1),(data_c3(:,2) - data_c3(id_fix,2)) ./ conc(3))
plot(data_c4(:,1),(data_c4(:,2) - data_c4(id_fix,2))./ conc(4) )
plot(data_c5(:,1),(data_c5(:,2) - data_c5(id_fix,2)) ./ conc(5))

%plot(data_blank(:,1), (data_blank(:,2) - data_blank(id_fix,2)) )

hold off

%-------------
wl1 = 650;
[~ , id1a] = min(abs(data_c2(:,1) - wl1));
wl2 = 570;
[~ , id2a] = min(abs(data_c2(:,1) - wl2));

a2 = data_c2(id1a:id2a,2) - data_c2(id_fix_c2, 2);
a2 = a2 ./ max(a2);


[~ , id1] = min(abs(data_c3(:,1) - wl1));
[~ , id2] = min(abs(data_c3(:,1) - wl2));

a3 = data_c3(id1:id2,2) - data_c3(id_fix,2);
a4 = data_c4(id1:id2,2) - data_c4(id_fix,2);
a5 = data_c5(id1:id2,2) - data_c5(id_fix,2);
a3 = a3 ./ max(a3);
a4 = a4 ./ max(a4);
a5 = a5 ./ max(a5);

plot(data_c2(id1a:id2a,1), a2)
hold on

plot(data_c3(id1:id2,1), a3)
plot(data_c3(id1:id2,1), a4)
plot(data_c3(id1:id2,1), a5)
hold off








