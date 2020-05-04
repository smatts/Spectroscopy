
close all
clear all

d80 = importdata('200427_BBI-AuColloid80nm.csv');
d80_2 = importdata('200427_BBI-AuColloid80nm_2.csv');

rod550 = importdata('200427_NP_rod_30-25-550-25.csv');
rod650 = importdata('200427_NP_rod_30-25-650-25.csv');
rod850 = importdata('200427_NP_rod_30-10-850-25.csv');
rod780 = importdata('200427_NP_rod_A12-40-780.csv');

np30 = importdata('200427_NP_sp_11-30-25.csv');
np50 = importdata('200427_NP_sp_12-50-25.csv');
np70 = importdata('200427_NP_sp_13-70-25.csv');

water = importdata('200427_DI-Water.csv');
water = water(1:601,:) ;
%plot(rod650(:,1), rod650(:,2)- water(:,2))

%plot(d80(:,1), d80(:,2) - water(:,2))


wl = np30(:,1);
np30 =  np30(:,2) - water(:,2);
np50 =  np50(:,2) - water(:,2);
np70 =  np70(:,2) - water(:,2);
np80 =  d80_2(:,2) - water(:,2);

hold on
 plot(wl, np30 ./ max(np30))
 plot(wl, np50 ./ max(np50))
 plot(wl, np70 ./ max(np70))
 plot(wl, 1.2 .* np80 ./ max(np80))
hold off

xlim([400, 1000])

figure

hold on
plot(rod550(:,1), rod550(:,2) - water(:,2))
plot(rod650(:,1), rod650(:,2)- water(:,2))
plot(rod780(:,1), rod780(:,2)- water(:,2))

 plot(rod850(:,1), rod850(:,2) - water(:,2))
 hold off
