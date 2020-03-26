%ensemble absorption and emission in PMMA film
load('TDI_absem.mat');
%typical single molecule tcspc trace
load('TDI_TCSPC.mat');

figure(1)
subplot(2,1,1)
hold on
plot(lambda_emission,emission,'Color',[1 0 0]);
plot(lambda_absorption, absorption,'Color',[0 0 1]);
axis([500 850 0 1.1]);
ylabel('absorption & emission');
xlabel('\lambda / nm');
box on
pbaspect([2 1 1]);

subplot(2,1,2)
plot(time_tcspc,emission_tcspc,'Color',[1 0 0]);
axis([-1 10 0 300]);
ylabel('counts / bin');
xlabel('time / ns');
box on
pbaspect([2 1 1]);
