clear;
clf();
close all;

%Konstanten
    NA = 6.022e23; %[1/mol] Avogadrokonstante
    c0 = 299792458;%[m/s] Lichtgeschwindigkeit
    n = 1;      %Brechungsindex Wasser
    c = c0/n;      %[m/s] Lichtgeschwindigkeit in Wasser
    h = 6.626e-34; %[Js] Planck-Konstante
    hq = h/(2*pi); %[Js]
    e = 1.602e-19; %[C] Elementarladung
    
%allgemeine Plot-Optionen
    col_plot = [0, 0, 0];
    col_fit = [0.8500, 0.3250, 0.0980];
    fontsize = 14;
    linewidth = 1;

%pump probe delay tau
    tau_min = -44e-12;
    tau_max = 12e-12;
    dtau = 1e-12;
    tau_array = tau_min:dtau:tau_max;
    
%Loop pump probe delay tau
for k=1:length(tau_array)
    tau = tau_array(k);
    disp(tau);
    %% Parameter
        S(1,1) = 1; %Bloch-Vektor (unnormiert)
        S(1,2) = 0;
        S(1,3) = 0;

        %Normierung
            sumS = sqrt(S(1,1)^2+S(1,2)^2+S(1,3)^2);
            S(1,1) = S(1,1)/sumS;
            S(1,2) = S(1,2)/sumS;
            S(1,3) = S(1,3)/sumS;

        dt = 1e-16;         %Zeitschritt
        tmax = 100e-12;     %Dauer der Simulation
        dRho = 0.5;         %durch Pump-Puls verursachte Störung von rho00
        dphi = -0.4*pi;     %Phasenshift zwischen E0 und Es
        tau_relax = 10e-12;  %Relaxationszeit des durch Pump-Puls angeregten Zustands

        lambL = 1519e-9;   %Laserwellenlänge
        lamb0 = 1519e-9;   %Resonanzwellenlänge (Wasser: ~1519 nm)
        pulseLength = 1e-12; %Dauer Probe Puls

        mu = 6.152e-30; %Übergangsdipolmoment
        EE = 1e-4;  %Feldstärke Probe

        wL = 2*pi*c/lambL;  %Laserfrequenz
        w0 = 2*pi*c/lamb0;  %Resonanzfrequenz
        
        %Zeit-Array
            t_array = 0:dt:tmax;
            N = length(t_array);
            N1 = round((tmax/2-tau)/dt);    %Index bei dem Pump-Puls auftritt
       
        E0 = EE.*gauss(t_array,tmax/2,pulseLength); %Feld des einfallenden Probe-Pulses


    %% Integration
        Sdot = [0,0,0];             %Init Bloch-Vektor
        M = [0, 0, (w0-wL)];        %M-Vektor
        dRho1 = dRho*dt/tau_relax;  %Hilfskonstante
        t0 = tmax/2-tau;         %Zeitpunkt Pump-Puls

        %Dgl Integration vor Pump-Puls
            for j=2:N1
                M(1) = -mu*E0(j-1)/hq;
                S(j,:) = S(j-1,:) + cross(M,S(j-1,:)).*dt;
            end

        %Pump-Puls anwenden
            S(j,3) = S(j,3)-dRho;
            if (S(j,3)<-1)
                S(j,3)=-1;
            end

        %Dgl Integration nach Pump-Puls
            for j=(N1+1):N
                M(1) = -mu*E0(j-1)/hq;
                S(j,:) = S(j-1,:) + cross(M,S(j-1,:)).*dt;

                S(j,3) = S(j,3) + exp(-(t_array(j)-t0)/tau_relax)*dRho1;
            end

    %% FFT
        rho01 = S(:,2); % =Es

        Es_w = fft(rho01./max(abs(rho01)));
        E0_w = fft(E0./max(E0).*exp(1i*dphi));
        
        Fs = 1/dt;                  % Sampling frequency
        f = Fs.*(0:dt:tmax)./tmax;  % Frequenz-array
        fshift = f-max(f)/2;
        
        dI(:,k) = fftshift(real(E0_w.*Es_w.'));  %Intensitätsänderung

    %% Vor-Plots
        %Frequenz-Achse in Energie-Achse wandeln (in meV)
            Eaxis = h.*fshift./e.*1e3;

%         %Probe-Puls
%             figure(1);
%             plot(t_array./1e-12,E0);
%             xlabel('t [ps]');
%             ylabel('Probe-Puls E0(t)');
% 
%         %rho01(t)
%             figure(2);
%             plot(t_array./1e-12,rho01);
%             xlabel('t [ps]');
%             ylabel('\rho_{01}(t) = E_S(t)');
% 
%         %dI(w) (detuning bei Verzögerung tau)
%             figure(3);
%             plot(Eaxis,dI(:,k)./max(dI(:,k)));
%             xlabel('detuning [meV]');
%             ylabel('dI');
%             xlim([-10,10]);
%             
%         %Bloch Vektor
%             figure(4);
%             hold on;
% %             plot(t_array./1e-12,S(:,1));
% %             plot(t_array./1e-12,S(:,2));
%             plot(t_array./1e-12,S(:,3));
%             hold off;
%             xlabel('t [ps]');
%             ylabel('Bloch-Vektor');
%             %legend('u','v','w');
            
end

%% Plot
    figure(5);
    surf(tau_array./1e-12, Eaxis, dI./max(max(dI)),'EdgeColor', 'None');%, 'facecolor', 'interp');
    view(2);
    ylim([-1,1]);
    xlabel('pump-probe-delay [ps]');
    ylabel('detuning [meV]');
    colormap Hot ;
    
