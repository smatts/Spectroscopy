close all
clear all

%Constants       % SI
h = 6.626E-34;        % J sec
hq = h/(2*pi);        % j sec
eps0 = 8.854E-12;       %  A sec / V m = A^2 sec / W m
NA = 6.022E23;             % 1 / mol
c = 299792458;           % m / sec
q_electron = 1.602E-19;  % C = A sec

debye = 0.208E-10 .* q_electron;  % m A sec

%Input
mode = 2;   % given:
            %1  eps
            %2  sigma
            %3  mu


% InGaAs
     mu = 34 .* debye;
     w_eV = 1.15;
     lambda = (h .* c ./ q_electron) / w_eV;
     omega = 2 .* pi .* c ./ lambda;

     d_omega_eV = 60e-3; % inhom. linewidth 0 60 meV
     d_omega = d_omega_eV .*  2 .* pi .* q_electron ./ h;
     
     % frequency integrated absorption cross section (m^2 / sec)
     sigma = pi .* omega .* mu.^2 ./ (1 .* hq .* c .* eps0);
     
     % approx peak sigma (m^2)
     sigma_peak = sigma ./ d_omega;
     sigma_peak .* 1e18;  % nm^2
     
     eps = NA * sigma_peak ./ log(10)  % m^2 /  Mol = 0.01 m^3 / (cm Mol)
                                        % = 0.01 1000 L / (cm Mol) 
                                        % = 10 / (cm M)
                      
     eps ./ 10;   % cm^-1 M^-1
     

% InGaAs V2

    alpha = 30 .* 100; % 30 / cm
    fwhm =  370E-9;  % m
    c_area = 3 .* 2e10 .* (100.^2);  % 2e10 / cm^2  % Teilchenzahldichte !!
                                % 3 layer of QD
    
    c_volume = c_area ./ fwhm;
    
    eps = alpha .* NA ./ (log(10) .* c_volume)
     
    lhs = eps ./ omega .* d_omega;
    mu = sqrt(lhs .* 1 .* log(10) .* hq .* c .* eps0 ./ ( pi .* NA));
    
    mu ./ debye
     
% CdSe

    eps = 2E5  .* 10 ; % 1/cm M  .* 10 = SI
    lambda = 550e-9;
    omega = 2 .* pi .* c ./ lambda;

    
    sigma_peak = eps .* log(10) ./ NA;
    sigma_peak .* 1e18  % nm^2
    
    
    
    mu_peak = sqrt( sigma_peak .* 3 .* hq .* c .* eps0 ./ (pi .* omega));
    mu_peak ./ debye
    
    d_lambda = 20e-9;
    d_omega = omega .* d_lambda ./ lambda;
    
    lhs = eps ./ omega .* d_omega;
    mu = sqrt(lhs .* 1 .* log(10) .* hq .* c .* eps0 ./ ( pi .* NA));
     mu ./ debye



% Xanthen dye
%     sigma = 0.44e-20;
%     lambda = 565e-9;

% MEH
%     sigma = 11.2e-18;
%     lambda = 815e-9%488e-9;
% 
% 
% %Frequ
% w = 2*pi*c/lambda;
% 
% 
% switch(mode)
%     case 1  %1  eps
%         sigma = log(10)*eps/NA;
%         mu = sqrt(3*hq*c*eps0*sigma/(pi*w));
%         
%     case 2  %2  sigma
%         eps = NA*sigma/log(10);
%         mu = sqrt(3*hq*c*eps0*sigma/(pi*w));
%         
%     case 3  %3  mu
%         sigma = pi*w*mu*mu/(3*hq*c*eps0);
%         eps = NA*sigma/log(10);
% end
% 
% 
% %Units
% eps = eps*1e-2;         %[cm^-1]
% sigma = sigma*1e18;     %[nm^2]
% mu = mu/(0.208e-10*q_electron);  %[D]
% 
