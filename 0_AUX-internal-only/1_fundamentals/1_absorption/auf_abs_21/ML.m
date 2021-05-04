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

% 
% % InGaAs
%      mu = 34 .* debye;
%      w_eV = 1.15;
%      lambda = (h .* c ./ q_electron) / w_eV;
%      omega = 2 .* pi .* c ./ lambda;
%          
%      tau = 1e-9; % sec
%      gamma = 2 .* pi / tau; % SI
%      
%      % peak absorption cross section, Eq. 1.30 (m^2 / sec)
%      sigma_peak = pi .* omega .* mu.^2 ./ (gamma .* hq .* c .* eps0);
%      sigma_peak .* 1e18  % nm^2  = 4.0053e+04 nm^2 = (200 nm)^2
%      
%      eps_peak = NA * sigma_peak ./ log(10)  % m^2 /  Mol = 0.01 m^3 / (cm Mol)
%                                         % = 0.01 1000 L / (cm Mol) 
%                                         % = 10 / (cm M)
%                       
%      eps_peak .* 10   % cm^-1 M^-1  = 1e11 /(cm M)

     
% 
% % InGaAs V2
% 
%     alpha = 30 .* 100; % 30 / cm
%     fwhm =  370E-9;  % m
%     c_area = 3 .* 2e10 .* (100.^2);  % 2e10 / cm^2  % Teilchenzahldichte !!
%                                 % 3 layer of QD
%     
%     c_volume = c_area ./ fwhm;
%     
%     eps = alpha .* NA ./ (log(10) .* c_volume) 
%         %  = 4.83e5 m^2 / Mol =4.8e4 / (cm M)
%    
%    %inhom. linewidth
%     d_omega_eV = 60e-3; % inhom. linewidth = 60 meV
%     d_omega = d_omega_eV .*  2 .* pi .* q_electron ./ h;
% 
%     % eq 1.33
%     lhs = eps ./ omega .* d_omega;
%     mu = sqrt(lhs .* 1 .* log(10) .* hq .* c .* eps0 ./ ( pi .* NA));
%     
%     mu ./ debye  % 28 D
%     
%     
%      
% CdSe

%     eps = 3E5  ./ 10 ; % 1/cm M  ./ 10 = SI
%     lambda = 550e-9;
%     omega = 2 .* pi .* c ./ lambda;
% 
%     
%     sigma_peak = eps .* log(10) ./ NA;
%     sigma_peak .* 1e18  % nm^2 = 7.5 nm^2 approx (3 nm)^2
%     
%     mu_peak = sqrt( sigma_peak .* 3 .* hq .* c .* eps0 ./ (pi .* omega));
%     mu_peak ./ debye  % = 7.3e-6 D
%     
%     d_lambda = 20e-9;
%     d_omega = omega .* d_lambda ./ lambda;
%     
%     lhs = eps ./ omega .* d_omega;
%     mu = sqrt(lhs .* 1 .* log(10) .* hq .* c .* eps0 ./ ( pi .* NA));
%      mu ./ debye   % = 47 D
% % 
% 
% % 
% % % %Xanthen dye
%     sigma_peak = 0.44e-16 .* 1e-4; % SI, m^2 = 0.44 A^2
%     
%     lambda = 565e-9; % Si m
%     omega = 2 .* pi .* c ./ lambda;
% 
%     eps = NA * sigma_peak / log(10); % SI
%     eps * 10  % =  115 /cm M %
% 
%      mu_peak = sqrt( sigma_peak .* 3 .* hq .* c .* eps0 ./ (pi .* omega));
%     mu_peak ./ debye  % = 1.7-e7D D
%     
%     d_lambda = 20e-9;
%     d_omega = omega .* d_lambda ./ lambda;
%     
%     lhs = eps ./ omega .* d_omega;
%     mu = sqrt(lhs .* 1 .* log(10) .* hq .* c .* eps0 ./ ( pi .* NA));
%      mu ./ debye   % = 1 D   
    

%MEH
   sigma_peak = 11.2e-18;
   lambda = 815e-9 ;%488e-9;
    omega = 2 .* pi .* c ./ lambda;

    eps = NA * sigma_peak / log(10); % SI
    eps * 10  % =  115 /cm M %

     mu_peak = sqrt( sigma_peak .* 3 .* hq .* c .* eps0 ./ (pi .* omega));
    mu_peak ./ debye  % = 1.7-e7D D
    
    d_lambda = 20e-9;
    d_omega = omega .* d_lambda ./ lambda;
    
    lhs = eps ./ omega .* d_omega;
    mu = sqrt(lhs .* 1 .* log(10) .* hq .* c .* eps0 ./ ( pi .* NA));
     mu ./ debye   % = 1 D   


