%Constants
h = 6.626e-34;
hq = h/(2*pi);
eps0 = 8.854e-12;
Na = 6.022e23;
c = 299792458;
e = 1.602e-19;


%Input
mode = 2;   % given:
            %1  eps
            %2  sigma
            %3  mu


% InGaAs
%     mu = 31*(0.208e-10*e);
%     lambda = h*c/(1.15*e);

% CdSe
%     eps = 1.3e7;
%     lambda = 350e-9;

% Xanthen dye
%     sigma = 0.44e-20;
%     lambda = 565e-9;

% MEH
    sigma = 11.2e-18;
    lambda = 815e-9%488e-9;


%Frequ
w = 2*pi*c/lambda;


switch(mode)
    case 1  %1  eps
        sigma = log(10)*eps/Na;
        mu = sqrt(3*hq*c*eps0*sigma/(pi*w));
        
    case 2  %2  sigma
        eps = Na*sigma/log(10);
        mu = sqrt(3*hq*c*eps0*sigma/(pi*w));
        
    case 3  %3  mu
        sigma = pi*w*mu*mu/(3*hq*c*eps0);
        eps = Na*sigma/log(10);
end


%Units
eps = eps*1e-2;         %[cm^-1]
sigma = sigma*1e18;     %[nm^2]
mu = mu/(0.208e-10*e);  %[D]

%Output
disp(['Wavelength: ', num2str(lambda*1e9),' nm']);
disp(['Absorptionskoeffizient eps: ', num2str(eps),' cm^-1']);
disp(['Absorptionsquerschnitt sigma: ', num2str(sigma),' nm^2']);
disp(['Übergangsdipolmoment mu: ', num2str(mu),' D']);
disp(['typische Länge L: ', num2str(mu*0.208e-1),' nm']);