
function [ext, sca] =  mie(IoR, wavelength, nenv, radius)
% calculates Mie scattering and extinction cross sections
%
% input    IoR: column vector of index of refraction
%          wavelength: column vector, in nm
%           nenv: scalar index of refraction of environment
%          radius: scalar radius of particle, in nm
%
% output     ext, sca: cross section in nm^2 at wavelength 

    m = (IoR ./ nenv);
    x = ( 2 .* pi .* radius .* nenv ./ wavelength);

    % order to which evaluate MIE series
    nmax = 10;

    [a, b] = mie_ab (x,m,nmax);


    n = (1:nmax);
    nhelp = repmat(2 .* n +1 , length(x),1);
    
    k = 2 .* pi .* nenv ./ wavelength;

    ext =  (2 .* pi)./ k.^2 .* sum(nhelp .* real( a + b), 2);
    sca =  (2 .* pi)./ k.^2 .* sum(nhelp .* (abs(a) .^2 + abs(b) .^2),2);        

%     % error correction .. isNAN -> zero
%     a(isnan(a)) = 0;
%     b(isnan(b)) = 0;


end

