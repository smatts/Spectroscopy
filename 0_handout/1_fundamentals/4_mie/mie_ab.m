function [mie_a, mie_b] = mie_ab (x, m, n_max)
% Mie a_n and b_n coefficients
%
% input   x: column vector of x = k a = reduced size \approx wavelength
%         m: column vector of m = n_particle / n_medium for each wavelength
%         n_nmax: maximum order  
%   
% output  a,b : matrix of size (length(x), n_max)

  n = (1:n_max);
    
  mhelp = repmat(m,1,n_max);

  [zj,  zh,  zjprime,  zhprime]  = riccati(n, x);
  [zjm, zhm, zjprimem, zhprimem] = riccati(n, m .* x);
  
  mie_a = (mhelp .* zjm .* zjprime - zj .* zjprimem) ./  ...
          (mhelp .* zjm .* zhprime - zh .* zjprimem);
       
  mie_b = (zjm .* zjprime - mhelp .* zj .* zjprimem) ./ ...
          (zjm .* zhprime - mhelp .* zh .* zjprimem);

end

%-----------
% helper functions

function [zj, zh, zjprime, zhprime] = riccati(n, x) 
% Ricatti-Bessel functions and their derivatives
% see Bohren / Huffman p.101

    xhelp = repmat(x,1,length(n));

    [zj, zh] =  spherical_bessel(n,x);
    [zjm1, zhm1] =  spherical_bessel(n - 1 ,x);
    [zjp1, zhp1] =  spherical_bessel(n + 1 ,x);

% for the derivatives see
% http://functions.wolfram.com/Bessel-TypeFunctions/SphericalBesselJ/20/ShowAll.html
    zjprime = xhelp .* ( (zjm1 - zjp1 + (zj ./ xhelp)) ./ 2);
    zhprime = xhelp .* ( (zhm1 - zhp1 + (zh ./ xhelp)) ./ 2);

    zj = xhelp .* zj;
    zh = xhelp .* zh;
end

function [zj, zh] = spherical_bessel(n, x) 
%  zj = spherical Bessel functions of fist kind
%  zh = spherical Bessel functions of third kind = Hankel function

    xhelp = repmat(x,1,length(n));
    nhelp = repmat(n,length(x),1);

    % we need _spherical_ Bessel functions j,h
    % but we have only Bessel functions J, Y
    % see  
    % Bohren / Huffman eq. 4.9
    %
    zj = sqrt(pi ./ (2 .* xhelp)) .* besselj( nhelp + 0.5, xhelp);
    zy = sqrt(pi ./ (2 .* xhelp)) .* bessely( nhelp + 0.5, xhelp);
    zh = zj + 1i .* zy;
    
    
end

