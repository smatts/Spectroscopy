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
    xhelp = repmat(x,1,length(n));

    [zj, zh] =  zjh(n,x);
    [zjm1, zhm1] =  zjh(n - 1 ,x);
    [zjp1, zhp1] =  zjh(n + 1 ,x);

    zjprime = sqrt(xhelp) .* ( (zjm1 - zjp1 + zj ./ xhelp) ./ 2);
    zhprime = sqrt(xhelp) .* ( (zhm1 - zhp1 + zh ./ xhelp) ./ 2);

    zj = sqrt(xhelp) .* zj;
    zh = sqrt(xhelp) .* zh;
end

function [zj, zh] = zjh(n, x) 

    xhelp = repmat(x,1,length(n));
    nhelp = repmat(n,length(x),1);

    % this is a Bessel J function, not a Bessel j function !!
    % See Wolfram research web page
    zj = besselj( nhelp + 0.5, xhelp);
    zh = besselj( nhelp + 0.5, xhelp) + 1i .* bessely( nhelp + 0.5, xhelp);
end

