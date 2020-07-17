function dy = ddember(t, y, p, f)
   dy = -p.rho*y + p.a*f(t)^2*(p.ys - y)/p.ys ; 
end