clear all

sweeppar = linspace(0,8,40) ; %logspace(-2,1,10) ;
sweeppar2 = 1 ; % 2.4:0.02:2.5 ; 

Ti = [0 2^13] ;
TUNIT = 1e-13 ; 
ps = 1e-12 ./ TUNIT ; 

sol = cell(length(sweeppar2), length(sweeppar));
del = cell(length(sweeppar));
kap = cell(length(sweeppar));
gam = cell(length(sweeppar));
m = zeros(length(sweeppar)) ; 
id = zeros(length(sweeppar)) ;
mkap = (zeros(length(sweeppar2),length(sweeppar))) ; 
mdel = (zeros(length(sweeppar2),length(sweeppar))) ; 
mgam = (zeros(length(sweeppar2),length(sweeppar))) ; 

delta = @(F) - 0.23 * (F + 17.86).^2 - 74.6 ;

kappa = @(F) -7.233e-6 * (F - 138.8).^2 + 7.233e-6*138.8^2 + 1 ;
kappa2 = @(F) 1.37e-3 * F + 1 ; 
gamma = @(F) -0.46 * F ;  

% f1 = figure ;
% ax1 = axes ;
%hold on ;
% 
% f2 = figure ; 
% ax2 = axes ; 
% hold on ;

f3 = figure ;
ax3 = axes ;
hold on ;

for n = 1:length(sweeppar2),
    for k = 1:length(sweeppar),
        
        fun = @(t,y) ddember(t,y,dember_par,@(x) gaussian(x,pump_par)) ;

        options = odeset('AbsTol',1e-12,'RelTol',1e-12,'InitialStep',1e-12,'Refine',1000,'Stats','off');

        dember_par.rho            = log(2)/(50*ps) ;                % duration
        dember_par.duration       = 2*ps / (2 * sqrt(2 * log(2)));
        dember_par.t0             = 50*ps ; 
        dember_par.ys             = 125 ; 
        dember_par.a              = 100 * sqrt(4*pi*dember_par.duration^2) / (2*pi)^2 ; %sqrt(dember_par.duration^2 /(4*pi^3)) ;

        e_pump                    = sweeppar(k)*[1 ; 0] ;                            

        A_pump = 1 / sqrt( 2 * pi * dember_par.duration^2) ; % integral of the pump   

        pump_pulse = @(t) sqrt(e_pump' * e_pump) * (pi * A_pump) ...
            * exp( - (t - dember_par.t0).^2 / (2 * dember_par.duration^2) ) ; 

        %% prepare dember field
        fun = @(t,y) ddember(t, y, dember_par, pump_pulse) ;    
        dember_options = odeset('AbsTol',1e-12,'RelTol',1e-12,'InitialStep',1e-12,'Refine',1000,'Stats','off');    
        sol{n,k} = ode113(fun,[dember_par.t0-5*dember_par.duration Ti(2)],0,dember_options) ;
        
        
        del{n,k} = delta(-sol{n,k}.y) ; 
        kap{n,k} = kappa(-sol{n,k}.y) ; 
        gam{n,k} = gamma(-sol{n,k}.y) ; 
% % 
       %plot(ax1, sol{n,k}.x/10-pump_par.t0/10, sol{n,k}.y) ; hold on ;
% 
%         plot(ax2, sol{k}.x-pump_par.t0, del{k}) ; hold on ;
% 
        plot(ax3, sol{n,k}.x-dember_par.t0, sqrt(kap{n,k})) ; hold on ;




        [m(n,k), tmp] = max(sol{n,k}.y) ; 
        
        mkap(n,k) = kappa(-m(n,k))  ;
        mdel(n,k) = delta(-m(n,k))-delta(0) ;
        mgam(n,k) = max(abs(gam{n,k}))  ;
        
        id(n,k) = sol{n,k}.x(tmp) ;     
    end
    
    
end

%plot(ax1,sol{k}.x/10-pump_par.t0/10, 1/5*gaussian(sol{k}.x,pump_par),'g') ; 
% ylim([1e-3 5])
% xlim([-30, 70])

%plot(ax1,id-pump_par.t0,m,'s');

