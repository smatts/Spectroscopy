clear all ;

% integration time
Ti = [0 2^12] ;
TUNIT = 1e-13 ; 
ps = 1e-12 ./ TUNIT ; 

%sweeppar = 50:10*2^13/1e4:2300 ;

% provide a vector that represents the change of a parameter
% e.g., pump-arrival time
sweeppar = (20:1:270)*ps ;


fname = './results/sol_' ; 

thetime = datestr(now,30) ;
mkdir(strcat(fname,thetime)); 

for n = 1:1:length(sweeppar)
    
    m = mod(n,10) ;
    
    if m == 0 
        m = 10 ; 
    end    
    
    % properties of the system
    P{m}.dw    = 0 * TUNIT * 1e-6 / 6.57e-16 ;      % overall detuning
    P{m}.w12   = 0 * TUNIT * 1e-6 / 6.57e-16 ;     % splitting frequency E/hbar 
    P{m}.g     = [ log(2)/ (230*ps) log(2)/ (90*ps) 0 ; 
                   log(2)/ (230*ps) log(2)/ (90*ps) 0  ] ; % damping of the transitions, g(3) is the cross damping ,1/s
        
    P{m}.wpp       = [0 0] * TUNIT * 1e-6 / 6.57e-16 ;           % detuning of pump and probe
    P{m}.widthpp   = [2*ps 2*ps] / (2*sqrt(2*log(2))) ;                                   % pulse duration of pump and probe in units of 0.1*pico seconds
    P{m}.tpp       = [50*ps 150*ps] ;                        % delay of pump and probe in units of 0.1*pico seconds
               
    % normalized jones vector of the exciton transitions representing selection rules              
    % you can put imaginary numbers here, e.g., [1 i] for circular
    % selection rules, [1 0] for linear selection rules. Remember
    % normalization. 
    P{m}.J_ex1  = [1 ; 0] ;                   
    P{m}.J_ex2  = [0 ; 1] ;
    
    % normalized jones vectors of light field
    % you can put complex numbers here, e.g., sqrt(2)^-1 * [1 i] for circularly
    % polarized light
    P{m}.j_pump    = [ 1 ; 0 ] ;                           % normalized jones vector of pump pulse
    P{m}.j_probe   = [ 0 ; 1 ] ;                             % normalized jones vector of probe pulse
        
    % pulse area of the light field in [x y] direction, units of pi
    P{m}.e_pump    = sweeppar(n)*[1 ; 0] ;                            % pulse area of the pump  
    P{m}.e_probe   = [0 ; 1/2] ;                            % pulse area of the probe  
    
    % dember field paramters 
    P{m}.dember_par.rho            = log(2)/(50*ps) ;                % duration
    P{m}.dember_par.duration       = P{m}.widthpp(1) ;
    P{m}.dember_par.t0             = P{m}.tpp(1) ; 
    P{m}.dember_par.ys             = 250 ; 
    P{m}.dember_par.a              = 150 * sqrt(4*pi*P{m}.dember_par.duration^2) / (2*pi)^2  ;  
    
    % to fit the data, set ys = 243.5 kV/cm and a = 134.5 kV/cm / (2pi
    % rad)^2 * sqrt(4*pi*sigma^2)
    
    A_pump = 1 / sqrt( 2 * pi * P{m}.widthpp(1)^2) ; % integral of the pump   
    
    pump_pulse = @(t) sqrt(P{m}.e_pump' * P{m}.e_pump) * (pi * A_pump) ...
        * exp( - (t - P{m}.tpp(1))^2 / (2 * P{m}.widthpp(1)^2) ) ; 
    
    %% prepare dember field
    fun = @(t,y) ddember(t, y, P{m}.dember_par, pump_pulse) ;    
    dember_options = odeset('AbsTol',1e-12,'RelTol',1e-12,'InitialStep',1e-12,'Refine',1000,'Stats','off');    
    P{m}.dember_par.sol = ode113(fun,[P{m}.dember_par.t0-5*P{m}.dember_par.duration Ti(2)],0,dember_options) ;
    
    plot(P{m}.dember_par.sol.x, P{m}.dember_par.sol.y) ; 
    pause(0.1) ;
    
    %% solve optical bloch equations

    options = odeset('AbsTol',1e-12,'RelTol',1e-12,'InitialStep',1e-12,'Refine',100,'Stats','off');           
  
    % the solver increases the time stepsize fast when the system does not
    % change, as it is before the light pulses arrive. It the misses the
    % arrival of the short pulses. Therefore the calculation time of the
    % solver starts just before the arrival of the first pulse. Each
    % solution has thus a different time intervall, so you have to match
    % them when evaluating.
    
    % get the arrival time of the first pulse 
    [ts,id] = min(P{m}.tpp) ;
    
    % as long as the pulse arrives far from time origin, cut the time
    % interval. Otherwise do nothing..
    if ts-10*P{m}.widthpp(id) > 0,
        ti = [ts - 5*P{m}.widthpp(id) Ti(2)] ;
    else
        ti = Ti ;
    end
    
    % solve with pump pulse
    fun = @(t,r) OBE3(t,r,P{m}) ;
    solP{m} = ode113(fun,ti,[1 0 0 0 0 0 0 0 0 ],options) ;

    % solve without pump pulse    
    % if probe pulse is far away from t == 0, split time intervall
    if n == 1,
        if P{m}.tpp(2)-10*P{m}.widthpp(2) > 0,
            ti = [P{m}.tpp(2)-5*P{m}.widthpp(2) Ti(2)] ;
        else
            ti = Ti ;
        end
        
        % temp save pump amplitude
        ep = P{m}.e_pump ;
    
        % pump amplitude zero
        P{m}.e_pump = 0 ;
    
        fun = @(t,r) OBE3(t,r,P{m}) ;
        sol{m} = ode113(fun,ti,[1 0 0 0 0 0 0 0 0 ],options) ;

        % get pump amp back for saving
        P{m}.e_pump = ep ;
        
        save(sprintf('%s%s/sol_probe.mat',fname,thetime),'sol','P','sweeppar','Ti') ;
    end
    
    disp(sprintf('sweeppar %i',sweeppar(n)))  
    
    % save solutions in many files, ten solutions each, to reduce memory
    % load and limit file size
    if m == 10,
        save(sprintf('%s%s/sol_%1.3i.mat',fname,thetime,n/10),'solP','P','sweeppar','Ti') ;
    end    
end 

if n < 10,
    save(sprintf('%s%s/sol_001.mat',fname,thetime),'solP','P','sweeppar','Ti') ;
end

disp(sprintf('sol_%s',thetime)) ;

clear all ;