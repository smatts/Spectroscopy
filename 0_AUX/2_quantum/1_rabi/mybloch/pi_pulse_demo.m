clear all 
close all

% integration time
Ti = [0 2^12] ;  % = 0 ... 409.6ps
TUNIT = 1e-13 ; 
ps = 1e-12 ./ TUNIT ; 


% properties of the system
P.dw    = 0 * TUNIT * 1e-6 / 6.57e-16 ;      % overall detuning E/hbar
P.widthpp   = 35*ps  / (2*sqrt(2*log(2))) ;   % pulse duration of pump and probe in units of 0.1*pico seconds
P.tpp       = 50*ps ;                        % delay of pump and probe in units of 0.1*pico seconds
P.e_pump    = 0.5;    % pulse area of the pump  
P.wpp(1) = 0 * TUNIT * 1e-6 / 6.57e-16 ;

%% solve optical bloch equations


% the solver increases the time stepsize fast when the system does not
% change, as it is before the light pulses arrive. It the misses the
% arrival of the short pulses. Therefore the calculation time of the
% solver starts just before the arrival of the first pulse. Each
% solution has thus a different time intervall, so you have to match
% them when evaluating.

% get the arrival time of the first pulse 
[ts,id] = min(P.tpp) ;

% as long as the pulse arrives far from time origin, cut the time
% interval. Otherwise do nothing..
if ts-10*P.widthpp(id) > 0
    ti = [ts - 5*P.widthpp(id) Ti(2)] ;
else
    ti = Ti ;
end


% solve equations
options = odeset('AbsTol',1e-12,'RelTol',1e-12,'InitialStep',1e-12,'Refine',100,'Stats','off');
fun = @(t,r) OBE_TLS_density_matrix(t,r,P) ;

solP = ode113(fun,ti,[1 0 0 0],options) ;


% Plot as Bloch vector
u1 =  2*real(solP.y(2,:)) ;
v1 = -2*imag(solP.y(2,:)) ;
w1 =  solP.y(1,:) - solP.y(4,:) ;

tau = (min(solP.x): ps./3 :max(solP.x));

u = interp1(solP.x,u1, tau);
v = interp1(solP.x,v1, tau);
w = interp1(solP.x,w1, tau);
% 
% Ns=1:length(u1);
% 
% grid off ;
% axis off ;                 
% xlim([-1 1]) 
% ylim([-1 1]) 
% zlim([-1 1]) 
% view(10,6) ;
% daspect([1 1 1]) ;
% hold on ;            
% plot3(u,v,w,'o-','markersize', 7, 'markerfacecolor',[0 0 1], 'markeredgecolor','none') ;
% alpha 0.5 ;
% phi = 0:0.01:2*pi ;                 
% plot3(sin(phi), cos(phi), zeros(1,length(phi)),'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
% plot3(cos(phi), zeros(1,length(phi)), sin(phi),'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
% plot3(zeros(1,length(phi)), sin(phi),cos(phi), 'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
% title('S_{01}') ;
% xlabel('u') ;
% ylabel('v') ;
% zlabel('w') ; 
% 
% data = [u,v,w]

ids = find( abs(tau - P.tpp(1)) < (P.widthpp(1)));
field = tau .* 0;
field(ids) = 1;

    
plot(tau, w)
hold on
plot(tau, field)
%plot(tau, u)
plot(tau, -v)
hold off
xlim([200 900])

%%-------------
% FT

spec_pump = FFT(field);



function dr = OBE_TLS_density_matrix(t,r,P) 
% Optical Bloch equations for 2x2 density matrix
    
    % reshape rho into 2x2 matrix
    r = transpose(reshape(r,2,[])) ;          
    
    %% calculate Rabi frequency

    % Gauss
%      wR1 = P.e_pump ./ (sqrt( 2 / pi) .*  P.widthpp(1))  .* exp( - (t - P.tpp(1))^2 / (2 * P.widthpp(1)^2) )  ...
%            * exp( -1i * ( P.wpp(1) ) * t );
% %       
    % rect
    if ( abs(t - P.tpp(1)) < (P.widthpp(1)) )
        wR1 = pi .* P.e_pump ./( 2.* P.widthpp(1)) * exp( -1i * ( P.wpp(1) ) * t );
    else
        wR1 = 0;
    end
    
    %% Equations of motion for density matrix elements
     
    H =  [0,   wR1 ;  conj(wR1), P.dw ]  ;
        
    dr = -1i * (r*H - H*r) ;       
         
    % transform back into column vector as required by ODE45 solver 
    dr = reshape(transpose(dr),[],1) ; 
end

