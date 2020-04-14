clear all ;
close all


                Ns = 1900:2150 ; 
                
thetime = datestr(now,30) ;

solname = './results/dember/sol_20160715T163348' ;
fouriername = 'fourier_20160616T144748.mat' ;
sweepparname = 'delay (ps)' ;

sweepparoffset = 0 ;
sweepparfactor = 1 ;

bd = 0 ; %burial depth in nm

numsol = 1:16 ; %1:10 ; % indices of files to include
nsp = 1:10 ;   % indices of soultions to include from each file

loadfourier = 0;    % load fourier data from previous runs
diagnostics = 0;    % diagnostic mode (load only one solution)
diag = [6 10] ;     % select a single solution for diagnostic mode [file, sol_idx]
plot_1D = 0;
subplot_switch = 0 ;
s = (1e-1 / 2.35) ; % spectrometer resolution

% calculate spectra only specific parameters
if diagnostics == 1,
   numsol = diag(1) ;
   nsp = diag(2) ;
end

% load probe solution (pump off) 
probe = load(sprintf('%s/sol_probe.mat',solname)) ;

flist = ls(sprintf('%s/sol_*.mat',solname)) ;

% time vector
TMAX = 2^12 ;
TNUM = 2^14 ; 
TUNIT = 1e-13 ; % time unit (0.1 ps)
ps =  10 ;

t  = linspace(0,TMAX,TNUM) ;   

% fourier transform frequency limit
fs = 1 / (TUNIT * TMAX / TNUM) ;  

% phase due to QD burial depth
% pha = 2 * pi * Ngaas * 2 * depth / la
pha = 2 * pi * 3.45 * 2 * bd / 720 ;

% rp contains rho(t) of the probe-only solution. the solution has a time 
% offset and we must add 'rho(t<t_offset) = init-state' by hand. ATM the
% initial state is the groundstate. 
%
% The format of rp is:
% rp = transp([rho11 rho12 rho13 rho21 rho22 rho23 rho31 rho32 rho33]) 

rpt = deval(probe.sol{1},t(t>probe.sol{1}.x(1))) ; 
rp  = [repmat(transpose([1 0 0 0 0 0 0 0 0]),1,length(t(t<=probe.sol{1}.x(1)))) rpt] ;    

fr_pumped = [] ;
fr_unpumped = [] ;
fr_probe = [] ; 
sweeppar = [] ;

if loadfourier==0,
    for k = numsol,
        %% load solution
        fname = sprintf('%s/%s',solname,flist(k,:)) ;
        pump = load(fname) ;

        %% fourier transform 
        for n = nsp,    
            % rpp contains rho(t) of the pump-probe solution. the solution has a time 
            % offset and we must add 'rho(t<t_offset) = init-state' by hand. ATM the
            % initial state is the groundstate. 
            %
            % The format of rpp is:
            % rpp = transp([rho11 rho12 rho13 rho21 rho22 rho23 rho31 rho32 rho33]) 
            rppt = deval(pump.solP{n},t(t>pump.solP{n}.x(1))) ; 
            rpp  = [repmat(transpose([1 0 0 0 0 0 0 0 0]),1,length(t(t<=pump.solP{n}.x(1)))) rppt] ;

            % get the off-diagonal elements
            % rho_xy(i,:) represents polarization of exciton i            
            rho_unpumped = imag(rp(2:3,:)) ;
            rho_pumped = imag(rpp(2:3,:)) ; 
            
            
            
            % define the probe electric field Jones vector
            Es = exp(-(t - pump.P{n}.tpp(2)).^2 / (2 * pump.P{n}.widthpp(2)^2) ) ... 
                .* exp( -1i * ( pump.P{n}.wpp(2) ) * t) ; 
            
            [Es, Js] =  meshgrid(Es, pump.P{n}.j_probe) ;  
            E_pulse = Js .* fftshift(fft(Es,[],2)) ;        
    
            fr_probe = [fr_probe ; E_pulse] ;
            
            % the next section circshifts in a very confusing and complicated
            % way the rho vectors such that the fft starts just when the
            % first incoming pulse is over. I have forgotten why it is so
            % complicated and confusing and keep it as it is.
            
            jex1 = pump.P{n}.J_ex1 ;
            jex2 = pump.P{n}.J_ex2 ;
            jpump = pump.P{n}.j_pump ;
            jprobe = pump.P{n}.j_probe ;
            
            % calculate fourier transforms
            % the format of fr_xy is:
            %
            % [ fft(rho12_sol1) ; 
            %   fft(rho13_sol1) ;
            %   fft(rho12_sol2) ; 
            %   fft(rho13_sol2) ;
            %      ...
            %   fft(rho12_solN) ; 
            %   fft(rho13_solN) ; ]
            
            fr_unpumped = [fr_unpumped ; fftshift(fft(rho_unpumped,[],2)/length(t),2)] ;
            fr_pumped   = [fr_pumped ; fftshift(fft(rho_pumped,[],2)/length(t),2)] ;

            % sweepparameter vector
            sweeppar = [sweeppar pump.sweeppar(n+10*(k-1))] ;

            % save circshift length
%             idpp1(n+10*(k-1)) = id_peak_1 ;
%             idpp2(n+10*(k-1)) = id_peak_2 ;

        end   

        disp(sprintf('file %i',k)) 
        
        %%
        for o=1  % dummy for-loop to hide diagnostic routines in matlab editor
            if diagnostics,
                
                
                xl = [0 700] ; 
                

                figure ; 
                
                for n = [1, 5, 9],
                    if subplot_switch ==1,
                        subplot(3,3,n) ;
                    else
                        figure 
                    end ; 
                    plot(t/ps,rpp(n,:),'r','linewidth',3) ;
                    xlim(xl); 
                    ylim([-0.05 1.05]) ; 
                    hold on ;
                    xlim(xl); 
                    ylim([0 1]) ; 
                end
                
                for n = setxor(1:9, [1,5,9]),
                    if subplot_switch ==1,
                        subplot(3,3,n) ;
                    else
                        figure
                    end
                    plot(t/ps,real(rpp(n,:)),'r','linewidth',3) ;
                    xlim(xl);
                    ylim([-0.5 0.5]) ; 
                    hold on ;
                    plot(t/ps,real(rp(n,:)),'g','linewidth',3) ;
                    xlim(xl); 
                    ylim([-0.5 0.5]) ; 
                end
                
                figure ;
                title('imaginary parts') ;
                for n = [1, 5, 9],
                    if subplot_switch ==1,
                        subplot(3,3,n) ;
                    else
                        figure 
                    end ; 
                    plot(t/ps,rp(n,:),'r','linewidth',3) ;
                    xlim(xl); 
                    ylim([-0.05 1.05]) ; 
                    hold on ;
                    xlim(xl); 
                    ylim([0 1]) ; 
                end
                
                for n = setxor(1:9, [1,5,9]),
                    if subplot_switch ==1,
                        subplot(3,3,n) ;
                    else
                        figure
                    end
                    plot(t/ps,imag(rpp(n,:)),'r','linewidth',3) ;
                    xlim(xl);
                    ylim([-0.5 0.5]) ; 
                    hold on ;
                    plot(t/ps,imag(rp(n,:)),'g','linewidth',3) ;
                    xlim(xl); 
                    ylim([-0.5 0.5]) ; 
                end

                figure 
                plot(real(fftshift(fft(real((rpp(2,:) + rpp(3,:)) - (rp(2,:) + rp(3,:)))))))
                %daspect([1 1 1])
                
                %figure ; 
                %hold on ;
                %plot(abs(transpose((fr_pumped)))) ;
                %plot(real(transpose(fr_unpumped)),'g') ; 
                %plot(real(transpose(fr_pumped - fr_unpumped)),'r') ;
                %xlim([4600 5400]) ;
        %%
                u1 =  2*real(rpp(2,:)) ;
                v1 = -2*imag(rpp(2,:)) ;
                w1 =  rpp(5,:) - (rpp(1,:) + rpp(9,:));

                u2 =  2*real(rpp(3,:)) ;
                v2 = -2*imag(rpp(3,:)) ;
                w2 =  rpp(9,:) - (rpp(1,:) +  rpp(5,:));

                u3 =  2*real(rpp(6,:)) ;
                v3 = -2*imag(rpp(6,:)) ;
                w3 =  rpp(9,:) - (rpp(5,:) +  rpp(1,:));

                u1p =  2*real(rp(2,:)) ;
                v1p = -2*imag(rp(2,:)) ;
                w1p =  rp(5,:) - (rp(1,:)  +  rp(9,:));

                u2p =  2*real(rp(3,:)) ;
                v2p = -2*imag(rp(3,:)) ;
                w2p =  rp(9,:) - (rp(1,:)  +  rp(5,:));

                u3p =  2*real(rp(6,:)) ;
                v3p = -2*imag(rp(6,:)) ;
                w3p =  rp(9,:) - (rp(5,:) +  rp(1,:));

%%
                if subplot_switch == 1,       
                    figure
                    subplot(1,2,1) ;                
                else 
                    figure ; 
                end 
                               
                grid off ;
                axis off ;                 
                xlim([-1 1]) 
                ylim([-1 1]) 
                zlim([-1 1]) 
                view(10,6) ;
                daspect([1 1 1]) ;
                hold on ;            
                plot3(u1(Ns),v1(Ns),w1(Ns),'o-','markersize', 7, 'markerfacecolor',[0 0 1], 'markeredgecolor','none') ;
                plot3(u1p(Ns),v1p(Ns),w1p(Ns),'or-','markersize', 7, 'markerfacecolor',[1 0 0], 'markeredgecolor','none') ;
                alpha 0.5 ;
                phi = 0:0.01:2*pi ;                 
                plot3(sin(phi), cos(phi), zeros(1,length(phi)),'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
                plot3(cos(phi), zeros(1,length(phi)), sin(phi),'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
                plot3(zeros(1,length(phi)), sin(phi),cos(phi), 'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
                title('S_{01}') ;
                xlabel('u') ;
                ylabel('v') ;
                zlabel('w') ;   
                
                if subplot_switch == 1 
                    subplot(1,2,2) ;                
                else
                    figure ; 
                end 
                
                grid off ;
                axis off ;                 
                xlim([-1 1]) 
                ylim([-1 1]) 
                zlim([-1 1]) 
                view(45,40) ;
                daspect([1 1 1]) ;
                hold on ;            
                plot3(u2,v2,w2,'o','markersize', 7, 'markerfacecolor',[0 0 1], 'markeredgecolor','none') ;
                plot3(u2p,v2p,w2p,'or','markersize', 7, 'markerfacecolor',[1 0 0], 'markeredgecolor','none') ;
                alpha 0.5 ;
                phi = 0:0.01:2*pi ;                 
                plot3(sin(phi), cos(phi), zeros(1,length(phi)),'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
                plot3(cos(phi), zeros(1,length(phi)), sin(phi),'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
                plot3(zeros(1,length(phi)), sin(phi),cos(phi), 'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
                title('S_{02}') ;
                xlabel('u') ;
                ylabel('v') ;
                zlabel('w') ;

%                 if subplot_switch == 1,
%                     subplot(1,3,3) ;                    
%                 else
%                   figure ; 
%                 end
%                 
%                 grid off ;
%                 axis off ;                 
%                 xlim([-1 1]) 
%                 ylim([-1 1]) 
%                 zlim([-1 1]) 
%                 view(45,40) ;
%                 daspect([1 1 1]) ;
%                 hold on ;            
%                 plot3(u3,v3,w3,'o','markersize', 7, 'markerfacecolor',[0 0 1], 'markeredgecolor','none') ;
%                 alpha 0.5 ;
%                 phi = 0:0.01:2*pi ;                 
%                 plot3(sin(phi), cos(phi), zeros(1,length(phi)),'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
%                 plot3(cos(phi), zeros(1,length(phi)), sin(phi),'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
%                 plot3(zeros(1,length(phi)), sin(phi),cos(phi), 'color',[0.7 0.7 0.7], 'linewidth', 3) ; 
%                 title('S_{12}') ;
%                 xlabel('u') ;
%                 ylabel('v') ;
%                 zlabel('w') ;
%                 
                %%
                figure 
                plot(u1, v1,'o','markersize', 7, 'markerfacecolor',[0 0 1], 'markeredgecolor','none') ;
                hold on ;
                plot(sin(phi), cos(phi),'color',[0.7 0.7 0.7], 'linewidth', 3) ;
                xlim([-1.2, 1.2]) ;
                ylim([-1.2, 1.2]) ;
                daspect([1 1 1]) ;

                figure 
                plot(v2, w2,'o','markersize', 7, 'markerfacecolor',[0 0 1], 'markeredgecolor','none') ;
                hold on ;
                plot(sin(phi), cos(phi),'color',[0.7 0.7 0.7], 'linewidth', 3) ;
                xlim([-1.2, 1.2]) ;
                ylim([-1.2, 1.2]) ;
                daspect([1 1 1]) ;
                                
            end
        end 

        clear rpp rppt ;
    end % for numsol

clear rp rpt ; 
end

if loadfourier == 1,
    load(sprintf('%s/%s',solname,fouriername)) ;
end 

if diagnostics == 0,
    if loadfourier==0,
        save(sprintf('%s/fourier_%s.mat',solname,thetime),'fr_pumped','fr_unpumped','fr_probe','t','sweeppar') ;
    end

       
    % these are the projections of the exciton axes to the x and y axes
    % J_ex = [x-proj y-proj] 
    jex1 = probe.P{1}.J_ex1 ; 
    jex2 = probe.P{1}.J_ex2 ;
           
    refl_factor = 10^3; % the reflectivity of the surface is by this factor stronger than the QD emission
       
    m1 = max(max(abs(fr_pumped),[],2)) ;
    m2 = max(max(abs(fr_unpumped),[],2)) ;
    
    qd_peak = max([m1 m2]) ;
    
    
    % calculate the scattered electric fields along x- and y-axis as a sum
    % of weighted contributions by exciton 1 and exciton 2:
    %
    % E_x = fft(ex1) * x-projection + fft(ex2) * x-projection
    % E_y = fft(ex1) * y-projection + fft(ex2) * y-projection
    %
    % for later use the 2D array is reshaped into a 1D-vector, where the
    % vectors of individual solutions are arranged in a sequence
    
    E_pumped_j(1,:) = reshape(transpose(fr_pumped(1:2:end,:)),1,[]) * jex1(1) ... 
                       + reshape(transpose(fr_pumped(2:2:end,:)),1,[]) * jex2(1) ;
    E_pumped_j(2,:) = reshape(transpose(fr_pumped(1:2:end,:)),1,[]) * jex1(2) ...
                       + reshape(transpose(fr_pumped(2:2:end,:)),1,[]) * jex2(2) ;          
    
    E_pulse_j(1,:) = reshape(transpose(fr_probe(1:2:end,:)),1,[]) ;   
    E_pulse_j(2,:) = reshape(transpose(fr_probe(2:2:end,:)),1,[]) ;
                   
    clear fr_pumped ;    
    
    % the scattered electric field is added to the spectrally flat surface
    % reflection of the probe beam.
    %E_pumped_j = refl_factor * qd_peak * repmat(probe.P{1}.j_probe,1,size(E_pumped_j,2)) + exp(1i * pha) * E_pumped_j ;
    E_pumped_j = E_pulse_j * qd_peak + exp(1i * pha)/refl_factor * E_pumped_j ;
    
    % calculate intensities and reshape into 2D-array
    I_pumped = sum(conj(E_pumped_j).*E_pumped_j,1) ;
    I_pumped = transpose(reshape(I_pumped,length(t),[])) ;
      
    clear E_pumped_j ;    
    
    % repeat for unpumped solution
    E_unpumped_j(1,:) = reshape(transpose(fr_unpumped(1:2:end,:)),1,[]) * jex1(1) ... 
                       + reshape(transpose(fr_unpumped(2:2:end,:)),1,[]) * jex2(1) ;
    E_unpumped_j(2,:) = reshape(transpose(fr_unpumped(1:2:end,:)),1,[]) * jex1(2) ...
                       + reshape(transpose(fr_unpumped(2:2:end,:)),1,[]) * jex2(2) ; 
                   
    clear fr_unpumped ;               
   
    %E_unpumped = refl_factor * qd_peak * repmat(probe.P{1}.j_probe,1,size(E_unpumped_j,2)) + exp(1i * pha) * E_unpumped_j ;
    E_unpumped = E_pulse_j* qd_peak + exp(1i * pha) / refl_factor * E_unpumped_j ;
    clear probe pump ;
    
    I_unpumped = sum(conj(E_unpumped).*E_unpumped,1) ;
    I_unpumped = transpose(reshape(I_unpumped,length(t),[])) ;
    
    clear E_unpumped_j ;    
    
    % frequency axis
    f = 1e3 * 4.13e-15 * fs/2 * linspace(-1,1,size(I_pumped,2)) ;

    
    % gaussian representing the spectrometer resolution
    ef2 = (max(f) - min(f))/length(f)/sqrt(2*pi*s^2)*exp( -f.^2 / (2 * s^2) ) ;
    
        
    
    % convolute difference spectra wih gaussian    
    for n = 1:size(I_pumped,1)
        dIoI = I_pumped(n,:)-I_unpumped(n,:) ; 
        b(n,:) = conv(dIoI, ef2, 'same');
    end

    sp = sweepparfactor*(sweeppar - sweepparoffset) ;

    
    
    frange = floor(length(f) * (1/2 - 0.025)) : ceil( length(f) * (1/2 + 0.025 )) ; 
    
    if plot_1D == 0,
        fig = figure ; 
        p = pcolor(sp,f(frange),b(:,frange)') ; shading flat ;
        ylabel('detuning (meV)') ;
        xlabel(sweepparname) ; 
        
        figure ; 
        plot(sp,max(b,[],2)) ;
    else
        fig = figure('Units','centimeters','position',[10 10 4.5 3.5]); 
        p = plot(f(frange),b(:,frange),'linewidth',3) ;
        set(gca,'fontname','Palatino Linotype','fontsize',9)
        ylim([-.5 .5]*1e-3) ;
        xlim([ceil(f(frange(1))) floor(f(frange(end)))])
        ylabel('dR/R') ;
        xlabel('detuning (meV)') ; 
    end


    
    saveas(fig,sprintf('%s/spec_%s.fig',solname,thetime)) ;
    save(sprintf('%s/out_%s.mat',solname,thetime),'sp','f','b') ;
end 










