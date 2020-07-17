function dr = OBE_TLS_density_matrix(t,r,P) 
% Optical Bloch equations for 2x2 density matrix
    
    % reshape rho into 2x2 matrix
    r = transpose(reshape(r,2,[])) ;          
    
    %% calculate Rabi frequency

    wR1 = P.e_pump ./ (sqrt( 2 / pi) .*  P.widthpp(1))  .* exp( - (t - P.tpp(1))^2 / (2 * P.widthpp(1)^2) )  ...
          * exp( -1i * ( P.wpp(1) ) * t );
      
    %% Equations of motion for density matrix elements
     
    H =  [0,   wR1 ;  conj(wR1), P.dw ]  ;
        
    dr = -1i * (r*H - H*r) ;       
         
    % transform back into column vector as required by ODE45 solver 
    dr = reshape(transpose(dr),[],1) ; 
end

