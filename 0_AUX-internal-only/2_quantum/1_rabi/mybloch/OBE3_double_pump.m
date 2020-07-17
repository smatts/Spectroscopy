function dr = OBE3(t,r,P) 
    lr = L(P) * r ;
    
    % reshape rho into 3x3 matrix
    r = transpose(reshape(r,3,[])) ;    

    %% calculate rabi frequency
    
    a_pump = 1 / (2 * P.widthpp(1)^2) ; 
    a_probe = 1 / (2 * P.widthpp(2)^2) ; 
    
        
    % construct jones vector of pump and probe and scale according to pulse
    % area
    J_Pump = P.j_pump .* P.e_pump * sqrt( pi * a_pump / 4 ) ;
    J_Probe = P.j_probe .* P.e_probe * sqrt( pi * a_probe / 4 ) ;
    
    % create time dependent rabi frequency and project to axes of the exciton transitions    
    % exciton 1
    wR1 = (P.J_ex1' * J_Pump) .* (exp( - (t - P.tpp(1))^2 / (2 * P.widthpp(1)^2) ) ... 
                + 1i*exp( - (t - (P.tpp(2) + 200))^2 / (2 * P.widthpp(1)^2) )) * exp( -1i * ( P.wpp(1) ) * t ) ...                
        + (P.J_ex1' * J_Probe) .* exp( - (t - P.tpp(2))^2 / (2 * P.widthpp(2)^2) ) ... 
                * exp( -1i * ( P.wpp(2) ) * t ) ;
          
    % eciton 2
    wR2 = P.J_ex2' * J_Pump * (exp( - (t - P.tpp(1))^2 / (2 * P.widthpp(1)^2) ) ... 
                + 1i*exp( - (t - (P.tpp(2) + 200))^2 / (2 * P.widthpp(1)^2) )) * exp( -1i * ( P.wpp(1) ) * t ) ...                
        + P.J_ex2' * J_Probe * exp( - (t - P.tpp(2))^2 / (2 * P.widthpp(2)^2) ) ... 
                * exp( -1i * ( P.wpp(2) ) * t ) ;        
     
    
    %% Equations of motion for density matrix elements
    % hamiltonian 
    % taken from Ficek (PHYSICAL REVIEW A 69, 023401 2004)
    H = (P.dw - P.w12)*A(2,2) + P.dw*A(3,3)...
            + ( wR1*A(2,1) + wR2*A(3,1) )...
            + ( wR1*A(2,1) + wR2*A(3,1) )' ;
        
    % optical bloch equation, drho/dt = -i*[rho,H]+L*rho
    % where H is the hamiltonian as defined above and L is the damping
    % taken from Ficek (PHYSICAL REVIEW A 69, 023401 2004)
    dr = -1i*(r*H - H*r) ; 
         

         
    % transform back into column vector as required by ODE45 solver 
    dr = reshape(transpose(dr),[],1) + lr ; 
end
