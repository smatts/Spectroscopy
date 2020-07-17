function dr = OBE3_ML(t,r,P) 
    

    
%    P.g(1,2) = P.g(1,2);
    
    %% calculate damping
   % lr = L(P) * r ;
    
    % reshape rho into 2x2 matrix
    r = transpose(reshape(r,2,[])) ;          
    
    %% calculate rabi frequency
    
    
    % create time dependent rabi frequency and project to axes of the exciton transitions    
    % exciton 1
    wR1 = P.e_pump ./ (sqrt( 2 / pi) .*  P.widthpp(1))  .* exp( - (t - P.tpp(1))^2 / (2 * P.widthpp(1)^2) )  ;
          
      
    %% Equations of motion for density matrix elements
    % hamiltonian 
    % taken from Ficek (PHYSICAL REVIEW A 69, 023401 2004)
    
    
    H =   P.dw *element(2,2) ...
            +   wR1*element(2,1) + ( wR1*element(2,1))' ;
        
         
    % optical bloch equation, drho/dt = -i*[rho,H]+L*rho
    % where H is the hamiltonian as defined above and L is the damping
    % taken from Ficek (PHYSICAL REVIEW A 69, 023401 2004)
    dr = -1i*(r*H - H*r) ;       
         
    % transform back into column vector as required by ODE45 solver 
    dr = reshape(transpose(dr),[],1) ; % + lr ; 
end


function [M] = element(i,j)
    M = zeros(2,2) ;
    M(i,j) = 1 ;
end