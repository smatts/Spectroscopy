function formula = make_time_domain_formula(track)
% build text string to calculate the NL polarisation in time domain
% corresponding to this Feynman diagram

    s = size(track.states);
    formula = '1';

    % iterate over all interactions
    for id = 2:s(1)
        
        column = track.interaction.column(id);
        
        % each interaction from the right makes a -1
        if (column == 1) 
            prefactor = '';
        else
            prefactor = ' .* -1 ';
        end

        oldstate = track.states(id-1, column);
        newstate = track.states(id, column);
        
        % dipole moment
        % we assume mu_ab = mu_ba and use only 'a<b' value
        dipole_states = sort ([oldstate, newstate]);
        dipole_term = sprintf('.* mu(%d,%d)',dipole_states);
        
        % phase
        leftstate = track.states(id, 1);
        rightstate = track.states(id, 2);
  
        % we use only omega_ab with 'a<b' and assume omega_ab = - omega_ba
        phase_states = sort ([leftstate, rightstate]);
        if (rightstate > leftstate)
            phase_sign = '+';
        else
            phase_sign = '-';
        end
        
        % rotating wave
        rotwave_term = '';
        if (id == 2)  rotwave_term = ' - omega21 '; end
        if (id == 4)  rotwave_term = ' - omega43 '; end
   
        % phase term
        if (leftstate ~= rightstate) 
            phase_term = sprintf('.* exp(%s1i .* (omega(%d,%d) %s ).* t%d)', phase_sign, phase_states, rotwave_term, id-1);
        else
            phase_term = '';
        end
        
        % decay term (the last state does not decay)
        if (id < s(1))
            decay_term =  sprintf('.* exp(-1 .* gamma(%d,%d) .* t%d)', phase_states, id-1);
        else
            decay_term = '';
        end
        
        formula = [formula, prefactor, dipole_term, phase_term, decay_term];
        
    end
    
end
