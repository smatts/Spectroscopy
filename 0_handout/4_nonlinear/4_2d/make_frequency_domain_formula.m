function formula = make_frequency_domain_formula(track)
% build text string to calculate the NL polarisation in frequency domain
% corresponding to this Feynman diagram

    % these intervals we take in freq. domain
    freqdomain = [2,4];

    s = size(track.states);
    formula = '1';

    % iterate over interactions
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
        
        % laser
        if (track.phasematch(id) == 1)
            freq_sign = '+';
        else
            freq_sign = '-';
        end
        
        
        % freq. domain terms
        % we assume all omega_laser identical, and alternating up and down 
        % so that we can skip the sum over omega_laser
        if (ismember(id, freqdomain))
            freq_term = sprintf('.* 1i ./(%s1 .* omega%d + %s1 .* omega(%d,%d) + 1i .* gamma(%d,%d) ) ',freq_sign, id-1, phase_sign, phase_states, phase_states );
        else
             freq_term = '';
        end
     
        % time domain terms
        time_term = '';
        if not(ismember(id, freqdomain))
            if (id < s(1))  % last state does not decay
                decay_term =  sprintf('.* exp(-1 .* gamma(%d,%d) .* t%d)', phase_states, id-1);
            else
                decay_term = '';
            end
                      
            if (leftstate ~= rightstate) 
                time_term = sprintf('.* exp(%s1i .* (omega(%d,%d)  ).* t%d)',phase_sign, phase_states, id-1);
            end
             
            time_term = [time_term, decay_term];
      
        end
        
         
        formula = [formula, prefactor, dipole_term, time_term, freq_term];
        
    end
    
end
