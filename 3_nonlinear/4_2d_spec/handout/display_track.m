function display_track(tracks)
% displays tracks as double-sided Feynman diagrams
% all tracks have to fit in one line !

    s = size(tracks(1).states);

    % print double sided Feynman diag.
    % from top to bottom = backwards in time
    for id = (s(1)):-1: 1
        for id2 = 1:length(tracks)
            
              interaction = sign_to_char(tracks(id2).interaction.direction(id));
        
              if (tracks(id2).interaction.column(id) == 1)
                  sl = interaction;  sr = ' ';
              else
                  sl = ' ';           sr = interaction;
              end
              
              % start quantum numbers from 0, although array index start
              % from 1
              fprintf('%s| %d %d |%s       ',sl, tracks(id2).states(id,:) -1, sr);

        end  

       fprintf('\n');
    end

    % print phase matching signs
     fprintf('   '); 
     for id2 = 1:length(tracks)
          
         for id = 2:s(1)
              fprintf('%s', sign_to_char(tracks(id2).phasematch(id)));
         end  
        fprintf('            ');   
     end
    fprintf('\n');
end


% Helper function -------

function s = sign_to_char(value)
% sign as char / string

    s = ' ';
    if (value > 0)   s = '+';  end
    if (value < 0)   s = '-';  end
end
    