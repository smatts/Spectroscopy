function display_tracks_by_phasematch(tracks, phasematch)
% display only those tracks that fullfuill the phasematching condition

    % select the tracks
    tr = tracks(select_by_phasematch(tracks, phasematch));

    % distribute tracks in lines of 4
    tracks_per_line = 4;
    n = fix(length(tr) ./ tracks_per_line);
    rem = length(tr) - n .* tracks_per_line;

    for id=0:(n-1)
        display_track(tr((id*tracks_per_line +1):(id*tracks_per_line + tracks_per_line)))
        fprintf('\n');
    end

    if (rem > 0)
        display_track(tr(end- rem +1:end))
    end
    
    fprintf('\n\n');
end