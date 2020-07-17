function ok = select_by_phasematch(tracks, phasematch)
% select those tracks that fullill phasematching condition
% phasematch: array of (+- 1)

    n = length(phasematch);
    
    for id=1:length(tracks)
         ok(id) = all(tracks(id).phasematch(2:2+n-1) == phasematch);
    end


end