function refframe = find_goodframe(ims, thresh, frac)
    
    numim = length(ims);
    
    % for each frame, compute the maximum pixel value and the
    % number of pixels in the connected component that is 
    % within the threshold of the max.
    for i=1:numim
        maxval(i) = max(ims{i}(:));
        
%%% NOTE: CHANGED TO EXPERIMENT WITH ABSOLUTE THRESHOLDS
        masked = ims{i}(:) > thresh;
%%%

        nummasked(i) = length(find(masked(:)));
    end
    
    % sort the maxval array
    [~,idx] = sort(maxval);
    
    % take the top frac% best frames to select from
    % assume frac is between 0 and 1.  A value like 0.25 is
    % ideal.
    
    cutoff = round(numim * frac);
    
    % find the biggest connected component in the frames we consider
    % to be in the best frac% of the frames.
    [~, idx2] = sort(nummasked((numim-cutoff):end));
    
    best = idx2(end)+(numim-cutoff-1);
    refframe = idx(best);
    