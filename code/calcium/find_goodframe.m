function refframe = find_goodframe(ims, thresh, frac)
    
    numim = length(ims);
    maxval = zeros(numim,1);
    nummasked = zeros(numim,1);
    
    % for each frame, compute the maximum pixel value and the
    % number of pixels that are above the provided threshold
    for i=1:numim
        maxval(i) = max(ims{i}(:));
        masked = ims{i}(:) > thresh;
        nummasked(i) = length(find(masked(:)));
    end
    
    % sort the maxval array to get an ordering of images by peak
    % brightness
    [~,idx] = sort(maxval);
    
    % take the top frac% best frames to select from
    % assume frac is between 0 and 1.  A value like 0.25 is
    % ideal (0.25 would capture the top 25%).
    cutoff = numim - round(numim * frac);
    
    % take the indices of the brightest frames
    bright_indices = idx(cutoff:end);
    
    % sort the frames indexed by bright_indices
    [~,idx2] = sort(nummasked(bright_indices));
    
    % best one is at the end
    best_idx2 = idx2(end);
    
    % since idx2 was computes from a SHORTER array than idx given that
    % we sorted only those frames from the bright index set, we need
    % to correct the best_idx2 value to map to the original index set.
    best = best_idx2 + cutoff;
    refframe = idx(best);
    
%   % broken code below
%
%     % Get the indices of the top frac of frames with the most pixels above
%     % threshold. 
%     [pixsort, idx2] = sort(nummasked);
%     topidx = idx2(cutoff:end);
%     refframe = idx2(end);
%     
%     % find the biggest connected component in the frames we consider
%     % to be in the best frac% of the frames.
%     %I can't figure out what this was supposed to be doing, but I
%     %think it was doing it wrong. KM
%     
% %     [~, idx2] = sort(nummasked((numim-cutoff):end));
% %     
% %     best = idx2(end)+(numim-cutoff-1);
% %     refframe = idx(best);
%     