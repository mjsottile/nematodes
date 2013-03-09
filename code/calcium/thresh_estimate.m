function [t,bt] = thresh_estimate(im)
    % estimate the thresholds to use for identifying the neuron of
    % interest as well as the background from the given image sequence.
    %
    % the criteria for selection are:
    %   - for each frame, use k-means clustering to partition the
    %     intensities into two groups (bright, dim)
    %   - for each frame, record the minimum intensity in the bright
    %     region and the cluster center for the dim region.
    %   - set the high-end threshold as the mean of the minimum
    %     bright region intensity minus one standard deviation of the
    %     minimum intensity over the sequence.
    %   - set the background threshold to 10% of the mean dim region
    %     cluster center.
    %
    numim = length(im);
    
    centers = zeros(1,numim);
    bright_min = zeros(1,numim);
    dim_ctr = zeros(1,numim);
    
    for i=1:numim
        [u,c] = kmeans(im{i},2,10);
        
        masked = (u==2).*im{i};
        bright_min(i) = min(find(masked(:)));
        
        dim_ctr(i) = c(1);

        if (mod(i,10)==0)
            disp(i)
        end
        
    end

    t = mean(bright_min) - std(bright_min);
    bt = mean(dim_ctr) * 0.1;
    
