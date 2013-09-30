function [ratio,yfp,cfp,refthresh,centx,centy, nangle] = calcium_process(frames, varargin)
%
% calcium imaging processing code
%
% This routine implements the method described in the WormBook
% found here in section 3.7, protocol 4:
%
%   http://www.wormbook.org/chapters/www_imagingneurons/imagingneurons.html
%
% Extraction of the cells of interest and the non-fluorescent background
% is based on crude threshold segmentation.  Threshold-based object
% selection is used, with two thresholds (the first is required, the
% second is optional): one for the high end to define bright regions 
% representing cells of interest, and one for the low end to define dim 
% regions representing background.  
%
% We are assuming that the frame is split down the middle
% and rely on a helper function (splitter.m) to do the frame splitting
% and registration of the halves to correct for distortions that result
% from the imaging aparatus.  The distortions are assumed to be
% nothing more than a composition of translations and rotations, with
% length scales being unmodified.  A radius is provided for a circle
% representing the region of interest around the identified bright cells.
%
% Example usage:
%
%   [ratio, yfp, cfp] = calcium_process(im,500,'bthresh',15,'circle',25)
%
% where im is the cell array of frames, 500 is the threshold above which
% we treat pixels as being part of the cell.  bthresh indicates that
% background selection is enabled with a threshold of 15.  circle indicates
% that a circle-based ROI is to be used with a radius of 25.
%
% input:
%   frames  : sequence of frames to process as a cell array
%
%   varargin: optional value/param arguments.  Legal optional arguments
%             include:
%                'thresh', val  : threshold for high end, where pixels 
%                                 above the threshold are considered
%                                 part of fluorescing cells of interest.
%
%                'bthresh', val : enable background removal with the given
%                                 background threshold value, where
%                                 pixels with intensity < val will be
%                                 considered background.
%
%                'circle', val  : toggle a circular ROI centered on the
%                                 centroid of the detected cell with radius
%                                 specified by val. If no param is
%                                 specified, largest connected component is
%                                 used
%
% output:
%   sig     : the signal obtained by computing the ratio of the yellow
%             channel over the cyan channel, with corresponding backgrounds
%             subtracted.
%   yfp     : the signal from the yellow channel with background removed
%   cfp     : the signal from the cyan channel with background removed
%
% Matthew Sottile / November 2012
% mjsottile@gmail.com
%

    %
    % argument handling
    %
    p = inputParser;
    addRequired(p, 'frames');
    addParamValue(p, 'thresh', -1);
    addParamValue(p, 'bthresh', -1);
    addParamValue(p, 'circle', -1);
    addParamValue(p, 'rig', 1);
    parse(p, frames, varargin{:});
    
    % frames
    im = p.Results.frames;

    % flag to see if we use the circle ROI or the largest connected
    % component. If not specified in call to function, connected component
    % only is used. If 'circle' flag is used in call, then a circle of the
    % radius specified by the user and centered at the cenroid of the connected 
    % component will be used as the ROI. 
    if (p.Results.circle == -1)
        use_circle = 0;
    else
        use_circle = 1;
        radius = p.Results.circle;
    end
    
    % background subtraction. Set to 1 for background, 0 for not. If not
    % specified in function call, default is to subtract background. 
    if (p.Results.bthresh == or(-1,1))
        handle_background = 1; 
    else
        handle_background = 0;
        bthresh = p.Results.bthresh;
    end
    
%     %%% NOTE: EXPERIMENTAL CODE
%     % thresh estimate toggle
%     use_thresh_estimate = 1;
%     if (use_thresh_estimate == 1)
%         [te_bright, te_dim] = thresh_estimate(frames);
%         thresh = te_bright;
%     end
%     %%% END EXPERIMENTAL CODE
    
%% the rest of the code

    % signal to return
    ratio = zeros(1,length(im));

    % signal for each side individually
    yfp = zeros(1,length(im));
    cfp = zeros(1,length(im));
    
    % find a good reference frame
    % 0.25 is hardcoded constant - could be adjusted.  
    % based on guessing a threshold from a small set of evenly spaced
    % images from the sequence.  Hopefully the majority of these would
    % be considerd "good" frames.
    refthresh = 0;
    numtests = 150;
    testindices = (1:numtests)*floor(length(im)/numtests);
    disp('Sampling images for thresholding.');
    
    if (p.Results.thresh ==-1)
        for i=1:numtests
        [testt,~] = thresh_estimate(im{testindices(i)});
        refthresh = refthresh + testt;
        end   
    refthresh = refthresh / numtests;
    else refthresh=p.Results.thresh;
    end
    
    
    refframe = find_goodframe(im, refthresh, 0.25);
    
    % register to obtain transform.  discard registered frames since we will
    % re-register anyway later using the tform object.
    [~,~,tform] = splitter(double(im{refframe}));
    framesize = size(im{refframe});

    disp(refthresh) % should we move this to an output param?
    disp('Processing.');
    for i=1:length(im)
        % compute threshold for this frame
              
        [lhs,rhs] = splitter2(im{i},tform);
        
        [thresh, bthresh] = thresh_estimate(lhs);
        
        % maximum intensity
        maxval = max(lhs(:));

        % binary image with all pixels within thresh of the max intensity
        BWmax = lhs > thresh;
        disp(length(find(BWmax==1)))

        % compute connected components of thresholded regions
        CCmax = bwconncomp(BWmax);

        % compute centroid of each CC
        Smax = regionprops(CCmax,'Centroid');

        % if we had more than one CC, pick the biggest
        if (length(Smax) > 1)
            largest_max = 0;
            largest_max_size = 0;
            for j=1:length(CCmax.PixelIdxList)
                len = length(CCmax.PixelIdxList{j});
                if (len > largest_max_size)
                    largest_max_size = len;
                    largest_max = j;
                end
            end
        else
            % if we only had 1, then index 1 is the biggest
            largest_max = 1;
        end
        
        % define the mask used to capture the neuron of interest
        if (use_circle == 1)
            % make the circle mask.
            [size_rows,size_cols] = size(lhs);
            [nc,nr] = ndgrid(1:size_rows,1:size_cols);
            lhs_mask = sqrt( ...
                (nc-Smax(largest_max).Centroid(2)).^2 + ...
                (nr-Smax(largest_max).Centroid(1)).^2) ...
                < radius;
            rhs_mask = lhs_mask;
        else
            % try to make a mask based on the largest connected
            % component.
            BWmax = double(BWmax);
            BWmax(CCmax.PixelIdxList{largest_max}) = 2;
            lhs_mask = double(BWmax == 2);
            strel_size = 15;
            se = strel('ball',strel_size,strel_size);
            lhs_mask = imdilate(lhs_mask,se) > strel_size;
            rhs_mask = lhs_mask;
        end
        
        % find all pixels that are below bthresh to define the background
        if (handle_background == 1)
          minval = min(lhs(:));
          background_mask = lhs < bthresh;
          
          lhs_background = lhs .* double(background_mask);
          rhs_background = rhs .* double(background_mask);
        
          % compute lbkg and rbkg as average value of pixels determined
          % to be within what we call the background region
          lbkg = sum(lhs_background(:))/length(find(background_mask));
          rbkg = sum(rhs_background(:))/length(find(background_mask));
        end          
      
        % mask the halves
        lhs_masked = lhs.*double(lhs_mask);
        rhs_masked = rhs.*double(rhs_mask);

        % count masked pixels in each side
        lhs_nnz = length(find(lhs_masked(:) > 0));
        rhs_nnz = length(find(rhs_masked(:) > 0));

        % do R computation.  Note that we take the mean of only the masked 
        % pixels to get the average intensity.  Subtract from each the
        % average of the background pixels determined above if
        % handle_bakground is set.
        if (p.Results.rig==1) % if rig 1, lhs is yfp. if rig 2, lhs is cfp
            if (handle_background == 1)
                yfp(i) = (sum(lhs_masked(:))/lhs_nnz) - lbkg;
                cfp(i) = (sum(rhs_masked(:))/rhs_nnz) - rbkg;
            else
                yfp(i) = (sum(lhs_masked(:))/lhs_nnz);
                cfp(i) = (sum(rhs_masked(:))/rhs_nnz);
            end
        else
            if (handle_background == 1)
                yfp(i) = (sum(rhs_masked(:))/rhs_nnz) - rbkg;
                cfp(i) = (sum(lhs_masked(:))/lhs_nnz) - lbkg;
            else
                yfp(i) = (sum(rhs_masked(:))/rhs_nnz);
                cfp(i) = (sum(lhs_masked(:))/lhs_nnz);
            end
        end
        ratio(i) = yfp(i)/cfp(i);
            
        % correct for R_CFP.  Paper recommends value of 0.6 if it wasn't
        % measured directly.  NOTE: we are not doing the f_bkg correction,
        % which is necessary to remove signal due to background noise versus
        % actual flourescence from the cell.
        ratio(i) = ratio(i)-0.6;

        %return centroid data
        centx(i) = Smax(largest_max).Centroid(1);
        centy(i) = Smax(largest_max).Centroid(2);

        % plotting stuff

        % LHS plot
        title(sprintf('%d',i));
        subplot(2,2,1);
        imagesc(lhs_masked);
        colormap gray;
        hold on;
        plot(Smax(largest_max).Centroid(1),Smax(largest_max).Centroid(2),'r+');
        hold off;

        % RHS plot    
        subplot(2,2,2);
        imagesc(rhs_masked);
        colormap gray;
        hold on;
        plot(Smax(largest_max).Centroid(1),Smax(largest_max).Centroid(2),'r+');
        hold off;

        % signal so far
        subplot(2,2,3:4);
        plot(ratio);

        drawnow;
    end
    orientation = p.Results.rig;
    [nangle] = neuron_angle(framesize, centx, centy,orientation);
end
