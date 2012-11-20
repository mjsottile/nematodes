function sig = calcium_process(thresh, radius, im)
    % signal to return
    sig = zeros(1,length(im));
    
    for i=1:length(im)
        [lhs,rhs] = splitter(double(im{i}));
        
        % maximum intensity
        maxval = max(lhs(:));

        % binary image with all pixels within thresh of the max intensity
        BWmax = lhs>(maxval-thresh);

        % compute connected components of thresholded regions
        CCmax = bwconncomp(BWmax);

        % compute centroid of each CC
        Smax = regionprops(CCmax,'Centroid');

        % if we had more than one CC, pick the biggest
        % MAX
        if (length(Smax) > 1)
            largest_max = 0;
            largest_max_size = 0;
            for j=1:length(CCmax.PixelIdxList)
                len = length(CCmax.PixelIdxList(j));
                if (len > largest_max_size)
                    largest_max_size = len;
                    largest_max = j;
                end
            end
        else
            % if we only had 1, then index 1 is the biggest
            largest_max = 1;
        end

        % make the circle mask.  
        [sr,sc] = size(lhs);
        [nc,nr] = ndgrid(1:sr,1:sc);
        circlemask_max = sqrt((nc-Smax(largest_max).Centroid(2)).^2 + ...
                          (nr-Smax(largest_max).Centroid(1)).^2) < radius;

        % mask the halves
        lhs_masked = lhs.*double(circlemask_max);
        rhs_masked = rhs.*double(circlemask_max);

        % do R computation.  Note that we take the mean of only the masked 
        % pixels to get the average intensity.
        sig(i) = ...
            (sum(lhs_masked(:))/length(find(circlemask_max))) / ...
            (sum(rhs_masked(:))/length(find(circlemask_max)));

        % correct for R_CFP.  Paper recommends value of 0.6 if it wasn't
        % measured directly.  NOTE: we are not doing the f_bkg correction,
        % which is necessary to remove signal due to background noise versus
        % actual flourescence from the cell.
        sig(i) = sig(i)-0.6;

        % plotting stuff

        % LHS plot
        subplot(1,3,1);
        imagesc(lhs_masked);
        colormap gray;
        hold on;
        plot(Smax(largest_max).Centroid(1),Smax(largest_max).Centroid(2),'r+');
        hold off;

        % RHS plot    
        subplot(1,3,2);
        imagesc(rhs_masked);
        colormap gray;
        hold on;
        plot(Smax(largest_max).Centroid(1),Smax(largest_max).Centroid(2),'r+');
        hold off;

        % signal so far
        subplot(1,3,3);
        plot(sig);

        drawnow;
    end