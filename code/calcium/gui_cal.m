%Calcium Imaging Analysis gui. The images are loaded from the chosen directory. The
%first image is diplayed and the progam waits for the user to choose the
%neurons of interest. The user clicks in the center of the neurons s/he
%would like to analyze and then presses return. 
%
%AUTHOR: Kat McCormick
%        kat.mccormick@gmail.com


%% Load the directory
ims = load_directory;

%% Get the transformation from left hand side (LHS) to right hand side (RHS)
refthresh = 0;
    numtests = 150;
    testindices = (1:numtests)*floor(length(ims)/numtests);
    disp('Sampling images for thresholding.');
    
 
    for i=1:numtests
        [testt,~] = thresh_estimate(ims{testindices(i)});
        refthresh = refthresh + testt;
    end
    refthresh = refthresh / numtests;
   
    
    
    refframe = find_goodframe(ims, refthresh, 0.25);
    
    % register to obtain transform.  discard registered frames since we will
    % re-register anyway later using the tform object.
    [~,~,tform] = splitter(double(ims{refframe}));
    framesize = size(ims{refframe});


%% Get user input on the neuron locations in the first frame
figure; 
[lhs,rhs] = splitter2(ims{1},tform);
imagesc(lhs);
disp('Choose neurons to analyze. Hit return when finished with selections.')
[xs,ys] = ginput;
seed_coords = [xs ys];
coord_nums = size(seed_coords)
track_points = ones(coord_nums(1), coord_nums(2), length(ims));
track_points(:,:,1) = seed_coords;


%% Find the neurons in the images. 
num_ims = length(ims);
for i=1:num_ims
        [lhs,rhs] = splitter2(ims{i},tform);     
        
        % kmeans segmentation: split into 3 segments, 10 iterations.
        % This is not convergent, but does put neurons in top segment on
        % the test data set
        [km,c] = kmeans(lhs,4,10);
        
        % pick the segment that holds the brightest pixels
        bright = km>=3;
        
        % compute centroid and area for each independent component in
        % the segmented image
        reg_lhs = regionprops(bright,lhs, 'centroid', 'area','meanintensity', 'pixelIdxList');
        reg_rhs = regionprops(bright,rhs, 'meanintensity', 'pixelIdxList');

  
        centroids = cat(1,reg_lhs.Centroid);
        areas = cat(1,reg_lhs.Area);
        lhs_intensities = cat(1, reg_lhs.MeanIntensity);
        rhs_intensities = cat(1, reg_rhs.MeanIntensity);

        % filter out tiny segments that aren't neurons
        centroids = centroids(areas>10,:);
        lhs_intensities = lhs_intensities(areas>10);
        rhs_intensities = rhs_intensities(areas>10);

        % compute distance of all track points to all centroids to find
        % the one nearest to each track point.
        dist_matrix = pdist2(track_points(:,:,i), centroids);

        % find the row value of the minimum distance
        [~,nearest_idx] = min(dist_matrix');
        
        % update track_point with centroid that was closest
        track_points(:,:,i+1) = centroids(nearest_idx,:);
        lhs_int(i,:) = lhs_intensities(nearest_idx,:);
        rhs_int(i,:) = rhs_intensities(nearest_idx,:);
        
        %get ratios for the spots. This setup is for rig 2. If rig 1, then
        % use lhs/rhs. 
        ratio(i,:) = rhs_int(i,:)./lhs_int(i,:);
        
        % DEBUG: plot
        cap = sprintf('Image %04d',i);
        subplot(2,1,1);
        imagesc(lhs);colormap(gray);
        hold on;
        plot(centroids(:,1),centroids(:,2),'y+');
        plot(track_points(:,1,i+1),track_points(:,2,i+1),'ro');
        hold off;
        title(cap);
        drawnow;
        subplot(2,1,2);
        plot(i, ratio(i,:))
        hold on;
    end
    

