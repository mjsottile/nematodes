%Calcium Imaging Analysis gui. The images are loaded from the chosen directory. The
%first image is diplayed and the progam waits for the user to choose the
%neurons of interest. The user clicks in the center of the neurons s/he
%would like to analyze and then presses return. 
%
%AUTHOR: Kat McCormick
%        kat.mccormick@gmail.com


%% Load the directory
%ims = ui_load_directory;


%% Get user input on the neuron locations
figure; 
imagesc(ims{1}(:,1:320));
disp('Choose neurons to analyze. Hit return when finished with selections.')
[xs,ys] = ginput;
seed_coords = [xs ys];
coord_nums = size(seed_coords)
track_points = ones(coord_nums(1), coord_nums(2), length(ims));
track_points(:,:,1) = seed_coords;


%% Find the neurons in the images. 
num_ims = length(ims);
for i=1:num_ims
        % convert to double
        img = double(ims{i});
       
        
        % kmeans segmentation: split into 4 segments, 35 iterations.
        % 35 was selected by observing when the algorithm appeared to
        % converge for the sample data I was provided.
        [km,c] = kmeans(img,5,35);
        
        % pick the segment that holds the brightest pixels
        bright = km>=4;
        
        % compute centroid and area for each independent component in
        % the segmented image
        s = regionprops(bright, 'centroid', 'area');
        
        centroids = cat(1,s.Centroid);
        areas = cat(1,s.Area);

        % filter out tiny segments that aren't neurons
        centroids = centroids(areas>10,:);

        % compute distance of all track points to all centroids to find
        % the one nearest to each track point.
        dist_matrix = pdist2(track_points(:,:,i), centroids);

        % find the row value of the minimum distance
        [~,nearest_idx] = min(dist_matrix');
        
        % update track_point with centroid that was closest
        track_points(:,:,i+1) = centroids(nearest_idx,:);
        
        % DEBUG: plot
        cap = sprintf('Image %04d',i);
        
        imagesc(img);colormap(gray);
        hold on;
        plot(centroids(:,1),centroids(:,2),'y+');
        plot(track_points(:,1,i),track_points(:,2,i),'ro');
        hold off;
        title(cap);
        drawnow;
    end
    

