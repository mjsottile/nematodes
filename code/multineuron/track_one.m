%
% be sure to addpath(genpath(/where/nematodes/code/is))
% to pick up the right kmeans
%
% call with image sequence and then a 1x2 vector [col, row] of the
% point of interest to track
%
% mjsottile@gmail.com
%
function track_one(im, seed_coords)
    n = length(im);

    % manually find column and row in first frame of neuron to
    % track
    % track_point = [507, 123];
    track_point = seed_coords;    

    for i=1:n
        % convert to double
        img = double(im{i});
       
        
        % kmeans segmentation: split into 4 segments, 35 iterations.
        % 35 was selected by observing when the algorithm appeared to
        % converge for the sample data I was provided.
        [km,c] = kmeans(img,4,35);
        
        % pick the segment that holds the brightest pixels
        bright = km==4;
        
        % compute centroid and area for each independent component in
        % the segmented image
        s = regionprops(bright, 'centroid', 'area');
        
        centroids = cat(1,s.Centroid);
        areas = cat(1,s.Area);

        % filter out tiny segments that aren't neurons
        big_enough = find(areas > 10);
        centroids = centroids(big_enough,:);

        % compute distance of track point to all centroids to find
        % the one nearest to it.
        d_centroids = (centroids - repmat(track_point,size(centroids,1),1)).^2;
        d_centroids = sqrt(sum(d_centroids,2));

        %disp(d_centroids)
        [~,nearest_idx] = min(d_centroids);
        
        % update track_point with centroid that was closest
        track_point = centroids(nearest_idx,:);
        
        % DEBUG: plot
        cap = sprintf('Image %04d',i);
        
        imagesc(img);colormap(gray);
        hold on;
        plot(centroids(:,1),centroids(:,2),'y+');
        plot(track_point(:,1),track_point(:,2),'ro');
        hold off;
        title(cap);
        drawnow;
    end
