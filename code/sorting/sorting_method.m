%
%
% prerequisite to running this: the images have been loaded.  For example:
%
%   d = dir('*.jpg');
%   for i=1:length(d)
%     im{i} = imread(d(i).name);
%   end
%
% once that is done, you have two variables of interest: d, the list of
% images in the directory, and im, the cell array of each image loaded
% off of disk.
%
%



% set up the bounding box for the sorting area.  this was manually
% extracted from the schematic diagram of the sorting device by loading
% it, displaying it with imshow, and using the coordinates tool to pick
% out a rectangle.
sorting_region_rows = [525, 616];
sorting_region_cols = [547, 1166];

% compute the width of the region
region_width = sorting_region_cols(2)-sorting_region_cols(1);

% what do we do if we get tons of regions?  this happened with one
% image where it looked like it was all noise with no worms visible.
too_many_region_limit = 10;

% spin through the directory contents
for i=1:length(d)
    % clip out region of interest
    q{i} = im{i}(sorting_region_rows(1):sorting_region_rows(2), ...
                 sorting_region_cols(1):sorting_region_cols(2), :);

    % just the green channel
    gchan{i} = q{i}(:,:,2);
    
    % compute a threshold for the channel and extract the corresponding
    % binary image
    level = graythresh(gchan{i});
    BW{i} = im2bw(gchan{i},level);
    
    % dilate the bw image to connect very close blobs together
    BW{i} = imdilate(BW{i},strel('disk',3,0));
    
    % compute areas and centroids
    rp = regionprops(BW{i}, 'Area', 'Centroid');
    
    % extract centroids
    n = length(rp);
    if (n > too_many_region_limit)
        disp('Probably garbage.');
        continue; % skip rest of loop body, go to next iteration
    end
    
    % pull the centroids out into a 2xN matrix, where we have N
    % centroids
    centroids = reshape([rp.Centroid], [2,n]);

    % compute the distance between all of the centroids
    xs = repmat(centroids(1,:),[n,1]);
    ys = repmat(centroids(2,:),[n,1]);
    dists = sqrt((xs-xs').^2 + (ys-ys').^2);
    
    % find the biggest distance
    biggest_dist = max(dists(:));
    
    % basic heuristic: males have at least one pair of regions that
    % are separated by more than 1/3 of the width of the region
    % where worms are sorted
    if (biggest_dist > region_width * 0.3)
        disp(sprintf('Male : %s',d(i).name));
    else
        disp(sprintf('Female : %s',d(i).name));
    end
        
%    subplot(3,3,i);
%    imagesc(BW{i});
%    title(sprintf('%f :: %s',biggest_dist,d(i).name));
end
