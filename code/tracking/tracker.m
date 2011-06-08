function [xs,ys,zs,adj] = tracker(ims, thresh, areaCutoff)

  close all;

  % how many images
  numim = length(ims);
  
  % split images in half
  [~,imwidth] = size(ims{1});
  for i=1:numim
    imLeft{i} = ims{i}(:,1:(imwidth/2));
    imRight{i} = ims{i}(:,((imwidth/2)+1):end);
  end
  
  % binarize images
  for i=1:numim
      % binarize right side with threshold
      binIms{i} = imLeft{i} > thresh;
      
      % plot
      subplot(1,2,1);imagesc(binIms{i});
      se = strel('disk',3);        
      binIms{i} = imerode(binIms{i},se);
            
      % label connected components, compute their area
      lbl = bwlabel(binIms{i});
      props = regionprops(lbl, 'Area');
      areas = [props.Area];
      
      % filter connected components by the area cutoff
      areaThresh = areas > areaCutoff;
      tmp = zeros(size(binIms{i}));
      for j=1:max(lbl(:))
          if (areaThresh(j) == 1)
              tmp = tmp + double(lbl==j);
          end
      end
      
      % store the new image that we computed with filtered components
      binIms{i} = tmp;
      subplot(1,2,2);
      imagesc(binIms{i});

      % find the centroid of the preserved components
      props = regionprops(bwlabel(binIms{i}), 'Centroid');
      centroids = [props.Centroid];

      % reshape the centroids to be an array where the first row is the
      % x-coordinates, second row is the y-coordinates, columns are the
      % individual centroids
      numcents = length(centroids) / 2;
      centroids = reshape(centroids,[2 numcents]);
      hold on;
      plot(centroids(1,:),centroids(2,:),'y.');
      hold off;
  
      % store centroid
      cents{i} = centroids;
      
      disp(100*(i/numim));

      drawnow;
  end
  close all;
  
  % triangulation stuff
  
  multiplier = 0.5;

  xs = cents{1}(1,:);
  ys = cents{1}(2,:);
  zs = ones(1,length(cents{1}(1,:))).*multiplier;
  
  for i=2:numim
    if (~isempty(cents{i}))
        xs = [xs cents{i}(1,:)];
        ys = [ys cents{i}(2,:)];
        zs = [zs (ones(1,length(cents{i}(1,:))).*i.*multiplier)];
    end
  end  
  
  tri = delaunay(xs,ys,zs);
  
  % tets to lines
  tri = [tri(:,1) tri(:,2); ...
         tri(:,2) tri(:,3); ...
         tri(:,3) tri(:,4); ...
         tri(:,4) tri(:,1)];
  
  % want to make sure lines are sorted such that the first point is always
  % of a lower Z value than the second point

  p1z = zs(tri(:,1));
  p2z = zs(tri(:,2));
  
  mask = p1z < p2z;
  notmask = not(mask);

  tri = [tri(mask,:);
         fliplr(tri(notmask,:))];

  % now, sort the lines in ascending order of z coordinate of first point
  tzs = zs(tri(:,1));
  [~,indices] = sort(tzs);
  tri = tri(indices,:);

  % ok, now we're set up with an array of lines that is tuned for quickly
  % creating the sparse adjacency matrix.
  lines = [0 0];

  [nsegs,~] = size(tri);
  
  adj = sparse(1,1);
  
  na = tri(:,1);
  nb = tri(:,2);
  xa = xs(na);
  ya = ys(na);
  za = zs(na);
  xb = xs(nb);
  yb = ys(nb);
  zb = zs(nb);

  dists = sqrt((xa-xb).^2 + (ya-yb).^2 + (za-zb).^2);
  for i=1:length(na)
      adj(na(i),nb(i))=dists(i);
  end

