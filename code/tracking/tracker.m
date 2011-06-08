function [xs,ys,zs,adj] = tracker(ims, thresh, areaCutoff)

  close all;

  % how many images
  numim = length(ims);
  
  % split images in half
  for i=1:numim
      imLeft{i} = ims{i}(:,1:168);
      imRight{i} = ims{i}(:,169:end);
  end
  
  % binarize images
  for i=1:numim
      % binarize right side with threshold
      binIms{i} = imRight{i} > thresh;
      
      % plot
      subplot(1,2,1);imagesc(binIms{i});
            
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

      drawnow;
  end
  close all;
  
  % triangulation stuff
  
  multiplier = 5;
  xs = cents{1}(1,:);
  ys = cents{1}(2,:);
  zs = ones(1,length(cents{1}(1,:)));
  
  for i=2:numim
    if (~isempty(cents{i}))
        xs = [xs cents{i}(1,:)];
        ys = [ys cents{i}(2,:)];
        zs = [zs (ones(1,length(cents{i}(1,:))).*i.*multiplier)];
    end
  end  
  
%  close all;
%  plot3(xs,ys,zs,'.');
%  drawnow;
  
  tri = delaunay(xs,ys,zs);
  
  % tets to lines
  tri = [tri(:,1) tri(:,2); ...
         tri(:,2) tri(:,3); ...
         tri(:,3) tri(:,4); ...
         tri(:,4) tri(:,1)];
       
  lines = [0 0];

  [nsegs,~] = size(tri);
  
  adj = sparse(numim,numim);
  
  h = waitbar(0);
  for i=1:nsegs
    na=tri(i,1);
    xa=xs(tri(i,1));
    ya=ys(tri(i,1));
    za=zs(tri(i,1));
    nb=tri(i,2);
    xb=xs(tri(i,2));
    yb=ys(tri(i,2));
    zb=zs(tri(i,2));    
    
    if (za ~= zb)
        if (za < zb)
            adj(na,nb) = sqrt((xa-xb)^2 + (ya-yb)^2 + (za-zb)^2);
        else
            adj(nb,na) = sqrt((xa-xb)^2 + (ya-yb)^2 + (za-zb)^2);
        end
    end
    
    waitbar(i/nsegs);
  end
  close(h);
