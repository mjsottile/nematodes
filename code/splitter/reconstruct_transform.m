function [theta,offsetx,offsety] = reconstruct_transform(regframe)

  % step 1: find some edges.  ideally we will have two dominant edges from
  %         the registration process, and can find those via the Hough
  %         transform.
  
  % create binary image where edges of registered frame are apparent
  BW = edge(regframe);

  % step 2: perform the hough transform.
  [H,theta,rho] = hough(BW);
  
  % step 3: extract the peaks from the Hough transform matrix that
  %         correspond to the strongest lines that were found in the
  %         image.
  P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));

  % step 4: find the lines that were identified by the houghpeaks
  %         functions in a form that we can use to analyze their
  %         position and derive the transform from them that was
  %         used for the registration
  lines = houghlines(BW,theta,rho,P,'FillGap',50,'MinLength',60);

  % now, find the two lines we care about - the longest mostly vertical
  % line and the longest mostly horizontal line.  we are assuming that
  % registration wasn't too extreme and is ONLY due to the imaging
  % effects that are known to be caused in the nematodes case we are
  % working with.  this algorithm is NOT general purpose for determining
  % the parameters of the transform found by imregister.  beware of reuse.
  
  % sort the lines by length
  linelen = zeros(1,length(lines));
  for i=1:length(lines)
      linelen(i) = norm(lines(i).point1 - lines(i).point2);
  end
  
  [~,idx] = sort(linelen,'descend');
  
  % take the two longest lines
  line1 = lines(idx(1));
  line2 = lines(idx(2));
  
  % for the longest line, compute its difference in x and y
  xy = [line1.point1; line1.point2];
  xdiff = abs(xy(1,1) - xy(2,1));
  ydiff = abs(xy(1,2) - xy(2,2));
  
  % if it ran longer in the x dimension than y, assume it is the
  % horizontal line and the other one is vertical.  otherwise,
  % assume the opposite.
  if (xdiff > ydiff)
      horizline = line1;
      vertline = line2;
  else
      vertline = line1;
      horizline = line2;
  end
  
  % compute the same angle for both.  average them to get the
  % theta we will use to remove pixel-level effects that cause
  % the angles to be slightly different.
  theta1 = atan(abs(vertline.point1(1)-vertline.point2(1)) / ...
                abs(vertline.point1(2)-vertline.point2(2)));
  theta2 = atan(abs(horizline.point1(2)-horizline.point2(2)) / ...
                abs(horizline.point1(1)-horizline.point2(1)));
            
  theta = (theta1+theta2)/2; 

  x1 = vertline.point1(1);
  x2 = vertline.point2(1);
  y1 = vertline.point1(2);
  y2 = vertline.point2(2);
  
  % want to do rotation in clockwise fashion, so if the rotation
  % was counter clockwise, change sign of theta to put it into
  % clockwise orientation
  if (y1<y2 && x1<x2)
      theta = -1.0 * theta;
  end
  
  x3 = horizline.point1(1);
  x4 = horizline.point2(1);
  y3 = horizline.point1(2);
  y4 = horizline.point2(2);
  
  % find where the two lines meet.  this will be the lower right
  % corner of the registered image inside the frame.  from this
  % we can approximate the translation that was performed.
  % TODO: this is not likely correct since rotation was probably
  % performed by imregister relative to the center of the rectangle,
  % not the corner.  Hopefully the error is small and will put our
  % target off by only a pixel or so, versus the 10 or so that we 
  % get due to the optics.
  %
  % line intersection equation : 
  %     http://en.wikipedia.org/wiki/Line-line_intersection
  px = round(((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4)) / ...
        ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)));
    
  py = round(((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4)) / ...
        ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)));
  
  [height,width]=size(regframe);

  offsetx = px-width;
  offsety = py-height;
