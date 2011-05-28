%
% code to perform a breadth first traversal from a computed
% starting position, and then select the longest path from that point to a
% tip of the skeleton.
%
% INPUTS:
%
%   I   : a binary image representing the skeletonized input image.
%
% OUTPUTS:
%
%   rmax      : row of the head
%   cmax      : column of the head
%   G         : the image of the centerline as pixels turned on
%   cline_seq : Mx2 matrix of rows and colums for all M points along
%               the center line
%
% AUTHOR : matt sottile
%          mjsottile@me.com // msottile@uoregon.edu
%
function [rmax,cmax,G,cline_seq] = dijk(I)

[rs,cs,~] = find(I);

% clip to remove excess background: force extra pixel boundary.
% assumes min/max row/col are not right at border...
Q = I(max([1 (min(rs)-1)]):max(rs)+1, min(cs)-1:max(cs)+1);

roffset = max([1 (min(rs)-1)]);
coffset = min(cs)-1;

% lazy.  get rows/cols again, but in clipped image
[rs,cs,~] = find(Q);

% sort
[rs,idx] = sort(rs);

%
% take some pixel in the lowest row, and the corresponding column
% this assumes that the worm is oriented verically, hanging DOWN
% from the clamp.  this is incorrect for other orientations.
%
r = rs(1);
c = cs(idx(1));

X = zeros(size(Q));
X(r,c) = 1;

% magic loop that performs the breadth first traversal in
% a vectorized form for high performance in matlab.
out = zeros(size(Q));
iter = 1;
while (isempty(find(X(:))) == false)
    out = out + X.*iter;
    iter = iter + 1;
    X = double(conv2(X,ones(3,3),'same') > 0).*(Q-double(out > 0));
end

% G now contains the skeleton, but not as binary but as a set of
% integers representing the distance of each pixel from the starting
% position (distance along skeleton, not euclidean or manhattan dist)
G = out;

% find all of the pixels that have a distance greater than 0
[rs,cs,v] = find(G);

% sort them
[~,idx] = sort(v);

% the head is the point furthest away on the graph.  add the row
% and column offsets computed above to compensate for clipping out
% the big sea of zeros around the skeleton.
rmax = rs(idx(end)) + roffset;
cmax = cs(idx(end)) + coffset;

% now backtrack along skeleton to point we came from finding shortest path
keepgoing = 1;

% un-adjust the head location back to being in the clipped region
rc = rmax-roffset;
cc = cmax-coffset;

% get the distance of the head
val = G(rc, cc);

% figure out the extent of the clipped region so we can make sure we
% don't accidentally try to walk out of it below.
[gr,gc] = size(G);

% mark value
vmark = val;

% first row of center line sequence is the head.
cline_seq = [rc cc];

% go until we decide to stop by some condition in the loop
while (keepgoing == 1)
    % create a 9 point stencil looking around the current row and col
    rsten = [rc-1 rc-1 rc-1 rc rc rc rc+1 rc+1 rc+1];
    csten = [cc-1 cc cc+1 cc-1 cc cc+1 cc-1 cc cc+1];
    
    % create a length 9 vector, one for each point in the stencil
    mask = zeros(1,9);
    for i=1:9
        % if the current stencil point is inside the image then...
        if (rsten(i) > 0 && rsten(i) < gr && ...
            csten(i) > 0 && csten(i) < gc)

            % set the mask value to a non-zero value if it is less than
            % the current distance value and greater than 0.
            mask(i) = G(rsten(i),csten(i)) * ...
                      double(G(rsten(i),csten(i)) < val) * ...
                      double(G(rsten(i),csten(i)) > 0);
        else
            % otherwise, mark this stencil point 0
            mask(i) = 0;
        end
    end
    
    % sort mask values in descending order
    [vals,idx] = sort(mask,'descend');
    
    % take the first one
    val = vals(1);
    
    % if it is zero, then we found nothing, so we're all done.
    if (val == 0)
        keepgoing = 0;
    else
        % otherwise, mark the location in G with the value to use for
        % marking
        G(rc,cc) = vmark;
        
        % now remember the new row and column from the stencil that
        % corresponds to the value we selected
        rc = rsten(idx(1));
        cc = csten(idx(1));
        
        % append a row with the next point we look at along the centerline
        cline_seq = [cline_seq; rc cc];
    end
end

% create a new matrix of the centerline, this time only including the
% points along the centerline versus the whole skeleton
newI = zeros(size(I))==1;
newI(roffset:(roffset+gr-1),coffset:(coffset+gc-1)) = (G==vmark);

% rename newI to G, which we return
G = newI;

