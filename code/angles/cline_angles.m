%
% head angles
%
% INPUT:
%
%   clines : a cell array for N frames, in which each cell represents the
%            centerline of M(i) points for frame i, stored as an Mx2
%            matrix
%   npts   : number of segments to split the center line into.  This
%            is usually small, like 7 or so.  That will yield
%            (npts-2) points connecting segments where we can compute an
%            angle.
%
% OUTPUT:
%
%   h : an N by (npts-2) matrix in which entry (i,j) represents the angle
%       at the segment joint between segment j and j+1 for frame number i
%
% AUTHOR: matt sottile
%         mjsottile@me.com // msottile@uoregon.edu
%
function h = hangles(clines,npts)

%
% compute angles in triangle formed by three points assuming the
% following point and angle naming scheme.  For reference, see:
%
% http://en.wikipedia.org/wiki/Triangle#Sine.2C_cosine_and_tangent
%
%            C
%
%            ^
%           / \
%          /   \
%         /gamma\
%      b /       \  a
%       /         \
%      /           \
%     /             \
%    / alpha    beta \
% A +-----------------+ B
%             c
%
% Matthew Sottile // mjsottile@me.com
%

n = length(clines);

h = zeros(n,npts-2);

for i=1:n
    l = length(clines{i});
    
    spacing = floor(l / npts);
    
    if (spacing < 1)
        continue;
    end
    
    for j = 1:(npts-2)
        % point A
        pA = clines{i}(1+((j-1)*spacing),:);
        
        % point B
        pB = clines{i}(1+(j*spacing),:);
        
        % point C
        pC = clines{i}(1+((j+1)*spacing),:);
        
        % PLOTTING CODE FOR DEBUGGING
%         plot([pA(1),pB(1),pC(1),pA(1)],[pA(2),pB(2),pC(2),pA(2)],'r');
%         hold on;
%         plot(pA(1),pA(2),'b+');
%         plot(pC(1),pC(2),'bo');
%         hold off;
%         drawnow;
        
        
        % side a
        a = sqrt(sum((pC-pB).^2));
        
        % side b
        b = sqrt(sum((pA-pB).^2));
        
        % side c
        c = sqrt(sum((pA-pC).^2));
        
        alpha = acos((b^2 + c^2 - a^2)/(2*b*c));
        beta = acos((a^2 + c^2 - b^2)/(2*a*c));
        gamma = acos((a^2 + b^2 - c^2)/(2*a*b));
        
        side = sign(det([(pC-pA)' (pB-pA)']));
        
        h(i,j) = side*(180-(180*(gamma/pi)));  
    end
        disp(i);
end
