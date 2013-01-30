%
% code to split an image down the center and return the two halves
% with the RHS registered with the LHS
%
function [lhs, rhsRegistered] = splitter2(im, tform)
  [rows,cols] = size(im);

  lhs = im(:,1:(cols/2));
  rhs = im(:,((cols/2)+1):end);

  rhsRegistered = imtransform(rhs,tform,...
      'XData',[1 (cols/2)],'YData',[1 rows]);
