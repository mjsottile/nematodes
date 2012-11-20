%
% code to split an image down the center and return the two halves
% with the RHS registered with the LHS
%
function [lhs, rhsRegistered] = splitter(im)
  [rows,cols] = size(im);

  lhs = im(:,1:(cols/2));
  rhs = im(:,((cols/2)+1):end);

%  subplot(1,2,1);
%  imshowpair(lhs,rhs,'Scaling','joint');
  
  [optimizer, metric] = imregconfig('monomodal');
  
  optimizer.MaximumIterations = 50;
  optimizer.MinimumStepLength = 5e-3;
  
  rhsRegistered = imregister(rhs, lhs, 'rigid', optimizer, metric);
  
%  subplot(1,2,2);
%  imshowpair(lhs,rhsRegistered,'Scaling','joint');
%  drawnow;