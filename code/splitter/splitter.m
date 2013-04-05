%
% code to split an image down the center and return the two halves
% with the RHS registered with the LHS
%
function [lhs, rhsRegistered, tform] = splitter(im)
  [~,cols] = size(im);

  lhs = im(:,1:(cols/2));
  rhs = im(:,((cols/2)+1):end);

%  subplot(1,2,1);
%  imshowpair(lhs,rhs,'Scaling','joint');
  
  [optimizer, metric] = imregconfig('monomodal');
  
  optimizer.MaximumIterations = 50;
  optimizer.MinimumStepLength = 5e-3;
  
  [tform_string,rhsRegistered] = evalc('imregister(rhs, lhs, ''rigid'', optimizer, metric, ''DisplayOptimization'', true)');
  
  % this is a horrendous hack!!!!
  strtok(tform_string,'imtransform.');
  [~,ss] = regexp(tform_string,'imtransform.','match','split');
  eval(strcat([ss{2},'; tform = T']));
  
  % find region that is zero due to registration fill in
  rhsMask = 1-double(rhsRegistered == 0);

  % zero out pixels in the registration region on the LHS as well
  lhs = lhs .* rhsMask;
  
%  subplot(1,2,2);
%  imshowpair(lhs,rhsRegistered,'Scaling','joint');
%  drawnow;
