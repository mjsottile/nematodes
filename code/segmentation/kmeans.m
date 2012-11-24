function [u,c] = kmeans(f,k,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [u,c] = kmeans(f,k,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function carries out k-means clustering procedure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   f = Given image to be segmented.
%   k = Number of clusters (regions) desired.
%   N = Number of iterations to take.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%   u = A function that takes k values, indicating the
%       different regions found.
%   c = A vector containing the cluster values (average
%       grayscale intensities) of the regions found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_default=2;
N_default=10;
if nargin<3;N=N_default;end
if nargin<2;k=k_default;end
if isempty(N);N=N_default;end
if isempty(k);k=k_default;end

% Initial choice of vector c: Uniform distribution of values.
fmax = max(max(f));
fmin = min(min(f));
c = 1:1:k;
c = fmin + (fmax-fmin).*c/k;

onez = ones(size(f));

% Initial allocation of pixels:
u = onez;
m = abs(f-c(1));
for i=2:k
   mnew = abs(f-c(i));
   u = u + ( (m-mnew)>0 );
   m = mnew;
end

% n = size(f,1)*size(f,2);
% F = sum(sum(f));

for t=1:N   % Main iteration starts here.

   % Calculate the new averages:
   for i=1:k
      ind = (u==i);
      U = sum(sum(ind));
      FU = sum(sum(ind.*f));
      c(i) = FU/(U+1e-10);
   end

	c = sort(c);

	% Relocate the pixels:
	u = onez;
	for i=2:k
		d = f-(c(i-1)+c(i))/2;
		d = d>0;
   	u = u + d;
	end

end         % Main iteration ends here.

