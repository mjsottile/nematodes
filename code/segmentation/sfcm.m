function [Qout, U, v] = sfcm(Qin, c, m, eps, p, q, width, sigma)
%
% spatial fuzzy c-means.
%
% This routine implements the spatial fuzzy c-means algorithm described
% by Chuang, et. al..  It requires six parameters, three of which are
% for traditional c-means, with the other three relevant for the spatial
% extension.  note that if q is zero, the algorithm defaults to traditional
% FCM.  Interesting choices of parameters are (m,p,q,width)=(2,1,1,5)
% or (2,0,2,5) or (2,1,0,5).  Please refer to either reference below, or
% just email Matt to figure it out.
%
% input:
%   Qin   : image data
%   c     : desired cluster count (DEFAULT=2)
%   m     : weighting parameter (DEFAULT=2)
%   eps   : epsilon for stopping (DEFAULT=0.1)
%   p     : pointwise membership importance parameter (DEFAULT=1)
%   q     : spatial membership importance parameter (DEFAULT=1)
%   width : width of the spatial weighting template.  must be odd!
%           (DEFAULT=5) 
%   sigma : variance^0.5 parameter for the 2d gaussian stencil.  a value 
%           of zero or lower will result in a flat stencil of 
%           ones(width,width) being used. (DEFAULT=0.5)
%
% output:
%   Qout : output clustering (crisp)
%   U    : membership matrix
%   v    : cluster centers
%
% references: Ross, Fuzzy logic + engineering applications, 2nd ed.
%
%             Chuang, K-S, et. al. 'Fuzzy c-means clustering with spatial
%              information for image segmentation.'  Computerized
%              Medical Imaging and Graphics, Vol. 30, pp 9-15, 2006.
%
% Matthew Sottile 
% matt@lanl.gov (5-6057) / April 2006
%

%%% defaults %%%
c_default = 2;
m_default = 2;
eps_default = 0.1;
p_default = 1;
q_default = 1;
width_default = 5;
sigma_default = 0.5;

if (nargin<8); sigma=sigma_default; end
if (nargin<7); width=width_default; end
if (nargin<6); q=q_default; end
if (nargin<5); p=p_default; end
if (nargin<4); eps=eps_default; end
if (nargin<3); m=m_default; end
if (nargin<2); c=c_default; end

if (isempty(sigma)); sigma=sigma_default; end
if (isempty(width)); width=width_default; end
if (isempty(q)); q=q_default; end
if (isempty(p)); p=p_default; end
if (isempty(eps)); eps=eps_default; end
if (isempty(m)); m=m_default; end
if (isempty(c)); c=c_default; end

%% done with defaults

if (width < 0)
    error('Width must be non-negative.');
end

if (width == (2*floor(width/2)))
    error('Width must be odd!');
end

%%
%% generate the kernel used by the convolution operator later for
%% computing the h_ij spatial function.  the Chuang paper employed
%% the ones(width,width) stencil, Matt prefers the 2D Gaussian.
%%
convborder = floor(width/2);
if (sigma <= 0)
    convkernel = ones(width,width);
else
    [x,y]=ndgrid(1:width,1:width);
    x=(x./width)-0.5;
    y=(y./width)-0.5;
    convkernel = (1./2.*pi.*sigma.*sigma).*exp(-((x).^2 + (y).^2)./ ...
                                               (2.*sigma.*sigma));
end

%%
%% step 0 : flatten image into a vector
%%
Q = Qin(:)';

%%
%% step 1 : initialize partition matrix.  k-means with a very low
%%          iteration count should get a good estimate.
%%
U = zeros(c,length(Q));
initial = kmeans(Q,c,5);

%% NOTE: some people suggest random assignment as a possible
%%       method.  this is a bad idea, as it is entirely possible that
%%       some combination of values will provide very bad segmentation
%%       and this will occur infrequently, and be a surprise.  to
%%       avoid this, avoid random initialization
%initial = ceil(rand(1,length(Q)).*c);

for i=1:length(Q)
    U(initial(i),i) = 1;
end

%
% two hacks to make sure we get through at least one iteration.
% first, we offset the centers to make sure the first distance
% computation is above eps, and we set dist to initially be big so
% the loop enters.
%
for i=1:c
    urow = U(i,:);
    v(i) = sum((urow.^m).*Q)/sum(urow.^m);
end
v = v+2;
dist = eps+1;

while (dist > eps)

    vold = v;
    %%
    %% step 2 : adjust cluster centers
    %%
    for i=1:c
        urow = U(i,:);
        v(i) = sum((urow.^m).*Q)/sum(urow.^m);
    end

    %%
    %% step 3 : update partition matrix
    %%
    for i=1:c
        d = abs(Q-v(i));
        accum = zeros(size(d));
        for j=1:c
            dj = abs(Q-v(j));
            quot = d./dj;
            quot = quot.^(2/(m-1));
            accum = accum + quot;
        end
        accum = 1./accum;
        U(i,:) = accum;
    end

    if (q ~= 0)
        %%
        %% step 4 : compute spatial membership function
        %%
        for i=1:c
            htmp = conv2(convkernel,reshape(U(i,:),size(Qin)));
            h(i,:) = reshape(htmp(convborder+1:end-convborder, ...
                convborder+1:end-convborder),1,length(Q));
        end

        tmp = ((U.^p).*(h.^q));
        hdenom = sum(tmp,1);
        Unew = tmp./repmat(hdenom,c,1);
        U=Unew;
        clear tmp;
    end

    dist = max(abs(vold-v));
end

%%
%% step 5: defuzzification and cluster reordering.
%%         clusters are reordered in increasing order of center
%%         values.
%%
[v,vidx] = sort(v);
U=U(vidx,:);
[val,idx] = max(U);

Qout = reshape(idx,size(Qin));
