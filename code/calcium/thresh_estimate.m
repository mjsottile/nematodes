function [t,bt] = thresh_estimate(im)

% compute thresholds for frame using kmeans for k=2, 10 iterations
[clusters, centers] = kmeans(im, 2, 10);

% not sure if safe to assume kmeans returns clusters in order of
% cluster center magnitude
[~,hi_cluster] = max(centers);
[~,lo_cluster] = min(centers);

bt_frac = 0.1; % choose value for bthresh that is the minimum
% from the dim cluster + 10% of the range.

lo_values = im(find(clusters==lo_cluster));
t = min(im(find(clusters==hi_cluster)));
bt = min(lo_values) + bt_frac * (max(lo_values)-min(lo_values));


%     numim = length(im);
%
%     bright_min = zeros(1,numim);
%     dim_ctr = zeros(1,numim);
%
%     for i=1:numim
%         [u,c] = kmeans(im{i},2,10);
%
%         masked = (u==2).*im{i};
%         masked(~masked) = inf;
%         bright_min(i) = min(masked(:));
%
%         dim_ctr(i) = c(1);
%     end
%
%     t = mean(bright_min) - std(bright_min);
%     bt = mean(dim_ctr) * 0.1;
%
