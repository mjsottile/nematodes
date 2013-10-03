function [t,bt] = thresh_estimate(im)
% Gives a threshold and background threshold for a single frame. 
% kmeans input of 3 clusters, 10 iterations is typical. The user can input
% hard-coded fractions of the top and bottom cluster intensities to be
% returned in the t, bt variables. 


% compute thresholds for frame using kmeans for k=3, 10 iterations
[clusters, centers] = kmeans(im, 3, 10);

% not sure if safe to assume kmeans returns clusters in order of
% cluster center magnitude
[~,hi_cluster] = max(centers);
[~,lo_cluster] = min(centers);

% fraction of the dimmest cluster that will be used. typical values is 0.6%
bt_frac = 0.6;
% fraction of the brightest cluster that will be used. typical value is
% 0.4%
t_frac = 1.0; 


lo_values = im(find(clusters==lo_cluster));
hi_values = im(find(clusters==hi_cluster));

% threshold is maximum value from brightest cluster - t_frac % of the cluster range. 
t = max(hi_values) - t_frac * (max(hi_values) - min(hi_values));
% background thresh is the minimum from the dim cluster 
% + bt_frac% of the cluster range.
bt = min(lo_values) + bt_frac * (max(lo_values)-min(lo_values));
