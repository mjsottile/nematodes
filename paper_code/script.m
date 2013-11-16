% script to generate pictures for the paper
datapath = 'C:\Users\Matt\Dropbox\Shared w Matt\Clamp Images\120730w03 ';

% clip out visible edge on rhs of worm tail, and 
% channel border on rhs where nose touches
clip_region = [1,   190, 190, 300; ...
               155, 242, 325, 400];

preclip_region = [170, 172,  290, 325; ...
                  200, 450, 367, 638]; % chop off tail

% threshold : < 90

[hr, hc, cr, cc, hangle_seq] = process(datapath,'.png', 1, 200, 10, 90, preclip_region, clip_region, 1);

% run theta.m to regenerate angle data
