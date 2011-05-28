fnamebase = '/Volumes/My Passport/worms/fall09_2/New4Matt ';
extn = '.png';
lo = 1;
hi = 1080;
thresh = 180;
% left, right, below, and two tight small lines to ensure work separation
% from clamp
premask = [[-1 -1 -1 285]; [-1 -1 490 -1]; [375 -1 -1 -1]; ...
    [120 265 401 409]; [180 265 375 380] ];
postmask = [[-1 185 -1 -1]];
output = 0;
outputdir = '/Users/matt/tmp/worms2';

[hr,hc] = process(fnamebase,extn,lo,hi,thresh,premask,postmask, ...
    output,outputdir);
