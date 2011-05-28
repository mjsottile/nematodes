fnamebase = '/Volumes/My Passport/worms/fall09_3/New4Matt2 ';
extn = '.png';
lo = 1;
hi = 7000;
thresh = 165;
premask = [[-1 -1 -1 200]; [-1 -1 440 -1]; [390 -1 -1 -1]];
postmask = [[-1 175 -1 -1]];
output = 2;
outputdir = '/Users/matt/tmp/worms3';

[hr,hc] = process(fnamebase,extn,lo,hi,thresh,premask,postmask, ...
    output,outputdir);

save nm3.mat
