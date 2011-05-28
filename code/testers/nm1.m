fnamebase = '/Volumes/My Passport/worms/fall09/fall09 ';
extn = '.png';
lo = 1;
hi = 200;
thresh = 160;
premask = [[-1 -1 -1 150]; [-1 -1 550 -1]; [375 -1 -1 -1]];
postmask = [[-1 170 -1 -1]];

[hr,hc] = process(fnamebase,extn,lo,hi,thresh,premask,postmask);