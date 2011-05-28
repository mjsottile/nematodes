fnamebase = '/Users/matt/Dropbox/100210w02/100210w02 ';
extn = '.png';

lo = 1;
hi = 9010;
refframe = 1;

thresh = 93;

premask = [[-1 -1 -1 140]; [-1 -1 475 -1]; [400 -1 -1 -1]; [-1 150 -1 -1];[150 257 308 311]; [150 251 331 340]; [170 195 397 405]; [200 227 364 367]];
postmask = [];

[hr,hc,cr,cc,hangle_seq] = process(fnamebase,extn,lo,hi,refframe,thresh,premask,postmask,1,'/tmp');