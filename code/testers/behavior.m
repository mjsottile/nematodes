worms = [16 17 18];
datestr = '100210';

for i = 1:length(worms)
    wormnum = worms(i);
    fnamebase = sprintf('/Users/lab/Documents/MATLAB/stills/%s/%sw%02d/%sw%02d ',...
        datestr,datestr,wormnum,datestr,wormnum);
    disp(fnamebase);

%fnamebase = '/Users/lab/Documents/MATLAB/stills/100210/100210w15/100210w15 ';
    extn = '.png';
    lo = 1;
    hi = 20;
    refframe = 1;
    thresh = 88;
    premask = [[-1 -1 -1 200]; [-1 -1 500 -1]; [375 -1 -1 -1]; [-1 150 -1 -1]; [218 229 363 369]; [200 208 382 390]; [237 250 345 347]]; 
    postmask = [[-1 150 -1 -1]; ];
    output = 1;
    outputdir = '/Users/lab/Documents/MatMovie';

    [hr,hc] = process(fnamebase,extn,lo,hi,refframe,thresh,premask,postmask, ...
        output,outputdir);

    savename = sprintf('%sw%02d.mat',datestr,wormnum);
    
    save(savename);
end
