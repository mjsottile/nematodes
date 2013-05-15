worms = [10:15];
datestr = '130125';


for i = 1:length(worms)
    wormnum = worms(i);
    fnamebase = sprintf('/Volumes/Internal 3/Calcium Images/%s/%sw%01d',...
        datestr,datestr,wormnum);
    disp(fnamebase);


    extn = '*.tif';

    [ims] = load_directory(fnamebase,extn);
    [ratio, yfp, cfp, nangle] = calcium_process(ims,'circle', 12);

    savename = sprintf('%sw%02d.mat',datestr,wormnum);
    
    clear ims 
    
    save(savename, 'ratio', 'yfp', 'cfp', 'nangle');
end
