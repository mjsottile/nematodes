worms = [17:18];
datestr = '130527';


for i = 1:length(worms)
    wormnum = worms(i);
    fnamebase = sprintf('/Volumes/Internal 3/Calcium Images/%s/%sw%01d',...
        datestr,datestr,wormnum);
    disp(fnamebase);


    extn = '*.tif';

    [ims] = load_directory(fnamebase,extn);
    [ratio,yfp,cfp,refthresh,centx,centy, nangle] = calcium_process(ims,'circle',10,'bthresh',1, 'rig', 2);
    
    clear ims 

    savename = sprintf('%sw%02d.mat',datestr,wormnum);
    save(savename, 'ratio', 'yfp', 'cfp', 'nangle','refthresh','centx', 'centy');
end
%