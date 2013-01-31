%
% load a directory of images
%
% INPUT:
% 
%   dname   : The directory name containing the images.
%   pattern : Filename pattern for files in the directory, such as '*.png'
%
% OUTPUT:
%
%  ims      : cell array of N images that matched the pattern.  The images
%             are returned as double precision gray scale, NOT rgb.
%
% AUTHOR: matt sottile
%         mjsottile@me.com // msottile@uoregon.edu
% MODIFIER: kat mccormick
%           kat.mccormick@gmail.com
function [ims] = load_directory(dname,pattern)
    d = dir(strcat(dname,'/',pattern));
    h = waitbar(0,'Loading...');
    t = regexp(dname,'/','split')
    fnamebase = t{end}
    for i = 1:length(d)
        fname = sprintf('%s/%s_%d.tif',dname,fnamebase,i)
        ims{i} = imread(fname);
        waitbar(i/length(d));
    end
    close(h);
    
    h = waitbar(0,'Converting...');
    for i = 1:length(d)
        if (size(ims{i},3)>1)
            ims{i} = double(rgb2gray(ims{i}));
        else
            ims{i} = double(ims{i});
        end
        waitbar(i/length(d));
    end
    
    close(h);
