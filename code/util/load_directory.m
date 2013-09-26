%
% Load a directory of images from user selected folder. 
%
% OUTPUT:
%
%  ims      : cell array of N images that matched the pattern.  The images
%             are returned as double precision gray scale, NOT rgb.
%
% AUTHORS: 
% 
% matt sottile
%         mjsottile@me.com // msottile@uoregon.edu
% kat mccormick
%         kat.mccormick@gmail.com

function [ims] = load_directory
    dname = uigetdir;
    d = dir(dname);
    names = {d.name};
    
    % Sometimes the directory contains empty items. Getting rid of
    % those here
    bytes = [d.bytes];
    names = names(bytes>0);
    
    h = waitbar(0,'Loading...');
    for i = 1:length(names)
        file_name = names{i};
        fname = sprintf('%s/%s',dname,file_name);
        ims{i} = imread(fname);
        waitbar(i/length(d));
    end
    close(h);
    h = waitbar(0,'Converting...');
    for i = 1:length(names)
        if (size(ims{i},3)>1)
            ims{i} = double(rgb2gray(ims{i}));
        else
            ims{i} = double(ims{i});
        end
        waitbar(i/length(names));
    end
    
    close(h);
    