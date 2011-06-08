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
%   
function [ims] = load_directory(dname,pattern)
    d = dir(strcat(dname,'/',pattern));
    h = waitbar(0,'Loading...');
    for i = 1:length(d)
        ims{i} = imread(strcat(dname,'/',d(i).name));
        waitbar(i/length(d));
    end
    close(h);
    
    h = waitbar(0,'Converting...');
    for i = 1:length(d)
        ims{i} = double(rgb2gray(ims{i}));
        waitbar(i/length(d));
    end
    
    close(h);
    
