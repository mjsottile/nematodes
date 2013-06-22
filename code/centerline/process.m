%
% usage:
% 
%  process(filenamebase, extn, lo, hi, refframe, thresh, premask, 
%          postmask, output, outputdir)
%
%  INPUTS:
%
%  filenamebase : fully qualified base filename, minus the 4 digit
%                 count
%  extn         : file extension (eg: '.png').  need the dot.
%  lo           : file number to start with
%  hi           : file number to count up to
%  refframe     : reference frame for adapting threshold
%  thresh       : threshold value
%  premask      : M by 4 matrix of masking regions that will be
%                 zeroed out before large CC extraction
%  postmask     : M by 4 matrix of masking regions that will be
%                 zeroed out after large CC extraction
%  output       : 0 (default) = none, 1 (debug) = screen, 2 = files
%  outputdir    : if output == 2, this is the directory to dump to
%
% each row in the mask represents :
%  [ minrow maxrow mincol maxcol ]
%
% these define rectangles for zeroing out regions of the image.
% note that -1 in any entry corresponds to the appropriate extrema
% (such as 1 or end).
%
% also, this code assumes a maximum of 10000 frames, since frame numbers
% are hard coded to four digits.
%
% OUTPUTS:
% 
% hr         : vector of row positions of the head
% hc         : vector of column positions of the head
% cr         : approximate row of the clamp
% cc         : approximate column of the clamp
% hangle_seq : sequence of head angle values
%
%
% AUTHOR: matt sottile
%         mjsottile@me.com // msottile@uoregon.edu
%
function [hr,hc,cr,cc,hangle_seq] = ...
    processcline(filenamebase,extn,lo,hi,refframe,thresh,premask,...
                 postmask,output,outputdir)

%Checking to make sure there are enoght arguments coming in
doOutput = 0;  
if (nargin > 8)  
    doOutput = output;
    if (doOutput == 2)
        if (nargin < 10)
            error('Not enough arguments');
        end
    end
end

%creating matricies of appropriate size filled with zeros
hr = zeros((hi-lo)+1,1); 
hc = zeros((hi-lo)+1,1);

avg_time = 0; %initialization and assignment?

im = imread(sprintf('%s%04d%s',filenamebase,refframe,extn));
if (length(size(im))==2)
    im = double(im);
else
    im = double(rgb2gray(im));
end
refmean = mean(im(:));

pixcount = zeros(size(im));

for i=lo:hi
    tic;
    im = imread(sprintf('%s%04d%s',filenamebase,i,extn));
    if (length(size(im))==2)
        im = double(im);
    else
        im = double(rgb2gray(im));
    end

    cap = sprintf('Image %04d',i);
    
    % threshold
    b = double(im < (thresh*(mean(im(:))/refmean)));
    
    bdims = size(b);
    % premask
    for j = 1:size(premask,1)
        pm = premask(j,:);
        if (pm(1) == -1)
            pm(1) = 1;
        end
        if (pm(2) == -1)
            pm(2) = bdims(1);
        end
        if (pm(3) == -1)
            pm(3) = 1;
        end
        if (pm(4) == -1)
            pm(4) = bdims(2);
        end
        b(pm(1):pm(2),pm(3):pm(4)) = 0;
    end
        
    % label
    lbl = bwlabel(b);
    
    % geometric measures
    s = regionprops(lbl, 'area');
    
    % extract arrays
    areas = cat(1,s.Area);

    % sort by area
    [~,idx] = sort(areas);    
    
    % screen output
    if (doOutput == 1)
        subplot(2,2,1);imagesc(b);
    end
    
    % replace thresholded with connected component of biggest area
    % post-clipping
    b = (lbl == idx(end));
    
    % postmask
    for j = 1:size(postmask,1)
        pm = postmask(j,:);
        if (pm(1) == -1)
            pm(1) = 1;
        end
        if (pm(2) == -1)
            pm(2) = bdims(1);
        end
        if (pm(3) == -1)
            pm(3) = 1;
        end
        if (pm(4) == -1)
            pm(4) = bdims(2);
        end
        b(pm(1):pm(2),pm(3):pm(4)) = 0;
    end


    % screen output
    if (doOutput == 1)
        subplot(2,2,2);
        imagesc(b);
    end
    
    % remember body pixels
    pixcount = pixcount + b;
    
    % compute skeleton
    b = bwmorph(b,'skel',Inf);    

    % find head tip
    [headR, headC, b, cline_seq] = dijk(b);
    
    hangle_seq{(i-lo)+1} = cline_seq;
    
    % screen output
    if (doOutput == 1)
        subplot(2,2,3);
        imagesc(b);
        colormap(gray);
        hold on; 
        plot(headC,headR,'yx');
        plot(headC,headR,'b+');
        plot(headC,headR,'ro');
        hold off;
    end

    % screen output
    if (doOutput == 1)
        subplot(2,2,4);
    end
    
    % both screen and file output
    if (doOutput > 0)
        imagesc(im);
        colormap(gray);
        hold on; 
        plot(headC,headR,'yx');
        plot(headC,headR,'b+');
        plot(headC,headR,'ro');
        hold off;
        title(cap);
        drawnow;
    end
    
    hr(i-(lo-1)) = headR;
    hc(i-(lo-1)) = headC;
    
    % file output
    if (doOutput == 2)
        f = getframe;
        imwrite(f.cdata,sprintf('%s/out%03d.png',outputdir,i),'PNG');
    end
    
    % do some timing calculations to make a crude estimate of how much
    % longer this will take...
    delt = toc;
    
    avg_time = (avg_time * (i-lo)/((i-lo)+1)) + (delt / ((i-lo)+1));
    remaining_time = ((hi-i)*avg_time)/60;
    
    disp(sprintf('%f minutes remaining. (Image %03d)',remaining_time,i));
end

% pixcount has count for how many frames each pixel was on
pixcount_mask = double(pixcount == max(pixcount(:)));

[rs,cs,junk] = find(pixcount_mask);

% sort
[rs,idx] = sort(rs);

% take some pixel in the lowest row, and the corresponding column
cr = rs(end);
cc = cs(idx(end));
