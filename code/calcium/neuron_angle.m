function [nangle] = neuron_angle(framesize, centx, centy, orientation) 
% Addendum to Calcium Imaging process code. This code identifies the angle
% that the neuron has from the center of the frame. Takes as input the size
% of the image, assuming that it will be split into YFP and CFP halves
% along the x axis, the centroid of the neuron in x and y, and the
% orientation of the prep. 
%   If orientation = 1, head is facing down
%   If orientation = 2, head is facing right

    ysize = framesize(1);
    xsize = framesize(2)/2;
    ymid = ysize/2;
    xmid = xsize/2;
    
    
    
    nangle = zeros(1, length(centx));
    
    if orientation ==1
        for i = 1:length(centx)
            xrel = centx(i)-xmid;
            yrel = centy(i) - ymid;
            if yrel > 0;
                nangle(i) = atan(xrel/yrel);
            elseif yrel < 0 && xrel > 0
                nangle(i) = atan(xrel/yrel)+ pi;
            elseif yrel < 0 && xrel < 0
                nangle(i) = atan(xrel/yrel)- pi;
            elseif yrel == 0 && xrel > 0 ;
                nangle(i) = pi/2;
            else
                nangle(i) = -pi/2;
            end
        end
    elseif orientation == 2
        for i = 1:length(centx)
            xrel = centx(i)-xmid;
            yrel = ymid - centy(i);
            nangle(i) = atan2(yrel, xrel);
        end
    else 
       disp('Error. Please provide valid orientation.')
    end
    
    
    nangle = nangle * 180/pi;
