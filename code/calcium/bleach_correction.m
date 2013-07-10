function [flat_ratio] =  bleach_correction (ratio)
% This function corrects for photo bleaching of ratiometric calcium
% indicators. The user inputs 4 point to fit the exponential trend, all
% points should be when calcium is not bound to probe. A double exponential
% is fit to these points and displayed so the user can evaluate goodness of
% fit. Then the fit is subtracted from the signal point by point and
% returned to the user. 
% Kathryn McCormick / 20 May 2013
% kat.mccormick@gmail.com
    figure; plot (ratio)
    [x,y] = ginput(4);
    f = fit(x,y,'exp2');
    hold on;
    plot(f,x,y);
    for i = 1:length(ratio);
        flat_ratio(i) = ratio(i) - f(i);
    end   
    figure; plot(ratio,'k');
    hold on; 
    plot(flat_ratio); 
end

    