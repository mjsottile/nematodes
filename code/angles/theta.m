%
% code to compute the angle from the head to the clamp
%

%% clamp row and column
cr=259;
cc=305;

%% assume hr and hc contain the row and column of the head
angle=atan((hc-cc)./(hr-cr));

mx=hc-cc;
my=hr-cr;

% some trig corrections based on where the head is relative to the
% clamp
corrected = ...
    angle.*double(sign(mx)<1).*double(sign(my)<1) + ...
    angle.*double(sign(mx)==1).*double(sign(my)<1) + ...
    (angle+pi).*double(sign(mx)<1).*(double(sign(my)==1)+double(sign(my)==0)) + ...
    (angle-pi).*double(sign(mx)==1).*(double(sign(my)==1)+double(sign(my)==0));

corrected = corrected .* double(sign(corrected)==1) + ...
            (corrected + 2*pi) .* double(sign(corrected)<1);
corrected=corrected-pi;

% convert to degrees from radians        
corrected=corrected.*(180.0/pi);
